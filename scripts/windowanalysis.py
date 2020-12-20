import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from scipy.stats import poisson
import math


def move_window(chrnum, window_size):
  bl = open(bl_folder+"/chr"+str(chrnum)+"_blacklist.txt", "r")
  f = open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt", "r")
  outputt = open(temp_folder+"/outputtpy.txt", "a+")
  # keep track of which line has been read in bl
  bl_pos = 0
  f_pos = 0
  x = 0
  y = x + window_size - 1
  # adjust window size wrt blacklist
  for linebl in bl:
    bl_start = int(linebl.split()[1])
    bl_end = int(linebl.split()[2])
    if x <= bl_start and bl_start <= y:
      y = bl_end + window_size - 1 - bl_start + x
      bl_pos += len(linebl)
    break
  count = 0
  for line in f:
    position = int(line.split()[4])
    # print("position: "+str(position)+"; x,y = "+str(x)+","+str(y))
    if position <= y:
      count += 1
      f_pos += len(line)
    else:
      outputt.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
      x = y + 1
      y = y + window_size
      count = 0
      f.seek(f_pos)
      bl.seek(bl_pos)
      linebl = bl.readline()
      if linebl:
        bl_start = int(linebl.split()[1])
        bl_end = int(linebl.split()[2])
        # print("x,y = "+str(x)+","+str(y))
        # print("blacklist: "+str(bl_start)+", "+str(bl_end))
        if x <= bl_start and bl_start <= y:
          y = bl_end + window_size - 1 - bl_start + x
          bl_pos += len(linebl)
  if x == 0:
    outputt.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
  outputt.close()
  f.close()
  bl.close()

  if no_control == '0': # with control:
    ref = open(temp_folder+"/outputtpy.txt", "r")
    f = open(temp_folder+"/chr"+str(chrnum)+"c_hitsfiltered.txt", "r")
    outputc = open(temp_folder+"/outputcpy.txt", "a+")
    bl_pos = 0
    f_pos = 0
    lineref = ref.readline()
    x = int(lineref.split()[1])
    y = int(lineref.split()[2])
    # adjust window size wrt blacklist
    for linebl in bl:
      bl_start = int(linebl.split()[1])
      bl_end = int(linebl.split()[2])
      if x <= bl_start and bl_start <= y:
        y = bl_end + window_size - 1 - bl_start + x
        bl_pos += len(linebl)
      break
    count = 0
    more = False
    for line in f:
      position = int(line.split()[4])
      if position <= y:
        count += 1
        f_pos += len(line)
      else:
        outputc.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
        lineref = ref.readline()
        more = False
        if lineref:
          more = True
          x = int(lineref.split()[1])
          y = int(lineref.split()[2])
          bl.seek(bl_pos)
          for linebl in bl:
            bl_start = int(linebl.split()[1])
            bl_end = int(linebl.split()[2])
            if x <= bl_start and bl_start <= y:
              y = bl_end + window_size - 1 - bl_start + x
              bl_pos += len(linebl)
            break
        count = 0
        f.seek(f_pos)
    if more:
      outputc.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t0\n")
      for lineref in ref:
        x = int(lineref.split()[1])
        y = int(lineref.split()[2])
        outputc.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t0\n")
    if x == 0:
      outputc.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
    outputc.close()
    f.close()
    ref.close()

  
# with control: hypergeometric p-value
# without control: poisson p-value
def calc_pval(chrnum, window_size):
  if no_control == '0': # with control:
    # N = total number in population = number of reads in the given chromosome for
    # both treated and non-treated samples
    # k = total number with condition in population = number of reads in the given
    # chromosome for the treated sample
    # m = number in subset = number of reads in the given window for treated and
    # non-treated samples
    # x = number with condition in subset = number of reads in a given window for
    # the treated sample
    n = sum(1 for line in open(temp_folder+"/chr"+str(chrnum)+"c_hitsfiltered.txt"))
    k = sum(1 for line in open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt"))
    N = k + n

    output = open(output_folder+'/chr'+str(chrnum)+'_pval.txt', 'a+')
    with open(temp_folder+"/outputtpy.txt") as ft, open(temp_folder+"/outputcpy.txt") as fc:
      for linet, linec in zip(ft, fc):
        info = linet.split()[0]+"\t"+linet.split()[1]+"\t"+linet.split()[2]
        x = int(linet.split()[3])
        c = int(linec.split()[3])
        if x == 0 and c == 0:
          output.write(info+"\t0\t0\t1.0\n")
          continue
        m = x + c
        #print('(N,k,m,x) = '+str(N)+','+str(k)+','+str(m)+','+str(x))
        p_val = hypergeom.pmf(x, N, m, k)
        #print(p_val)
        output.write(info+'\t'+str(x)+'\t'+str(c)+'\t'+str(p_val)+'\n')
    os.remove(temp_folder+"/outputtpy.txt")
    os.remove(temp_folder+"/outputcpy.txt")
  else:
    n = sum(1 for line in open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt"))
    num_windows = len(open(temp_folder+"/outputtpy.txt").readlines())
    mean_hits = n/num_windows
    # print("mean_hits = "+str(mean_hits))
    N = 2913022398 # Effective genome size of GRCh38, not used here
    chr_size = 0
    with open(temp_folder+"/outputtpy.txt", 'r') as ft:
      lines = ft.read().splitlines()
      chr_size = int(lines[-1].split()[2])
      mu_chr = n/chr_size

      output = open(output_folder+'/chr'+str(chrnum)+'_pval.txt', 'a+')
      index = 0 # linet = lines[index]
      for linet in lines:
        x = int(linet.split()[3])
        x3 = 0
        x3_window_size = 0
        x11 = 0
        x11_window_size = 0
        for i in range(max(0, index-1), min(index+2, num_windows)):
            x3 += int(lines[i].split()[3])
            x3_window_size += window_size
        for i in range(max(0, index-5), min(index+6, num_windows)):
            x11 += int(lines[i].split()[3])
            x11_window_size += window_size
        index += 1
 
        mu_window1 = x/window_size
        mu_window3 = x3/(x3_window_size)
        mu_window11 = x11/(x11_window_size)
        # mu = expected number of reads in the window = max(mu_window1, mu_window3, mu_window11, mu_chr)
        mu = max(mu_window1, mu_window3, mu_window11, mu_chr)
        # print('mu: ',mu_window1, mu_window3, mu_window11, mu_chr)
        p_val = poisson_pval(x-1, mu) # upper tail, P(X >= x)
        # p_val = poisson.pmf(math.floor(mean_hits), mu)
        # print(p_val)
        output.write(linet.replace('\n','')+'\t--\t'+str(p_val)+'\t'+str(mu)+'\n')
    os.remove(temp_folder+"/outputtpy.txt")


def poisson_pval(t, mu): # P(X > t); p-value = 1 - cdf(t, mu), this function avoids precision error from floating point number
  result = 0.0
  for k in range(t+1, 50):
    result = result + math.exp(-mu)*((mu)**k)/math.factorial(k)
  return result


# to find an appropriate window size
def find_window_size():
  sizes = [10, 20, 40, 60, 80, 100, 200, 500, 700, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000, 100000, 500000]
  
  indexes = np.array(xs).astype('str')
  out = pd.DataFrame(0, index=indexes, columns=['A', 'T', 'C', 'G'])

  for size in sizes:
    for x in xs:
      move_window(x, size)
      calc_pval(x, size)



#print("Usage: python windowanalysis.py temp output 10000 outputtest0825 blacklist_files $no_control")
print("Analyzing sensitive windows......")
# main program:
temp_folder = sys.argv[1]
output_folder = sys.argv[2]
window_size = int(sys.argv[3])
bl_folder = sys.argv[5]
no_control = sys.argv[6]



xs = list(range(1,23))
xs.append('X')
xs.append('Y')
xs.append('M')

for x in xs:
  move_window(x, window_size)
  calc_pval(x, window_size)


# # <- to find an appropriate window size -> ##
# sizes = [10, 20, 40, 60, 80, 100, 200, 500, 700, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000, 100000, 500000]
 
# indexes = np.array(xs).astype('str')
# out = pd.DataFrame(0, index=indexes, columns=['A', 'T', 'C', 'G'])

# for size in sizes:
#   for x in xs:
#     move_window(x, size)
#     calc_pval(x, size)

  

#   for x in xs:
#     os.remove(output_folder+'/chr'+str(x)+'_pval.txt')



