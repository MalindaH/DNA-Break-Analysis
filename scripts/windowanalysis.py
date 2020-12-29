import sys
import os
import math
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from collections import defaultdict
from decimal import *



# Displays or updates a console progress bar: <0 = 'halt'; >=1 = 100%
def update_progress(progress, chrnum):
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress) 
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done!                                      \r\n"
    if progress >= 0 and progress < 1:
        status = "processing chromosome "+str(chrnum)+"..."
    block = int(round(barLength*progress))
    text = "\r  [{}] {:3.2f}% {}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


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
def calc_pval(chrnum, window_size, output_folder):
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
        x5 = 0
        x5_window_size = 0
        x11 = 0
        x11_window_size = 0
        for i in range(max(0, index-1), min(index+2, num_windows)):
            x3 += int(lines[i].split()[3])
            x3_window_size += window_size
        for i in range(max(0, index-2), min(index+3, num_windows)):
            x5 += int(lines[i].split()[3])
            x5_window_size += window_size
        for i in range(max(0, index-5), min(index+6, num_windows)):
            x11 += int(lines[i].split()[3])
            x11_window_size += window_size
        index += 1
 
        mu_window1 = x/window_size
        mu_window3 = x3/(x3_window_size)
        mu_window5 = x5/(x5_window_size)
        mu_window11 = x11/(x11_window_size)
        mu = mu_chr
        # topological domain average = 25kb
        if window_size <= 1000:
          # mu = expected number of reads in the window = max(mu_window1, mu_window3, mu_window5, mu_window11, mu_chr)
          mu = max(mu_window1, mu_window3, mu_window5, mu_window11, mu_chr)/mu_chr
        elif window_size <= 5000: 
          mu = max(mu_window1, mu_window3, mu_window5, mu_chr)/mu_chr
        elif window_size > 5000:
          mu = max(mu_window1, mu_window3, mu_chr)/mu_chr
        p_val = poisson_pval(x-1, mu) # upper tail, P(X >= x)
        # p_val = poisson.pmf(math.floor(mean_hits), mu)
        # print(p_val)
        output.write(linet.replace('\n','')+'\t--\t'+str(p_val)+'\t'+str(mu)+'\n')
    os.remove(temp_folder+"/outputtpy.txt")


# private method used in calc_pval()
# P(X > t); p-value = 1 - cdf(t, mu), this function avoids precision error from floating point number
def poisson_pval(t, mu): 
  result = Decimal(0.0)
  a = Decimal(mu)
  for k in range(t+1, t+50):
    b = Decimal(k)
    result = result + ((-a).exp())*((a)**b)/math.factorial(int(b))
    # result = result + math.exp(-mu)*((mu)**k)/math.factorial(k)
  return float(result)


# to find an appropriate window size: maximize variance of p-values
def find_window_size():
  print("Finding best window size...")
  # takes longer for smaller window sizes; more than 30000bp is too sparse
  sizes = [1000, 2000, 4000, 6000, 8000, 10000, 15000, 20000, 25000, 30000]
  # sizes = [500, 750, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 25000, 50000]

  df = pd.DataFrame(0.0, index=sizes, columns=['pval_variance', 'pval_numzeros'])

  window_size = -1

  for size in sizes:
    print('Trying window size: ',size,'bp...')
    i = 0
    for x in xs:
      update_progress(i/25, x)
      move_window(x, size)
      calc_pval(x, size, temp_folder)
      i += 1
    update_progress(1, 0)

    # collect variance and sum from output files
    find_var_sum(size, df, temp_folder) 

    # delete output files
    for x in xs:
      os.remove(temp_folder+'/chr'+str(x)+'_pval.txt')

    if df.at[size,'pval_variance'] < df.at[sizes[sizes.index(size)-1],'pval_variance']: # max variance is at sizes[sizes.index(size)-1]
      window_size = sizes[sizes.index(size)-1]
      break

  print(df)  
  if window_size == -1:
    window_size = df['pval_variance'].idxmax()
  print('Window size with highest variance of p-values is',window_size,'bp, use this size to analyze sensitive windows...')
  i = 0
  for x in xs:
    update_progress(i/25, x)
    move_window(x, window_size)
    calc_pval(x, window_size, output_folder)
    i += 1
  update_progress(10, -2)


# private method to be used in find_window_size()
def find_var_sum(size, df, output_folder):
    # collect variance and sum from output files
    zerosum = 0
    pvariance = []
    for x in xs:
      columns = defaultdict(list)
      fp = open(output_folder+'/chr'+str(x)+'_pval.txt', 'r')
      for line in fp:
          for (i,v) in enumerate(line.replace('\n','').split(sep='\t')):
              columns[i].append(v)
      
      pvals = np.array(columns[5], dtype=np.float32)
      pvariance.append(np.var(pvals))

      counter = 0
      for pval in columns[5]:
        if pval == '0.0': # need to use string here because floating point bad precision
          counter += 1
      zerosum += counter
    
    df.at[size, 'pval_variance'] = sum(pvariance)/len(pvariance)
    df.at[size, 'pval_numzeros'] = zerosum
    
    print('average variance of p-values: ',sum(pvariance)/len(pvariance))
    # print('sum: ',zerosum) 


#print("Usage: python windowanalysis.py temp output 10000 outputtest0825 blacklist_files $no_control")
print("\n-> Analyzing sensitive windows......")
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

if window_size == -1: # no window size provided
  find_window_size()
elif window_size > 0: # window size provided
  i = 0
  for x in xs:
    update_progress(i/25, x)
    move_window(x, window_size)
    calc_pval(x, window_size, output_folder)
    i += 1
  update_progress(1, 0)

  df = pd.DataFrame(0.0, index = range(1), columns=['pval_variance', 'pval_numzeros'])
  find_var_sum(window_size, df, output_folder)



