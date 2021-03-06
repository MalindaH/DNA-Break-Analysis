import sys
import os
import numpy as np
from scipy.stats import hypergeom
from scipy.stats import poisson
import math


def move_window(chrnum):
  bl = open(bl_folder+"/chr"+str(chrnum)+"_blacklist.txt", "r")
  f = open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt", "r")
  outputt = open(temp_folder+"/outputtpy.txt", "a+")
  # keep track of which line has been read in bl
  bl_pos = 0
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
    if position <= y:
      count = count + 1
    elif position > y:
      outputt.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
      x = y + 1
      y = y + window_size
      bl.seek(bl_pos)
      for linebl in bl:
        bl_start = int(linebl.split()[1])
        bl_end = int(linebl.split()[2])
        if x <= bl_start and bl_start <= y:
          y = bl_end + window_size - 1 - bl_start + x
          bl_pos += len(linebl)
        break
      count = 0
  if x == 0:
    print("write")
    outputt.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
  outputt.close()
  f.close()
  bl.close()

  if no_control == '0': # with control:
    ref = open(temp_folder+"/outputtpy.txt", "r")
    f = open(temp_folder+"/chr"+str(chrnum)+"c_hitsfiltered.txt", "r")
    outputc = open(temp_folder+"/outputcpy.txt", "a+")

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
        count = count + 1
      elif position > y:
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
    #print(str(x)+" and "+str(y))
    if more:
      #print('more')
      outputc.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t0\n")
      for lineref in ref:
        x = int(lineref.split()[1])
        y = int(lineref.split()[2])
        #print(str(x)+" and "+str(y))
        outputc.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t0\n")
    if x == 0:
      #print("write")
      outputc.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
    outputc.close()
    f.close()
    ref.close()

  
# with control: hypergeometric p-value
# without control: poisson p-value
def calc_pval(chrnum):
  if no_control == '0': # with control:
    # N = total number in population = number of reads in the given chromosome for
    # both treated and non-treated samples
    # k = total number with condition in population = number of reads in the given
    # chromosome for the treated sample
    # m = number in subset = number of reads in the given window for treated and
    # non-treated samples
    # x = number with condition in subset = number of reads in a given window for
    # the treated sample

    # (not needed) check outputtpy.txt and outputcpy.txt have same number of lines
    '''
    a = sum(1 for line in open(folder+"/outputtpy.txt"))
    b = sum(1 for line in open(folder+"/outputcpy.txt"))
    if a != b:
      print("Something went wrong in move_window(). Program stopped.")
      sys.exit(12)
    '''
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
    N = 2913022398 # Effective genome size of GRCh38

    output = open(output_folder+'/chr'+str(chrnum)+'_pval.txt', 'a+')
    with open(temp_folder+"/outputtpy.txt") as ft:
      for linet in ft:
        info = linet.split()[0]+"\t"+linet.split()[1]+"\t"+linet.split()[2]
        x = int(linet.split()[3])
        if x == 0:
          output.write(info+"\t0\t--\t1.0\n")
          continue
        # print('(N,n,x) = '+str(N)+','+str(n)+','+str(x))
        mu = (4000*x)/n
        # print(mu)
        p_val = poisson_pval(x-1, mu) # upper tail, P(X>=x)
        # p_val = poisson.pmf(mean_hits, mu)
        # print(p_val)
        output.write(info+'\t'+str(x)+'\t--\t'+str(p_val)+'\t'+str(mu)+'\n')
    os.remove(temp_folder+"/outputtpy.txt")


def poisson_pval(t, mu): # p-value = 1 - cdf(t, mu), this function avoids precision error from floating point number
  result = 0.0
  for k in range(t+1, 50):
    result = result + math.exp(-mu)*((mu)**k)/math.factorial(k)
  return result


# Benjamini-Hochberg correction, q-values are corrected p-values
# not working well, so not run
def calc_qval(chrnum):
  plist = []
  with open(output_folder+'/chr'+str(chrnum)+'_pval.txt', 'r') as pfile:
    # sorted(pfile, key = lambda line: line.split()[5])
    for line in pfile:
      plist.append(line.split()[5])
  qlist = p_adjust_bh(plist)
  output = open(output_folder+'/'+sys.argv[3]+'chr'+str(chrnum)+'.txt', 'a+')
  index = 0
  with open(output_folder+'/chr'+str(chrnum)+'_pval.txt', 'r') as pfile:
    for line in pfile:
      #info = linet.split()[0]+"\t"+linet.split()[1]+"\t"+linet.split()[2]+"\t"+linet.split()[3]+"\t"+linet.split()[4]+"\t"+linet.split()[5]
      output.write(line.rstrip("\n")+'\t'+str(qlist[index])+'\n')
      index += 1
  output.close()
  os.remove(output_folder+'/chr'+str(chrnum)+'_pval.txt')


# Benjamini-Hochberg p-value correction for multiple hypothesis testing 
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    #q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    q = np.minimum(1, steps * p[by_descend])
    return q[by_orig]

'''
def delete_p1(chrnum, direction):
  if direction == '+':
    output = open(folder+'/chr'+str(chrnum)+'_fwdpval2.txt', 'a+')
    with open(folder+'/chr'+str(chrnum)+'_fwdpval.txt', 'r') as pfile:
      for line in pfile:
        if float(line.split()[5]) != 1.0:
          output.write(line)
    output.close()
    #os.remove(folder+'/chr'+str(chrnum)+'_fwdpval.txt')
  else:
    output = open(folder+'/chr'+str(chrnum)+'_revpval2.txt', 'a+')
    with open(folder+'/chr'+str(chrnum)+'_revpval.txt', 'r') as pfile:
      for line in pfile:
        if float(line.split()[5]) != 1.0:
          output.write(line)
    output.close()
    os.remove(folder+'/chr'+str(chrnum)+'_revpval.txt')
'''

#print("Usage: python windowanalysis.py temp output 48000 outputtest0825 blacklist_files $no_control")
print("Analyzing sensitive windows......")
# main program:
temp_folder = sys.argv[1]
output_folder = sys.argv[2]
window_size = int(sys.argv[3])
bl_folder = sys.argv[5]
no_control = sys.argv[6]
# check if output file exits
# if os.path.exists(output_folder+'/'+sys.argv[4]+'.txt'):
#   print(output_folder+'/'+sys.argv[4]+'.txt already exists. Program stopped. Please change output filename.')
#   exit(11)



xs = list(range(1,23))
xs.append('X')
xs.append('Y')
xs.append('M')

for x in xs:
  move_window(x)
  calc_pval(x)
  # calc_qval(x)



  
# alpha = allowed false discovery rate (DFA), for calculating q-values
#alpha = 0.05

# test bh correction
'''
x = np.array([0.18562413450620835, 0.2112989052874284, 0.29661616404894264, 0.3606172669310967, 0.11795313342899175, 0.05321566053559549])
y = np.array([0.74879, 0.289455, 0.982681, 0.68646, 0.997208, 0.141006, 0.839462, 0.0144635, 0.899394, 0.22398, 0.0531287, 0.392512, 0.036908, 0.540398, 0.533397, 0.287477, 0.223107, 0.0943217, 0.863734, 0.118825])
print(p_adjust_bh(x))
print(p_adjust_bh(y))
'''

