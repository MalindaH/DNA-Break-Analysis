import sys
import os
import numpy as np
from scipy.stats import hypergeom


def move_window(chrnum, direction):
  f = open(folder+"/chr"+str(chrnum)+"t"+direction+"hitsfiltered.txt", "r")
  bl = open(bl_folder+"/chr"+str(chrnum)+"_blacklist.txt", "r")
  outputt = open(folder+"/outputtpy.txt", "a")
  # keep track of which line has been read in bl
  bl_pos = 0
  
  x = 0
  y = x + window_size - 1
  for linebl in bl:
    bl_start = int(linebl.split()[1])
    bl_end = int(linebl.split()[2])
    if x <= bl_start and bl_start <= y:
      y = bl_end + window_size - 1 - bl_start + x
      bl_pos += len(linebl)
  count = 0
  for line in f:
    position = int(line.split()[4])
    if position <= y:
      count = count + 1
    elif position > y:
      outputt.write("chr"+str(chrnum)+"\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
      #print("chr1\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
      x = y + 1
      y = y + window_size
      bl.seek(bl_pos)
      for linebl in bl:
        bl_start = int(linebl.split()[1])
        bl_end = int(linebl.split()[2])
        if x <= bl_start and bl_start <= y:
          y = bl_end + window_size - 1 - bl_start + x
          bl_pos += len(linebl)
      count = 0
  outputt.close()
  f.close()

  f = open(folder+"/chr"+str(chrnum)+"c"+direction+"hitsfiltered.txt", "r")
  outputc = open(folder+"/outputcpy.txt", "a")

  count = 0
  x = 0
  y = x + window_size - 1
  for line in f:
    #print(line)
    position = int(line.split()[4])
    
    if position <= y:
      count = count + 1
    elif position > y:
      outputc.write("chr1\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
      #print("chr1\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
      x = y + 1
      y = y + window_size
      count = 0
  outputc.close()
  f.close()

  
# N = total number in population = number of reads in the given chromosome for
# both treated and non-treated samples
# k = total number with condition in population = number of reads in the given
# chromosome for the treated sample
# m = number in subset = number of reads in the given window for treated and
# non-treated samples
# x = number with condition in subset = number of reads in a given window for
# the treated sample
def calc_pval(chrnum, direction):
  # (optional) check outputtpy.txt and outputcpy.txt have same number of lines

  n = sum(1 for line in open(folder+"/chr"+str(chrnum)+"c"+direction+"hitsfiltered.txt"))
  k = sum(1 for line in open(folder+"/chr"+str(chrnum)+"t"+direction+"hitsfiltered.txt"))
  N = k + n
  if direction == '+':
    output = open(folder+'/chr'+str(chrnum)+'_fwdpval.txt', 'a+')
    with open(folder+"/outputtpy.txt") as ft, open(folder+"/outputcpy.txt") as fc:
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
  elif direction == '-':
    output = open(folder+'/chr'+str(chrnum)+'_revpval.txt', 'a+')
    with open(folder+"/outputtpy.txt") as ft, open(folder+"/outputcpy.txt") as fc:
      for linet, linec in zip(ft, fc):
        info = linet.split()[0]+"\t"+linet.split()[1]+"\t"+linet.split()[2]
        x = int(linet.split()[3])
        c = int(linec.split()[3])

        if x == 0 and c == 0:
          output.write(info+"\t0\t0\t1.0\n")
          continue

        m = x + c
        p_val = hypergeom.pmf(x, N, m, k)
        output.write(info+'\t'+str(x)+'\t'+str(c)+'\t'+str(p_val)+'\n')
  os.remove(folder+"/outputtpy.txt")
  os.remove(folder+"/outputcpy.txt")



# Benjamini-Hochberg correction, q-values are corrected p-values
def calc_qval(chrnum, direction):
  plist = []
  if direction == '+':
    with open(folder+'/chr'+str(chrnum)+'_fwdpval.txt', 'r') as pfile:
      # sorted(pfile, key = lambda line: line.split()[5])
      for line in pfile:
        plist.append(line.split()[5])
    qlist = p_adjust_bh(plist)
    output = open('output/'+sys.argv[3]+'chr'+str(chrnum)+'_fwd.txt', 'a+')
    index = 0
    with open(folder+'/chr'+str(chrnum)+'_fwdpval.txt', 'r') as pfile:
      for line in pfile:
        #info = linet.split()[0]+"\t"+linet.split()[1]+"\t"+linet.split()[2]+"\t"+linet.split()[3]+"\t"+linet.split()[4]+"\t"+linet.split()[5]
        output.write(line.rstrip("\n")+'\t'+str(qlist[index])+'\n')
        index += 1
    output.close()
    os.remove(folder+'/chr'+str(chrnum)+'_fwdpval.txt')
  elif direction == '-':
    with open(folder+'/chr'+str(chrnum)+'_revpval.txt', 'r') as pfile:
      # sorted(pfile, key = lambda line: line.split()[5])
      for line in pfile:
        plist.append(line.split()[5])
    qlist = p_adjust_bh(plist)
    output = open('output/'+sys.argv[3]+'chr'+str(chrnum)+'_rev.txt', 'a+')
    index = 0
    with open(folder+'/chr'+str(chrnum)+'_revpval.txt', 'r') as pfile:
      for line in pfile:
        #info = linet.split()[0]+"\t"+linet.split()[1]+"\t"+linet.split()[2]+"\t"+linet.split()[3]+"\t"+linet.split()[4]+"\t"+linet.split()[5]
        output.write(line.rstrip("\n")+'\t'+str(qlist[index])+'\n')
        index += 1
    output.close()
    os.remove(folder+'/chr'+str(chrnum)+'_revpval.txt')



# Benjamini-Hochberg p-value correction for multiple hypothesis testing 
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    #q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    q = np.minimum(1, steps * p[by_descend])
    return q[by_orig]


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


#print("Usage: python windowanalysis.py temp 48000 outputtest0825 blacklist_files")
print("Calculating p and q values......")
# main program:
folder = sys.argv[1]
window_size = int(sys.argv[2])
bl_folder = sys.argv[4]
# check if output file exits
if os.path.exists(sys.argv[3]+'_fwd.txt') or os.path.exists(sys.argv[3]+'_rev.txt'):
  print(sys.argv[3]+'_fwd.txt and/or '+sys.argv[3]+'_rev.txt already exists. Program stopped. Please change output filename.')
  exit(1)


'''
move_window(10, '+')
calc_pval(10, '+')
delete_p1(10, '+')
calc_qval(10, '+')
#move_window(2, '-')
#calc_pval(2, '-')
#calc_qval(2, '-')
'''

x = 1
while x <= 22:
  move_window(x, '+')
  calc_pval(x, '+')
  calc_qval(x, '+')
  move_window(x, '-')
  calc_pval(x, '-')
  calc_qval(x, '-')
  x += 1

move_window('X', '+')
calc_pval('X', '+')
calc_qval('X', '+')
move_window('X', '-')
calc_pval('X', '-')
calc_qval('X', '-')

move_window('Y', '+')
calc_pval('Y', '+')
calc_qval('Y', '+')
move_window('Y', '-')
calc_pval('Y', '-')
calc_qval('Y', '-')

move_window('M', '+')
calc_pval('M', '+')
calc_qval('M', '+')
move_window('M', '-')
calc_pval('M', '-')
calc_qval('M', '-')

  
# alpha = allowed false discovery rate (DFA), for calculating q-values
#alpha = 0.05

# test bh correction
'''
x = np.array([0.18562413450620835, 0.2112989052874284, 0.29661616404894264, 0.3606172669310967, 0.11795313342899175, 0.05321566053559549])
y = np.array([0.74879, 0.289455, 0.982681, 0.68646, 0.997208, 0.141006, 0.839462, 0.0144635, 0.899394, 0.22398, 0.0531287, 0.392512, 0.036908, 0.540398, 0.533397, 0.287477, 0.223107, 0.0943217, 0.863734, 0.118825])
print(p_adjust_bh(x))
print(p_adjust_bh(y))
'''

