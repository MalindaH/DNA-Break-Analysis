import sys
import os
from scipy.stats import hypergeom

print("Usage: python windowanalysis.py temp 48000 outputhhhhhh")

def move_window(chrnum, direction):
  f = open(folder+"/chr"+str(chrnum)+"t"+direction+"hitsnew.txt", "r")
  outputt = open(folder+"/outputtpy.txt", "a")

  count = 0
  x = 0
  y = x + window_size - 1
  for line in f:
    #print(line)
    position = int(line.split()[4])
    
    if position <= y:
      count = count + 1
    elif position > y:
      outputt.write("chr1\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
      #print("chr1\t"+str(x)+"\t"+str(y)+"\t"+str(count)+"\n")
      x = y + 1
      y = y + window_size
      count = 0

  outputt.close()
  f.close()

  f = open(folder+"/chr"+str(chrnum)+"c"+direction+"hitsnew.txt", "r")
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
# both treated and non-treated sample
# k = total number with condition in population = number of reads in the given
# chromosome for the treated sample
# m = number in subset = number of reads in the given window for treated and
# non-treated samples
# x = number with condition in subset = number of reads in a given window for
# the treated sample
def calc_pval(chrnum, direction):
  # check outputtpy.txt and outputcpy.txt have same number of lines

  n = sum(1 for line in open(folder+"/chr"+str(chrnum)+"c"+direction+"hitsnew.txt"))
  k = sum(1 for line in open(folder+"/chr"+str(chrnum)+"t"+direction+"hitsnew.txt"))
  N = k + n
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
  os.remove(folder+"/outputtpy.txt")
  os.remove(folder+"/outputcpy.txt")



folder = sys.argv[1]
window_size = int(sys.argv[2])
# check if output file exits
if os.path.exists(sys.argv[3]+'_fwd.txt') or os.path.exists(sys.argv[3]+'_rev.txt'):
  print(sys.argv[3]+'_fwd.txt and/or '+sys.argv[3]+'_rev.txt already exists. Program stopped. Please change output filename.')
  #exit(1)

output = open(sys.argv[3]+'_fwd.txt', 'a+')
print('Scanning chr1 + direction......')
#move_window(1,"+")
print('Calculating p-values......')
#calc_pval(1,"+")
print('chr1 + direction done :)')


output.close()

output = open(sys.argv[3]+'_rev.txt', 'a+')
print('Scanning chr1 - direction......')
move_window(1,"-")
print('Calculating p-values......')
calc_pval(1,"-")
print('chr1 - direction done :)')


output.close()

