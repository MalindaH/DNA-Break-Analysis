import sys
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

# Usage: python hypergeometrictest.py 422332 337059 20 12
# input (N, k, m, x)
# N = total number in population = sys.argv[1] = number of reads in the given 
# chromosome for both treated and non-treated sample
# k = total number with condition in population = sys.argv[2] = number of reads in the 
# given chromosome for the treated sample
# m = number in subset = sys.argv[3] = number of reads in the given window for treated 
# and non-treated samples
# x = number with condition in subset = sys.argv[4] = number of reads in a given 
# window for the treated sample

# cumulative probabilities:
#print ('p-value <= ' + sys.argv[4] + ': ' + str(stats.hypergeom.cdf(int(sys.argv[4]) ,int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))))
#print ('p-value >= ' + sys.argv[4] + ': ' + str(stats.hypergeom.sf(int(sys.argv[4]) - 1,int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))))
#print()

[N, k, m, x] = [int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])]
pval = stats.hypergeom.pmf(x, N, m, k)
print(pval)

# plot hypergeometric distribution: (use Jupyter Notebook on Windows to show)
#rv = stats.hypergeom(N, m, k)
#xrange = np.arange(0, m+1)
#pmf_alignments = rv.pmf(xrange)
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(xrange, pmf_alignments, 'bo')
#ax.vlines(xrange, 0, pmf_alignments, lw=2)
#ax.set_xlabel('# of alignments in the current window of the treated chromosome')
#ax.set_ylabel('hypergeometric PMF (p-value)')
#plt.show()