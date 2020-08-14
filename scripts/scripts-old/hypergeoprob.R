#!/usr/bin/env Rscript
# Usage: Rscript hypergeoprob.R q m n k

# q = success in sample = args[1] = number of reads in a given 
# window for the treated sample
# c = number of reads in a given window for the non-treated sample
# m = sample size = args[2] = number of reads in the given 
# window for treated and non-treated samples = q + c
# n = failure in population = args[3] = number of reads in the 
# given chromosome for the non-treated sample
# k = sample in population = args[4] = number of reads in the 
# given chromosome for the treated sample
# N = n+k

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  print("Usage: Rscript hypergeoprob.R q m n k")
}

q <- strtoi(args[1])
m <- strtoi(args[2])
n <- strtoi(args[3])
k <- strtoi(args[4])
N <- n+k

# print(q)
# print(m)
# print(n)
# print(k)

#dhyper(12,20,85273,337059, log = TRUE)
#phyper(12,20,85273,337059, lower.tail = TRUE, log.p = FALSE)

fold.enrichment <- (q/k)/(m/N)
print(fold.enrichment)

hypergeometric <- dhyper(q, m, n, k, log = FALSE)
#hypergeometric <- phyper(q, m, n, k, lower.tail = FALSE, log.p = TRUE)
cat(hypergeometric)

q=12
m=20
n=85273
k=33705
x.range <- 0:min(k,m) 
## Compute the distribution of density P(X=x)
dens <- dhyper(x=x.range, m=m, n=n, k=k)
## Plot the distributon of hypergeometric densities
plot (x.range, dens, type="h", lwd=2, col="blue", main="Hypergeometric density", xlab="x = marked elements in the selection", ylab="density = P(X=x)", ylim=c(0, max(dens)*1.25))

