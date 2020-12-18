import sys
import os
import csv
import numpy as np
from scipy.stats import hypergeom
from scipy.stats import poisson
import math

# align to genes from processed gtf file (gencode.v34.annotation.gtf)
def alignto_genes(chrnum):
    with open(anno_folder+"/chr"+str(chrnum)+"_annotation-genes.gtf", "r") as genes:
        with open(output_folder+"/"+output_name+"chr"+str(chrnum)+"_peaks.txt", "r") as peaks:
            outputf = open(output_folder+"/chr"+str(chrnum)+"_sensi-genes.txt", "a+")
            counter = 0
            p_pos = 0
            for lineg in genes:
                p_val_sum = 0
                start = int(lineg.split()[3])
                end = int(lineg.split()[4])
                #print("start = "+str(start)+", end = "+str(end))
                for linep in peaks:
                    p_start = int(linep.split()[1])
                    p_end = int(linep.split()[2])
                    #print("p_start = "+str(p_start)+", p_end = "+str(p_end))
                    if p_start < start and p_end <= start:
                        p_pos += len(linep)
                        continue
                    elif (p_end > start and p_end <= end) or (p_start >= start and p_start < end):
                        counter += 1
                        # p_val_sum += float(linep.split()[5])
                        p_val_sum += float(linep.split()[3])
                        p_pos += len(linep)
                        continue
                    elif p_start >= end:
                        if counter > 0:
                            p_val_avg = p_val_sum/float(counter)
                            p_val_log = -math.log(p_val_avg, 10)
                            outputf.write(str(counter)+"\t%.10f\t"%(p_val_log)+lineg)
                        # go back one line before break
                        peaks.seek(p_pos)
                        break
                if peaks.readline() == "":
                    break
                else:
                    peaks.seek(p_pos)
                    counter = 0
            outputf.close()


# edited output: column 1 - number of peaks within gene; column 2 - -log_10(average p-value); column 3 - chromosome number; 
# column 4&5 - gene start & end positions; column 6 - gene name
def edit_output(chrnum):
    first = "gene_name \""
    last = "\"; "
    with open(output_folder+"/chr"+str(chrnum)+"_sensi-genes.txt", "r") as f:
        outputf = open(output_folder+"/chr"+str(chrnum)+"_sensitive-genes.txt", "a+")
        for s in f:
            try:
                start = s.rindex( first ) + len( first )
                end = s.index( last, start )
                outputf.write(s.split()[0]+"\t"+s.split()[1]+"\t"+s.split()[2]+"\t"+s.split()[4]+"\t"+s.split()[5]+"\t"+s[start:end]+"\n")
                # print(s[start:end])
            except ValueError: # should not have an error here
                print("error?!")
        outputf.close()


# process csv file from Cancer gene census (Census_all-Sep17-2020.csv)
# output files: column 1 - chromosome number; column 2&3 - start&end positions; column 4 - gene symbol; 
# column 5 - gene name; column 6 - tumor type (somatic); column 7 - mutation
def process_cancer_genes():
    with open(anno_cancer_file, "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                genomic_location = row[3]
                colon = genomic_location.index(":")
                dash = genomic_location.index("-")
                if(dash > colon+1): # disregard those without position info
                    chr_num = genomic_location[0:colon]
                    start = genomic_location[colon+1:dash]
                    end = genomic_location[dash+1:]
                    # print(f'chrnum = {chr_num}, start = {start}, end = {end}')
                    with open(anno_folder+"/chr"+str(chr_num)+"_cancer-genes.txt", "a+") as outputf:
                        outputf.write(f'chr{chr_num}\t{start}\t{end}\t{row[0]}\t{row[1]}\t{row[9]}\t{row[15]}\n')
                line_count += 1
    x = 1
    y = ''
    while x <= 22 or y=='X':
        if y == 'X':
            x = 'X'
        with open(anno_folder+"/chr"+str(x)+"_cancer-genes.txt", "r") as unsorted:
            sorted = open(anno_folder+"/chr"+str(x)+"_cancer-genes-sorted.txt", "a+")
            l = []
            arr = []
            for row in unsorted:
                l.append(row)
                arr.append(row.split()[1])
            arr = np.asarray(arr).astype(int)
            by_ascend = arr.argsort()[::1]
            l = np.asarray(l)
            l_sorted = l[by_ascend]
            for content in l_sorted:
                sorted.write(content)
            sorted.close()
            os.remove(anno_folder+"/chr"+str(x)+"_cancer-genes.txt")
        if(x == 22):
            y = 'X' # only chrX has cancer genes in this file
        if(x=='X'):
            break
        x += 1



# align peaks to csv file from Cancer gene census (Census_all-Sep17-2020.csv)
# output files: column 1 - number of peaks, column 2 - -log_10(average p-value)
def alignto_cancer_genes_old(chrnum): # deprecated
    with open(anno_folder+"/chr"+str(chrnum)+"_cancer-genes-sorted.txt", "r") as genes:
        with open(output_folder+"/"+output_name+"chr"+str(chrnum)+"_peaks.txt", "r") as peaks:
            outputf = open(output_folder+"/chr"+str(chrnum)+"_sensitive-cancer-genes.txt", "a+")
            counter = 0
            p_pos = 0
            for lineg in genes:
                p_val_sum = 0
                start = int(lineg.split()[1])
                end = int(lineg.split()[2])
                #print("start = "+str(start)+", end = "+str(end))
                for linep in peaks:
                    p_start = int(linep.split()[1])
                    p_end = int(linep.split()[2])
                    #print("p_start = "+str(p_start)+", p_end = "+str(p_end))
                    if p_start < start and p_end <= start:
                        p_pos += len(linep)
                        continue
                    elif (p_end > start and p_end <= end) or (p_start >= start and p_start < end):
                        counter += 1
                        # p_val_sum += float(linep.split()[5])
                        p_val_sum += float(linep.split()[3])
                        p_pos += len(linep)
                        continue
                    elif p_start >= end:
                        if counter > 0:
                            p_val_avg = p_val_sum/float(counter)
                            p_val_log = -math.log(p_val_avg, 10)
                            outputf.write(str(counter)+"\t%.10f\t"%(p_val_log)+lineg)
                        # go back one line before break
                        peaks.seek(p_pos)
                        break
                if peaks.readline() == "":
                    break
                else:
                    peaks.seek(p_pos)
                    counter = 0
            outputf.close()


# move windows of csv file from Cancer gene census (Census_all-Sep17-2020.csv)
def alignto_cancer_genes(chrnum):
    with open(anno_folder+"/chr"+str(chrnum)+"_cancer-genes-sorted.txt", "r") as genes:
        f = open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt", "r")
        outputf = open(temp_folder+"/outputtpy.txt", "a+")
        # keep track of which line has been read in genes
        genes_pos = 0
        for lineg in genes:
            p_val_sum = 0
            start = int(lineg.split()[1])
            end = int(lineg.split()[2])
            #print("start = "+str(start)+", end = "+str(end))
            count = 0.0
            for line in f:
                position = int(line.split()[4])
                if position <= end:
                    count = count + 1
                elif position > end:
                    outputf.write(str(count/(end-start))+'\t'+lineg)
                    break
            if start == 0:
                outputf.write(str(count/(end-start))+'\t'+lineg)
        outputf.close()
        f.close()
    if no_control == '0': # with control:
        ref = open(temp_folder+"/outputtpy.txt", "r")
        f = open(temp_folder+"/chr"+str(chrnum)+"c_hitsfiltered.txt", "r")
        outputc = open(temp_folder+"/outputcpy.txt", "a+")

        lineref = ref.readline()
        x = int(lineref.split()[2])
        y = int(lineref.split()[3])
        
        count = 0
        more = False
        for line in f:
            position = int(line.split()[4])
            if position <= y:
                count = count + 1
            elif position > y:
                outputc.write(str(count/(end-start))+'\t'+lineg)
                lineref = ref.readline()
                more = False
                if lineref:
                    more = True
                    x = int(lineref.split()[2])
                    y = int(lineref.split()[3])
                count = 0
        #print(str(x)+" and "+str(y))
        if more:
            #print('more')
            outputc.write('0\t'+lineg)
            for lineref in ref:
                x = int(lineref.split()[2])
                y = int(lineref.split()[3])
                #print(str(x)+" and "+str(y))
                outputc.write('0\t'+lineg)
        if x == 0:
            #print("write")
            outputc.write(str(count/(end-start))+'\t'+lineg)
        outputc.close()
        f.close()
        ref.close()


# output files: column 1 - number of peaks, column 2 - -log_10(average p-value)
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
        n = sum(1 for line in open(temp_folder+"/chr"+str(chrnum)+"c_hitsfiltered.txt"))
        k = sum(1 for line in open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt"))
        N = k + n
        output = open(output_folder+"/chr"+str(chrnum)+"_sensitive-cancer-genes.txt", "a+")
        with open(temp_folder+"/outputtpy.txt") as ft, open(temp_folder+"/outputcpy.txt") as fc:
            for linet, linec in zip(ft, fc):
                # info = linet.split()[0]+"\t"+linet.split()[1]+"\t"+linet.split()[2]
                x = float(linet.split()[0])
                c = float(linec.split()[0])
                if x == 0 and c == 0:
                    output.write(linet+"\t0\t0\t1.0\n")
                    continue
                m = x + c
                #print('(N,k,m,x) = '+str(N)+','+str(k)+','+str(m)+','+str(x))
                p_val = hypergeom.pmf(x, N, m, k)
                #print(p_val)
                output.write(linet+'\t'+str(x)+'\t'+str(c)+'\t'+str(p_val)+'\n')
        output.close()
        os.remove(temp_folder+"/outputtpy.txt")
        os.remove(temp_folder+"/outputcpy.txt")
    else:
        n = sum(1 for line in open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt"))
        num_windows = len(open(temp_folder+"/outputtpy.txt").readlines())
        mean_hits = n/num_windows
        # print("mean_hits = "+str(mean_hits))
        N = 2913022398 # Effective genome size of GRCh38

        output = open(output_folder+"/chr"+str(chrnum)+"_sensitive-cancer-genes.txt", "a+")
        with open(temp_folder+"/outputtpy.txt") as ft:
            for linet in ft:
                # info = linet.split()[0]+"\t"+linet.split()[1]+"\t"+linet.split()[2]
                x = float(linet.split()[0])
                if x == 0.0:
                    output.write("0\t--\t1.0\t"+linet)
                    continue
                # print('(N,n,x) = '+str(N)+','+str(n)+','+str(x))
                mu = (4000*x)/n
                # print(mu)
                p_val = poisson_pval(x-1, mu) # upper tail, P(X>=x)
                # p_val = poisson.pmf(mean_hits, mu)
                # print(p_val)
                output.write(str(x)+'\t--\t'+str(p_val)+'\t'+str(mu)+'\t'+linet)
        output.close()
        os.remove(temp_folder+"/outputtpy.txt")


def poisson_pval(t, mu): # p-value = 1 - cdf(t, mu), this function avoids precision error from floating point number
  result = 0.0
  for k in range(math.floor(t+1), 50):
    result = result + math.exp(-mu)*((mu)**k)/math.factorial(k)
  return result



#print("Usage: python geneanalysis.py output annotation_files outputtest0827 ../genome-annotation/Census_all-Sep17-2020.csv temp $no_control")
print("Analyzing sensitive genes......")
output_folder = sys.argv[1]
anno_folder = sys.argv[2]
output_name = sys.argv[3]
anno_cancer_file = sys.argv[4]
# window_size = sys.argv[5]
temp_folder = sys.argv[5]
# bl_folder = sys.argv[6]
no_control = sys.argv[6]

# x = 1
# while x <= 22:
#     alignto_genes(x)
#     x += 1

# alignto_genes("X")
# alignto_genes("Y")
# alignto_genes("M")

# x = 1
# while x <= 22:
#     edit_output(x)
#     x += 1

# edit_output("X")
# edit_output("Y")
# edit_output("M")

# only need to run once
# process_cancer_genes()

# alignto_cancer_genes(1)
# calc_pval(1)

x = 1
while x <= 22:
    alignto_cancer_genes(x)
    calc_pval(x)
    x += 1
alignto_cancer_genes("X") # only chrX has cancer genes in cancer file
calc_pval("X")