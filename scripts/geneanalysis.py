import sys
import os
import csv
import numpy as np
from decimal import *
import math
import operator as op
from functools import reduce


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
                print("ValueError?!")
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


# output files: column 1 - relative number of peaks, column 4 - p-value
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
                if lineref:
                    more = True
                    x = int(lineref.split()[2])
                    y = int(lineref.split()[3])
                else:
                    more = False
                    break
                count = 0
        #print(str(x)+" and "+str(y))
        if more:
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
                start = int(linet.split()[2])
                end = int(linet.split()[3])
                # info = linet.split()[0]+"\t"+linet.split()[1]+"\t"+linet.split()[2]
                x = int(round(float(linet.split()[0])*(end - start)))
                c = int(round(float(linec.split()[0])*(end - start)))
                if x == 0 and c == 0:
                    output.write("0\t0\t1.0"+'\t--\t'+linet)
                    continue
                m = x + c
                # print('(N,k,m,x) = '+str(N)+','+str(k)+','+str(m)+','+str(x))
                p_val = hypergeom_pval(N,k,m,x)
                # if str(p_val) == "0.0":
                #     print('(N,k,m,x) = '+str(N)+','+str(k)+','+str(m)+','+str(x))
                # print(p_val)
                output.write(linet.split()[0]+'\t'+linec.split()[0]+'\t'+str(p_val)+'\t--\t'+linet)
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


# private methods
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom

def hypergeom_pval(N,k,m,x):
    result = Decimal(ncr(m,x))*Decimal(ncr(N-m, k-x))/Decimal(ncr(N, k))
    return float(result)
    

# here p-value is not used for plotting / determining sensitive genes
# private method used in calc_pval()
# P(X > t); p-value = 1 - cdf(t, mu), this function avoids precision error from floating point number
def poisson_pval(t, mu): 
  result = Decimal(0.0)
  a = Decimal(mu)
  for k in range(math.floor(t+1), math.floor(t+1)):
    b = Decimal(k)
    result = result + ((-a).exp())*((a)**b)/math.factorial(int(b))
    # result = result + math.exp(-mu)*((mu)**k)/math.factorial(k)
  return float(result)


def process_refgene(xs):
    for x in xs:
        with open(anno_refgene, "r") as csv_file:
            outputf = open(anno_folder+"/chr"+str(x)+"_refgene.txt", "a+")
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                if row[2] == 'chr'+str(x): # or row[2].startswith('chr'+str(x)+'_'):
                    outputf.write(row[1]+'\t'+row[2]+'\t'+row[3]+'\t'+row[4]+'\t'+row[5]+'\t'+row[6]+'\t'+row[7]+'\t'+row[8]+'\n')
            outputf.close()
            lines = open(anno_folder+"/chr"+str(x)+"_refgene.txt", 'r').readlines()
            outputf2 = open(anno_folder+"/chr"+str(x)+"_refgene_sorted.txt", "a+")
            for line in sorted(lines, key=lambda line: int(line.split()[3])):
                outputf2.write(line)
            outputf2.close()
            os.remove(anno_folder+"/chr"+str(x)+"_refgene.txt")       


def alignto_refgene(chrnum):
    if no_control == '0': # with control:
        with open(anno_folder+"/chr"+str(chrnum)+"_refgene_sorted.txt", "r") as refg:
            outputf = open(output_folder+"/chr"+str(chrnum)+"_refgene_counts.txt", "a+")
            # name, chr, direction, TSS, TTS, windows 1-20, windows 21-180, windows 181-200 (1-20 before TSS, 181-200 after TTS)
            win_size1 = 250
            TSS_prev = 0
            TTS_prev = 0
            linefi_pos = 0
            for liner in refg:
                if liner.split()[2]=='+':
                    TSS = int(liner.split()[3])
                    TTS = int(liner.split()[4])
                    if TSS == TSS_prev and TTS == TTS_prev:
                        TSS_prev = TSS
                        TTS_prev = TTS
                        continue
                    TSS_prev = TSS
                    TTS_prev = TTS
                    temp_line = liner.split()[0]+'\t'+liner.split()[1]+'\t'+liner.split()[2]+'\t'+str(TSS)+'\t'+str(TTS)
                    win_size2 = abs(TTS-TSS)/160
                    start = max(TSS - 5000, 1)
                    end = start + win_size1
                    temp = ''
                    wrote1 = 0
                    wrote2 = 0
                    linefi_pos_updated = 0

                    while start < TTS + 5000:
                        # f = open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt", "r")
                        fi = open(temp_folder+"/chr"+str(x)+"_comparedhits_repeats_combined.txt", "r")
                        fi.seek(linefi_pos)
                        count = 0.0
                        for line in fi:
                            position = int(line.split()[1])
                            if position <= start and not linefi_pos_updated:
                                linefi_pos += len(line)
                            elif position <= end and position > start:
                                count += int(line.split()[2])
                                linefi_pos_updated = 1
                            elif position > end:
                                temp += str(count/(end-start))+','
                                break
                        fi.close()

                        start = end + 1
                        if start < TSS:
                            end = start + win_size1
                        elif start >= TSS and start < TTS:
                            end = start + win_size2
                            if not wrote1:
                                if temp == '':
                                    continue
                                wrote1 = 1
                                outputf.write(temp_line+'\t'+temp)
                                temp = ''
                        elif start >= TTS:
                            end = start + win_size1
                            if not wrote2:
                                wrote2 = 1
                                outputf.write('\t'+temp)
                                temp = ''
                    outputf.write('\t'+temp+'\n')
                else: # direction is '-'
                    TSS = int(liner.split()[4])
                    TTS = int(liner.split()[3])
                    if TSS == TSS_prev and TTS == TTS_prev:
                        TSS_prev = TSS
                        TTS_prev = TTS
                        continue
                    TSS_prev = TSS
                    TTS_prev = TTS
                    temp_line = liner.split()[0]+'\t'+liner.split()[1]+'\t'+liner.split()[2]+'\t'+str(TSS)+'\t'+str(TTS)
                    win_size2 = (TSS-TTS)/160
                    start = TSS + 5000
                    end = start - win_size1
                    temp = ''
                    wrote1 = 0
                    wrote2 = 0
                    linefi_pos_updated = 0

                    while end > TTS - 5000:
                        # f = open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt", "r")
                        fi = open(temp_folder+"/chr"+str(x)+"_comparedhits_repeats_combined.txt", "r")
                        fi.seek(linefi_pos)
                        count = 0.0
                        for line in fi:
                            position = int(line.split()[1])
                            if position <= end and not linefi_pos_updated:
                                linefi_pos += len(line)
                            if position <= start and position > end:
                                count += int(line.split()[2])
                                linefi_pos_updated = 1
                            elif position > start:
                                temp += str(count/(start-end))+','
                                break
                        fi.close()

                        start = end - 1
                        if start >= TSS:
                            end = start - win_size1
                        elif start >= TTS and start < TSS:
                            end = start - win_size2
                            if not wrote1:
                                if temp == '':
                                    continue
                                wrote1 = 1
                                outputf.write(temp_line+'\t'+temp)
                                temp = ''
                        elif start < TTS:
                            end = start - win_size1
                            if not wrote2:
                                wrote2 = 1
                                outputf.write('\t'+temp)
                                temp = ''
                        if end < 1:
                            break
                    outputf.write('\t'+temp+'\n')
            outputf.close()
        # os.remove(temp_folder+"/chr"+str(x)+"_comparedhits_repeats_combined.txt")
    else: # without control
        with open(anno_folder+"/chr"+str(chrnum)+"_refgene_sorted.txt", "r") as refg:
            outputf = open(output_folder+"/chr"+str(chrnum)+"_refgene_counts.txt", "a+")
            # name, chr, direction, TSS, TTS, windows 1-20, windows 21-180, windows 181-200 (1-20 before TSS, 181-200 after TTS)
            win_size1 = 250
            TSS_prev = 0
            TTS_prev = 0
            linef_pos = 0
            for liner in refg:
                if liner.split()[2]=='+':
                    TSS = int(liner.split()[3])
                    TTS = int(liner.split()[4])
                    if TSS == TSS_prev and TTS == TTS_prev:
                        TSS_prev = TSS
                        TTS_prev = TTS
                        continue
                    TSS_prev = TSS
                    TTS_prev = TTS
                    temp_line = liner.split()[0]+'\t'+liner.split()[1]+'\t'+liner.split()[2]+'\t'+str(TSS)+'\t'+str(TTS)
                    # outputf.write(liner.split()[0]+'\t'+liner.split()[1]+'\t'+liner.split()[2]+'\t'+str(TSS)+'\t'+str(TTS))
                    win_size2 = abs(TTS-TSS)/160
                    start = max(TSS - 5000, 1)
                    end = start + win_size1
                    temp = ''
                    wrote1 = 0
                    wrote2 = 0
                    linef_pos_updated = 0

                    while start < TTS + 5000:
                        f = open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt", "r")
                        f.seek(linef_pos)
                        count = 0.0
                        for line in f:
                            position = int(line.split()[4])
                            if position <= start and not linef_pos_updated:
                                linef_pos += len(line)
                            elif position <= end and position > start:
                                count += 1
                                linef_pos_updated = 1
                            elif position > end:
                                # outputf.write('\t'+str(count/(end-start)))
                                temp += str(count/(end-start))+','
                                break
                        f.close()

                        start = end + 1
                        if start < TSS:
                            end = start + win_size1
                        elif start >= TSS and start < TTS:
                            end = start + win_size2
                            if not wrote1:
                                if temp == '':
                                    continue
                                wrote1 = 1
                                outputf.write(temp_line+'\t'+temp)
                                temp = ''
                        elif start >= TTS:
                            end = start + win_size1
                            if not wrote2:
                                wrote2 = 1
                                outputf.write('\t'+temp)
                                temp = ''
                    outputf.write('\t'+temp+'\n')
                else: # direction is '-'
                    TSS = int(liner.split()[4])
                    TTS = int(liner.split()[3])
                    if TSS == TSS_prev and TTS == TTS_prev:
                        TSS_prev = TSS
                        TTS_prev = TTS
                        continue
                    TSS_prev = TSS
                    TTS_prev = TTS
                    temp_line = liner.split()[0]+'\t'+liner.split()[1]+'\t'+liner.split()[2]+'\t'+str(TSS)+'\t'+str(TTS)
                    win_size2 = (TSS-TTS)/160
                    start = TSS + 5000
                    end = start - win_size1
                    temp = ''
                    wrote1 = 0
                    wrote2 = 0
                    linef_pos_updated = 0

                    while end > TTS - 5000:
                        f = open(temp_folder+"/chr"+str(chrnum)+"t_hitsfiltered.txt", "r")
                        f.seek(linef_pos)
                        count = 0.0
                        for line in f:
                            position = int(line.split()[4])
                            if position <= end and not linef_pos_updated:
                                linef_pos += len(line)
                            if position <= start and position > end:
                                count += 1
                                linef_pos_updated = 1
                            elif position > start:
                                # outputf.write('\t'+str(count/(start-end)))
                                temp += str(count/(start-end))+','
                                break
                        f.close()

                        start = end - 1
                        if start >= TSS:
                            end = start - win_size1
                        elif start >= TTS and start < TSS:
                            end = start - win_size2
                            if not wrote1:
                                if temp == '':
                                    continue
                                wrote1 = 1
                                outputf.write(temp_line+'\t'+temp)
                                temp = ''
                        elif start < TTS:
                            end = start - win_size1
                            if not wrote2:
                                wrote2 = 1
                                outputf.write('\t'+temp)
                                temp = ''
                        if end < 1:
                            break
                    outputf.write('\t'+temp+'\n')
            outputf.close()
        

def rank_output():
    if os.path.exists('output/chr1_sensitive-cancer-genes.txt'):
        filenames = ['output/chr1_sensitive-cancer-genes.txt', 'output/chr2_sensitive-cancer-genes.txt', 'output/chr3_sensitive-cancer-genes.txt',
            'output/chr4_sensitive-cancer-genes.txt', 'output/chr5_sensitive-cancer-genes.txt', 'output/chr6_sensitive-cancer-genes.txt', 
            'output/chr7_sensitive-cancer-genes.txt', 'output/chr8_sensitive-cancer-genes.txt', 'output/chr9_sensitive-cancer-genes.txt', 
            'output/chr10_sensitive-cancer-genes.txt', 'output/chr11_sensitive-cancer-genes.txt', 'output/chr12_sensitive-cancer-genes.txt', 
            'output/chr13_sensitive-cancer-genes.txt', 'output/chr14_sensitive-cancer-genes.txt', 'output/chr15_sensitive-cancer-genes.txt', 
            'output/chr16_sensitive-cancer-genes.txt', 'output/chr17_sensitive-cancer-genes.txt', 'output/chr18_sensitive-cancer-genes.txt', 
            'output/chr19_sensitive-cancer-genes.txt', 'output/chr20_sensitive-cancer-genes.txt', 'output/chr21_sensitive-cancer-genes.txt',
            'output/chr22_sensitive-cancer-genes.txt', 'output/chrX_sensitive-cancer-genes.txt']
        with open('output/allchr_sensitive-cancer-genes-temp.txt', 'w+') as tfile:
            for fname in filenames:
                with open(fname) as infile:
                    tfile.write(infile.read())
        for fname in filenames:
            os.remove(fname)
        if no_control == '0': # with control:
            tempfile = open('output/allchr_sensitive-cancer-genes-temp.txt', 'r')
            tempfile2 = open('output/allchr_sensitive-cancer-genes-sorted-temp.txt', 'a+')
            for line in sorted(tempfile, key=lambda line: float(line.split()[0]), reverse=True):
                tempfile2.write(line)
            tempfile.close()
            os.remove('output/allchr_sensitive-cancer-genes-temp.txt')
            tempfile2.close()
            tempfile3 = open('output/allchr_sensitive-cancer-genes-sorted-temp.txt', 'r')
            tempfile4 = open('output/allchr_sensitive-cancer-genes-sorted-temp-bigpval.txt', 'a+')
            outfile = open('output/allchr_sensitive-cancer-genes-sorted.txt', 'a+')
            for line in tempfile3:
                if float(line.split()[2]) > 0.00001:
                    tempfile4.write(line)
                else:
                    outfile.write(line)
            tempfile4.close()
            tempfile5 = open('output/allchr_sensitive-cancer-genes-sorted-temp-bigpval.txt', 'r')
            outfile.write(tempfile5.read())
            outfile.close()
            os.remove('output/allchr_sensitive-cancer-genes-sorted-temp.txt')
            os.remove('output/allchr_sensitive-cancer-genes-sorted-temp-bigpval.txt')
        else:
            tempfile = open('output/allchr_sensitive-cancer-genes-temp.txt', 'r')
            outfile = open('output/allchr_sensitive-cancer-genes-sorted.txt', 'a+')
            for line in sorted(tempfile, key=lambda line: float(line.split()[0])):
                outfile.write(line)
            tempfile.close()
            os.remove('output/allchr_sensitive-cancer-genes-temp.txt')
            outfile.close()

    if os.path.exists('output/chr1_refgene_counts.txt'):
        filenames = ['output/chr1_refgene_counts.txt', 'output/chr2_refgene_counts.txt', 'output/chr3_refgene_counts.txt', 
            'output/chr4_refgene_counts.txt', 'output/chr5_refgene_counts.txt', 'output/chr6_refgene_counts.txt', 
            'output/chr7_refgene_counts.txt', 'output/chr8_refgene_counts.txt', 'output/chr9_refgene_counts.txt', 
            'output/chr10_refgene_counts.txt', 'output/chr11_refgene_counts.txt', 'output/chr12_refgene_counts.txt', 
            'output/chr13_refgene_counts.txt', 'output/chr14_refgene_counts.txt', 'output/chr15_refgene_counts.txt', 
            'output/chr16_refgene_counts.txt', 'output/chr17_refgene_counts.txt', 'output/chr18_refgene_counts.txt', 
            'output/chr19_refgene_counts.txt', 'output/chr20_refgene_counts.txt', 'output/chr21_refgene_counts.txt', 
            'output/chr22_refgene_counts.txt', 'output/chrX_refgene_counts.txt', 'output/chrY_refgene_counts.txt']
        with open('output/chrAll_refGene_counts.txt', 'w+') as tfile:
            for fname in filenames:
                with open(fname) as infile:
                    tfile.write(infile.read())
        for fname in filenames:
            os.remove(fname)


#print("Usage: python geneanalysis.py output annotation_files outputtest0827 ../genome-annotation/Census_all-Sep17-2020.csv temp $no_control ../genome-annotation/refGene.txt")
print("\n-> Analyzing sensitive genes......")
output_folder = sys.argv[1]
anno_folder = sys.argv[2]
output_name = sys.argv[3]
anno_cancer_file = sys.argv[4]
temp_folder = sys.argv[5]
no_control = sys.argv[6]
anno_refgene = sys.argv[7]

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



## <- find sensitive cancer genes -> ##

# only need to run this method once
# process_cancer_genes()

xs = list(range(1,23))
xs.append('X') # only chrX has cancer genes in cancer file

# i=0
# for x in xs:
#     update_progress(i/23, x)
#     alignto_cancer_genes(x)
#     calc_pval(x)
#     i += 1
# update_progress(1, 0)


print("Analyzing break density of genes...")

## <- find distribution wrt all genes -> ##
xs.append('Y')

# only need to run this method once
# process_refgene(xs)

i=0
for x in xs:
    update_progress(i/24, x)
    alignto_refgene(x)
    i += 1
update_progress(1, 0)

rank_output()
