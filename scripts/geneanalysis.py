import sys
import os
import csv
import numpy as np

# align to genes from processed gtf file (gencode.v34.annotation.gtf)
def alignto_genes(chrnum):
    with open(anno_folder+"/chr"+str(chrnum)+"_annotation-genes.gtf", "r") as genes:
        with open(output_folder+"/"+output_name+"chr"+str(chrnum)+"_peaks.txt", "r") as peaks:
            outputf = open(output_folder+"/chr"+str(chrnum)+"_sensi-genes.txt", "a+")
            counter = 0
            p_pos = 0
            for lineg in genes:
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
                        p_pos += len(linep)
                        continue
                    elif p_start >= end:
                        if counter > 0:
                            outputf.write(str(counter)+"\t"+lineg)
                        # go back one line before break
                        peaks.seek(p_pos)
                        break
                if peaks.readline() == "":
                    break
                else:
                    peaks.seek(p_pos)
                    counter = 0
            outputf.close()


# edited output: column 1 - number of peaks within gene; column 2 - chromosome number; 
# column 3&4 - gene start & end positions; column 5 - gene name
def edit_output(chrnum):
    first = "gene_name \""
    last = "\"; "
    with open(output_folder+"/chr"+str(chrnum)+"_sensi-genes.txt", "r") as f:
        outputf = open(output_folder+"/chr"+str(chrnum)+"_sensitive-genes.txt", "a+")
        for s in f:
            try:
                start = s.rindex( first ) + len( first )
                end = s.index( last, start )
                outputf.write(s.split()[0]+"\t"+s.split()[1]+"\t"+s.split()[4]+"\t"+s.split()[5]+"\t"+s[start:end]+"\n")
                # print(s[start:end])
            except ValueError: # should not have an error here
                print("error?!")
        outputf.close()


# process csv file from Cancer gene census (Census_allThu Sep 17 03_10_33 2020.csv)
# output file: column 1 - chromosome number; column 2&3 - start&end positions; column 4 - gene symbol; 
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
        if(x == 22):
            y = 'X' # only chrX has cancer genes in this file
        if(x=='X'):
            break
        x += 1
        





#print("Usage: python geneanalysis.py output annotation_files outputtest0827 ../genome-annotation/Census_all-Sep17-2020.csv")
print("Analyzing sensitive genes......")
output_folder = sys.argv[1]
anno_folder = sys.argv[2]
output_name = sys.argv[3]
anno_cancer_file = sys.argv[4]

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


#process_cancer_genes()

