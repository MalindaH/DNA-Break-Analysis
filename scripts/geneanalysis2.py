import sys
import os
import csv
import numpy as np
from decimal import *
import math
import operator as op


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

# precess gene expression file
def process_files():
    with open(anno_folder+"/annotation-genes.txt", "r") as f:
        allgenes = f.readlines()
    num_lines = sum(1 for line in open(anno_gene_exp)) - 1
    errfile = open(temp_folder+"/genes_not_found.txt", "a+")
    with open(anno_gene_exp) as expf:
        next(expf)
        counter = 1
        for geneline in sorted(expf, key=lambda geneline: geneline.split()[3], reverse = True):
            if counter <= 0.1*num_lines: # top 10% most expressed genes
                index = -1
                for i, l in enumerate(allgenes): # search for this gene's name in annotation-genes.txt file
                    if '-'+geneline.split()[0]+';' in l: 
                        index = i
                if index >= 0:
                    write_position(allgenes[index], "top", geneline)
                else:
                    errfile.write(geneline)
            elif counter >= 0.9*num_lines: # bottom 10% most expressed genes
                index = -1
                for i, l in enumerate(allgenes): # search for this gene's name in annotation-genes.txt file
                    if '-'+geneline.split()[0]+';' in l: 
                        index = i
                if index >= 0:
                    write_position(allgenes[index], "bottom", geneline)
                else:
                    errfile.write(geneline)
            counter += 1
    errfile.close()

# private method used in process_files(), retrieve chromosome number and position from a line 
# of annotation-genes.txt file, and write to file for each chromosome
def write_position(line, tb, geneline):
    lines = line.rstrip().split('\t')
    genelines = geneline.rstrip().split('\t')
    chrnum = -1
    chrn = lines[0]
    if chrn == "NC_000001.11":
            chrnum = 1
    elif chrn == "NC_000002.12":
            chrnum = 2
    elif chrn == "NC_000003.12":
            chrnum = 3
    elif chrn == "NC_000004.12":
            chrnum = 4
    elif chrn == "NC_000005.10":
            chrnum = 5
    elif chrn == "NC_000006.12":
            chrnum = 6
    elif chrn == "NC_000007.14":
            chrnum = 7
    elif chrn == "NC_000008.11":
            chrnum = 8
    elif chrn == "NC_000009.12":
            chrnum = 9
    elif chrn == "NC_000010.11":
            chrnum = 10
    elif chrn == "NC_000011.10":
            chrnum = 11
    elif chrn == "NC_000012.12":
            chrnum = 12
    elif chrn == "NC_000013.11":
            chrnum = 13
    elif chrn == "NC_000014.9":
            chrnum = 14
    elif chrn == "NC_000015.10":
            chrnum = 15
    elif chrn == "NC_000016.10":
            chrnum = 16
    elif chrn == "NC_000017.11":
            chrnum = 17
    elif chrn == "NC_000018.10":
            chrnum = 18
    elif chrn == "NC_000019.10":
            chrnum = 19
    elif chrn == "NC_000020.11":
            chrnum = 20
    elif chrn == "NC_000021.9":
            chrnum = 21
    elif chrn == "NC_000022.11":
            chrnum = 22
    elif chrn == "NC_000023.11":
            chrnum = "X"
    elif chrn == "NC_000024.10":
            chrnum = "Y"
    elif chrn == "NC_012920.1":
            chrnum = "M"
    # else:
    #     print("error retriving chrnum "+chrn)
    if chrnum != -1 and tb == "top":
        topfile = open(temp_folder+'/chr'+str(chrnum)+'_top_expressed.txt', 'a+')
        topfile.write(f'{genelines[0]}\tchr{chrnum}\t{lines[6]}\t{lines[3]}\t{lines[4]}\t{genelines[3]}\n')
        topfile.close()
    elif chrnum != -1 and tb == "bottom":
        btnfile = open(temp_folder+'/chr'+str(chrnum)+'_btm_expressed.txt', 'a+')
        # genename  chrnum    +-   startpos    endpos   
        btnfile.write(f'{genelines[0]}\tchr{chrnum}\t{lines[6]}\t{lines[3]}\t{lines[4]}\t{genelines[3]}\n')
        btnfile.close()    
    # return chrnum+'\t'+line.split()[3]+'\t'+line.split()[4]+'\n'


def alignto_geneexp(chrnum, tb):
    if os.path.isfile(temp_folder+"/chr"+str(chrnum)+"_"+tb+"_expressed.txt"):
        # if no_control == '0': # with control:
            with open(temp_folder+"/chr"+str(chrnum)+"_"+tb+"_expressed.txt", "r") as refg:
                outputf = open(output_folder+"/chr"+str(chrnum)+"_geneexp_counts_"+tb+".txt", "a+")
                # name, chr, direction, TSS, TTS, windows 1-20, windows 21-180, windows 181-200 (1-20 before TSS, 181-200 after TTS)
                # 5000bp upstream of TSS ~ 5000bp downstream of TTS
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
                        # print(f'size1 {win_size1} size2 {win_size2}')

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
                os.remove(temp_folder+"/chr"+str(chrnum)+"_"+tb+"_expressed.txt")
        # else: # without control

def process_output():
    # if os.path.exists('output/chr1_geneexp_counts_top.txt'):
        filenames = ['output/chr1_geneexp_counts_top.txt', 'output/chr2_geneexp_counts_top.txt', 'output/chr3_geneexp_counts_top.txt', 
            'output/chr4_geneexp_counts_top.txt', 'output/chr5_geneexp_counts_top.txt', 'output/chr6_geneexp_counts_top.txt', 
            'output/chr7_geneexp_counts_top.txt', 'output/chr8_geneexp_counts_top.txt', 'output/chr9_geneexp_counts_top.txt', 
            'output/chr10_geneexp_counts_top.txt', 'output/chr11_geneexp_counts_top.txt', 'output/chr12_geneexp_counts_top.txt', 
            'output/chr13_geneexp_counts_top.txt', 'output/chr14_geneexp_counts_top.txt', 'output/chr15_geneexp_counts_top.txt', 
            'output/chr16_geneexp_counts_top.txt', 'output/chr17_geneexp_counts_top.txt', 'output/chr18_geneexp_counts_top.txt', 
            'output/chr19_geneexp_counts_top.txt', 'output/chr20_geneexp_counts_top.txt', 'output/chr21_geneexp_counts_top.txt', 
            'output/chr22_geneexp_counts_top.txt', 'output/chrX_geneexp_counts_top.txt', 'output/chrY_geneexp_counts_top.txt', 
            'output/chrM_geneexp_counts_top.txt', ]
        filenames1 = []
        for f in filenames:
            if os.path.exists(f):
                filenames1.append(f)
        with open('output/chrAll_geneexp_counts_top.txt', 'w+') as tfile:
            for fname in filenames1:
                with open(fname) as infile:
                    tfile.write(infile.read())
        for fname in filenames1:
            os.remove(fname)
        
        filenames = ['output/chr1_geneexp_counts_btm.txt', 'output/chr2_geneexp_counts_btm.txt', 'output/chr3_geneexp_counts_btm.txt', 
            'output/chr4_geneexp_counts_btm.txt', 'output/chr5_geneexp_counts_btm.txt', 'output/chr6_geneexp_counts_btm.txt', 
            'output/chr7_geneexp_counts_btm.txt', 'output/chr8_geneexp_counts_btm.txt', 'output/chr9_geneexp_counts_btm.txt', 
            'output/chr10_geneexp_counts_btm.txt', 'output/chr11_geneexp_counts_btm.txt', 'output/chr12_geneexp_counts_btm.txt', 
            'output/chr13_geneexp_counts_btm.txt', 'output/chr14_geneexp_counts_btm.txt', 'output/chr15_geneexp_counts_btm.txt', 
            'output/chr16_geneexp_counts_btm.txt', 'output/chr17_geneexp_counts_btm.txt', 'output/chr18_geneexp_counts_btm.txt', 
            'output/chr19_geneexp_counts_btm.txt', 'output/chr20_geneexp_counts_btm.txt', 'output/chr21_geneexp_counts_btm.txt', 
            'output/chr22_geneexp_counts_btm.txt', 'output/chrX_geneexp_counts_btm.txt', 'output/chrY_geneexp_counts_btm.txt', 
            'output/chrM_geneexp_counts_btm.txt', ]
        filenames1 = []
        for f in filenames:
            if os.path.exists(f):
                filenames1.append(f)
        with open('output/chrAll_geneexp_counts_btm.txt', 'w+') as tfile:
            for fname in filenames1:
                with open(fname) as infile:
                    tfile.write(infile.read())
        for fname in filenames1:
            os.remove(fname)



#print("Usage: python geneanalysis2.py output annotation_files outputtestname temp $no_control ../genome-annotation/U2OS_gene_expression.txt")
print("\n-> Analyzing break density wrt gene expressions......")
# for a in sys.argv:
#     print(a)
output_folder = sys.argv[1]
anno_folder = sys.argv[2]
output_name = sys.argv[3]
temp_folder = sys.argv[4]
no_control = sys.argv[5]
anno_gene_exp = sys.argv[6]


xs = list(range(1,23))
xs.append('X')
xs.append('Y')
xs.append('M')

process_files()

i=0
for x in xs:
    update_progress(i/24, x)
    alignto_geneexp(x, "top")
    alignto_geneexp(x, "btm")
    i += 1
update_progress(1, 0)

process_output()

# if no_control == '0': # with control:
#     for x in xs:
#         os.remove(temp_folder+"/chr"+str(x)+"_comparedhits_repeats_combined.txt")

for x in xs:
        os.remove(temp_folder+"/chr"+str(x)+"_comparedhits_repeats_combined.txt")