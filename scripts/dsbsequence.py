import sys
import os
import collections
import pandas as pd
import numpy as np
from collections import Counter

xs = list(range(1,23))
xs.append('X')
xs.append('Y')
xs.append('M')

def compare_tc(chrnum):
    # delete sequences slightly divergent owing to sequencing errors
    if no_control == '0':
        with open(temp_folder+'/chr'+str(chrnum)+'t_hitsfiltered.txt', 'r') as tfile, open(temp_folder+'/chr'+str(chrnum)+'c_hitsfiltered.txt', 'r') as cfile:
            outputf = open(temp_folder+"/chr"+str(chrnum)+"_comparedhits.txt", "a+")
            f_pos = 0
            for lc in cfile:
                pc = int(lc.split()[4])
                for lt in tfile:
                    pt = int(lt.split()[4])
                    if pt > pc + 12:
                        # go back one line before break
                        tfile.seek(f_pos)
                        break
                    elif abs(pt - pc) <= 12: # don't add this read
                        f_pos += len(lt)
                        break
                    else:
                        outputf.write(lt.split()[3]+'\t'+lt.split()[4]+'\t'+lt.split()[5]+'\t'+str(len(lt.split()[5]))+'\n')
                        f_pos += len(lt)
                        continue
            # print remaining lines
            for lt in tfile:
                outputf.write(lt.split()[3]+'\t'+lt.split()[4]+'\t'+lt.split()[5]+'\t'+str(len(lt.split()[5]))+'\n')
            outputf.close()
    else:
        with open(temp_folder+'/chr'+str(chrnum)+'t_hitsfiltered.txt', 'r') as tfile:
            outputf = open(temp_folder+"/chr"+str(chrnum)+"_comparedhits.txt", "a+")
            for lt in tfile:
                outputf.write(lt.split()[3]+'\t'+lt.split()[4]+'\t'+lt.split()[5]+'\t'+str(len(lt.split()[5]))+'\n')
            outputf.close()

    repetitive_reads(chrnum)
    slight_divergent(chrnum)


def repetitive_reads(chrnum):
    with open(temp_folder+"/chr"+str(chrnum)+"_comparedhits.txt", 'r') as tfile:
        outputf1 = open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats1.txt", "a+")
        counts = collections.Counter(l.strip() for l in tfile)
        for line, count in counts.most_common():
            outputf1.write(line+'\t'+str(count)+'\n')
        outputf1.close()
    # sort again by position
    outputf2 = open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats2.txt", "a+")
    lines = open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats1.txt", 'r').readlines()
    for line in sorted(lines, key=lambda line: int(line.split()[1])):
        outputf2.write(line)
    outputf2.close()


def slight_divergent(chrnum):
    with open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats2.txt", 'r') as tfile:
        outputf = open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats_combined.txt", "a+")
        prev_line = ''
        temp = []
        n = []
        for current_line in tfile: 
            if prev_line != '':
                if int(prev_line.split()[1])+5 > int(current_line.split()[1]):
                    temp.append(prev_line.split()[1])
                    n.append(int(prev_line.split()[4]))
                else:
                    if n: # if n not empty
                        temp.append(prev_line.split()[1])
                        n.append(int(prev_line.split()[4]))
                        index = n.index(max(n))
                        outputf.write(prev_line.split()[0]+'\t'+temp[index]+'\t'+str(sum(n))+'\n')
                        temp = []
                        n = []
                    else:
                        outputf.write(prev_line.split()[0]+'\t'+str(prev_line.split()[1])+'\t1\n')
            prev_line = current_line
        # deal with the last line
        if n: # if n not empty
            temp.append(prev_line.split()[1])
            n.append(int(prev_line.split()[4]))
            index = n.index(max(n))
            outputf.write(prev_line.split()[0]+'\t'+temp[index]+'\t'+str(sum(n))+'\n')
            temp = []
            n = []
        else:
            outputf.write(prev_line.split()[0]+'\t'+str(prev_line.split()[1])+'\t1\n')
        outputf.close()


# need to run this if hg19.fa contains sequence of all chromosomes, this method separates into each chromosome
# don't need to run if .fa files come in individual chromosomes
def generate_hg_files():
    if not os.path.isdir(hg_filepath):
        with open(hg_filepath, 'r') as hgfile:
            outputf = open(temp_folder+'/>chrtemp_hgfile.txt', 'a')
            copy = False
            for line in hgfile:
                if line.startswith('>chr'):
                    if line.strip().find('_') == -1:
                        copy = True
                        # print(line)
                        outputf.close()
                        outputf = open(temp_folder+'/'+line.strip()+'_hgfile.txt', 'a+')
                        outputf.write(line)
                    else:
                        # print(line)
                        copy = False
                elif copy:
                    outputf.write(line)
            outputf.close()
            os.remove(temp_folder+'/>chrtemp_hgfile.txt')
        hg_file_generated = True


# try to not be limited by 5-base window (ie. window_size = 2) as in generate_chr_sequences()
# window_size = # bp to count before and after the break
def generate_chr_sequences2(chrnum, window_size):
    if hg_file_generated:
        hg_file_individual = temp_folder+'/>chr'+str(chrnum)+'_hgfile.txt'
    else:
        hg_file_individual = hg_filepath+'/chr'+str(chrnum)+'.fa'
    with open(hg_file_individual, 'r') as hgfile:
        line = hgfile.readline().replace('\n', '')
        genome = hgfile.read().replace('\n', '')
        if not line.startswith('>chr'):
            line += genome
            genome = line
        with open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats_combined.txt", 'r') as hitsfile:
            indexes = np.array(range(-window_size, window_size+1)).astype('str')
            # print(indexes)
            out = pd.DataFrame(0, index=indexes, columns=['A', 'T', 'C', 'G'])
            for lineh in hitsfile:
                # pr1 is at position 0
                pr0 = int(lineh.split()[1]) - 1
                sequence = genome[pr0-window_size : pr0+window_size+1]
                add_count2(sequence.upper(), out, indexes)
    return out

# helper method for generate_chr_sequences2()
def add_count2(bases, df, indexes):
    if bases.find('N') != -1: # if bases sequence contains N
        return df
    index = 0
    for i in indexes:
        df.at[str(i), bases[index]] += 1
        index += 1
    return df       


def count_sample(output_csv, window_size):
    df = generate_chr_sequences2(1, window_size)
    for x in xs:
        if x == 1:
            df.to_csv(path_or_buf=output_csv+'-chr1.csv', mode='a')
        else:
            df_temp = generate_chr_sequences2(x, window_size)
            df_temp.to_csv(path_or_buf=output_csv+'-chr'+str(x)+'.csv', mode='a')
            df = df.add(df_temp)
    print(df)
    df.to_csv(path_or_buf=output_csv+'-chrAll.csv', mode='a')
    return df


def count_ref_chr2(chrnum):
    if hg_file_generated:
        hg_file_individual = temp_folder+'/>chr'+str(chrnum)+'_hgfile.txt'
    else:
        hg_file_individual = hg_filepath+'/chr'+str(chrnum)+'.fa'
    with open(hg_file_individual, 'r') as hgfile:
        # remove telomere sequences
        genome = ""
        fronttel_over = False
        endtel_met = False
        last_line_temp = ""
        after_endtel_temp = ""
        for linef in hgfile:
            if linef.startswith('>chr'):
                continue
            if not endtel_met and "TTAGGGTTAGGGTTAGGG" in linef.upper(): # end telomere
                last_line_temp = linef[:-1:].upper().replace("TTAGGG","zzzzzz")
                endtel_met = True
                after_endtel_temp += linef[:-1:]
                continue
            if endtel_met: # the end telomere was a false one
                after_endtel_temp += linef[:-1:]
                if not "TTAGGGTTAGGGTTAGGG" in linef.upper():
                    endtel_met = False
                    genome += after_endtel_temp
                    after_endtel_temp = ""
                    last_line_temp = ""
            if not fronttel_over and "CCCTAACCCTAACCCTAA" in linef.upper(): # front telomere
                genome += linef[:-1:].upper().replace("CCCTAA","zzzzzz")
                continue
            else:
                if "NNNNNN" not in linef:
                    fronttel_over = True
                genome += linef[:-1:] 
        genome += last_line_temp

        # remove blacklist regions of genome
        with open(bl_folder+"/chr"+str(chrnum)+"_blacklist.txt", "r") as bl:
            counter = Counter("")
            start_ref = 0
            for linebl in bl:
                start = int(linebl.split()[1])-1
                counter += Counter(genome[start_ref:start:])
                start_ref = int(linebl.split()[2])
            counter += Counter(genome[start_ref::])
            out = pd.DataFrame(0, index=['pos'], columns=['A', 'T', 'C', 'G'])
            out.at['pos', 'A'] += counter['A'] + counter['a']
            out.at['pos', 'T'] += counter['T'] + counter['t']
            out.at['pos', 'C'] += counter['C'] + counter['c']
            out.at['pos', 'G'] += counter['G'] + counter['g']
            return out


def count_hg(output_csv_hg):
    df = count_ref_chr2(1)
    for x in xs:
        if x == 1:
            print('chr'+str(x))
            print(df)
            df.to_csv(path_or_buf=output_csv_hg, mode='a+')
        else:
            print('chr'+str(x))
            df_temp = count_ref_chr2(x)
            df_temp.to_csv(path_or_buf=output_csv_hg, mode='a')
            df = df.add(df_temp)
            print(df)
    print(df)
    df.to_csv(path_or_buf=output_csv_hg, mode='a')
    return df



#print("Usage: python dsbsequence.py temp ../bowtie-files/hg19.fa output blacklist_files $no_control DSB-count-1126 2")
#print("Usage: python dsbsequence.py temp bowtie-files/GRCh38 output blacklist_files $no_control DSB-count-1126 10")
print("\n-> Analyzing break sequence bias......")
temp_folder = sys.argv[1]
hg_filepath = sys.argv[2]
output_folder = sys.argv[3]
bl_folder = sys.argv[4]
no_control = sys.argv[5]
hg_file_generated = False
output_csv = output_folder+'/'+sys.argv[6]
output_csv_hg = output_folder+'/'+sys.argv[6]+'-hg.csv'
window_size = int(sys.argv[7])


# for x in xs:
#     compare_tc(x)

# print('......generating counts......')
df = count_sample(output_csv, window_size)

# only need to do this step once (not needed for GRCh38)
# generate_hg_files()

df_hg = count_hg(output_csv_hg)

for x in xs:
    os.remove(temp_folder+"/chr"+str(x)+"_comparedhits.txt")
    os.remove(temp_folder+"/chr"+str(x)+"_comparedhits_repeats1.txt")
    os.remove(temp_folder+"/chr"+str(x)+"_comparedhits_repeats2.txt")
    os.remove(temp_folder+"/chr"+str(x)+"_comparedhits_repeats_combined.txt")


