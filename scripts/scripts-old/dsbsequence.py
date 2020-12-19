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


# (deprecated)
# output 4x5 array (dataframe): out[0,:] is count of 4 bases (A,T,C,G) for position -2
# out[1,:] is count of 4 bases (A,T,C,G) for position -1
# out[2,:] is count of 4 bases (A,T,C,G) for position 0 (break position)
# out[3,:] is count of 4 bases (A,T,C,G) for position +1
# out[4,:] is count of 4 bases (A,T,C,G) for position +2
# out[:,0] is count of A, out[:,1] is count of T, out[:,2] is count of C, out[:,3] is count of G
def generate_chr_sequences(chrnum):
    with open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats_combined.txt", 'r') as hitsfile:
        if hg_file_generated:
            hg_file_individual = temp_folder+'/>chr'+str(chrnum)+'_hgfile.txt'
        else:
            hg_file_individual = hg_filepath+'/chr'+str(chrnum)+'.fa'
        with open(hg_file_individual, 'r') as hgfile:
            pos = 0
            f_pos = 0
            bug_count = 0
            out = pd.DataFrame(0, index=['pos -2', 'pos -1', 'pos 0', 'pos +1', 'pos +2'], columns=['A', 'T', 'C', 'G'])
            for lineh in hitsfile:
                # pr1 is at position +1
                pr1 = int(lineh.split()[1]) - 1
                fivebases = ''
                for lineg in hgfile:
                    if lineg.startswith('>chr'):
                        continue
                    if pos + len(lineg)-1 < pr1-3:
                        pos += len(lineg)-1 # -1 is for the '\n' character
                        f_pos += len(lineg)
                        continue
                    elif pos < pr1-3 and pos + len(lineg)-1 >= pr1-3:
                        if pr1+1 <= pos + len(lineg)-1:
                            fivebases += lineg.strip()[pr1-4-pos : pr1+1-pos]
                            hgfile.seek(f_pos)
                            break
                        else:
                            fivebases += lineg.strip()[pr1-4-pos : pos + len(lineg)]
                            pos += len(lineg)-1 # -1 is for the '\n' character
                            f_pos += len(lineg)
                            continue
                    else:
                        fivebases += lineg.strip()[ : pr1+1-pos]
                        hgfile.seek(f_pos)
                        break
                # print("five bases: "+fivebases+' at position: '+str(pr1+1))
                if len(fivebases) == 5:
                    add_count(fivebases.upper(), out)
                    # print(out)
                else:
                    bug_count += 1
                    # print("!!!!!!five bases: "+fivebases+' at position: '+str(pr1+1))
            if bug_count > 0:
                print("total bug count: "+str(bug_count))
    # print('chr'+str(chrnum)+':')
    # print(out)
    # print()
    return out


# (deprecated)
# helper method for generate_chr_sequences()
def add_count(fivebases, df):
    if fivebases.find('N') != -1:
        return df
    df.at['pos -2', fivebases[0]] += 1
    df.at['pos -1', fivebases[1]] += 1
    df.at['pos 0', fivebases[2]] += 1
    df.at['pos +1', fivebases[3]] += 1
    df.at['pos +2', fivebases[4]] += 1
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

def count_sample_old(output_csv):
    df1 = generate_chr_sequences(1)
    df1.to_csv(path_or_buf=output_csv, mode='a')
    df2 = generate_chr_sequences(2)
    df2.to_csv(path_or_buf=output_csv, mode='a')
    df3 = generate_chr_sequences(3)
    df3.to_csv(path_or_buf=output_csv, mode='a')
    df4 = generate_chr_sequences(4)
    df4.to_csv(path_or_buf=output_csv, mode='a')
    df5 = generate_chr_sequences(5)
    df5.to_csv(path_or_buf=output_csv, mode='a')
    df6 = generate_chr_sequences(6)
    df6.to_csv(path_or_buf=output_csv, mode='a')
    df7 = generate_chr_sequences(7)
    df7.to_csv(path_or_buf=output_csv, mode='a')
    df8 = generate_chr_sequences(8)
    df8.to_csv(path_or_buf=output_csv, mode='a')
    df9 = generate_chr_sequences(9)
    df9.to_csv(path_or_buf=output_csv, mode='a')
    df10 = generate_chr_sequences(10)
    df10.to_csv(path_or_buf=output_csv, mode='a')
    df11 = generate_chr_sequences(11)
    df11.to_csv(path_or_buf=output_csv, mode='a')
    df12 = generate_chr_sequences(12)
    df12.to_csv(path_or_buf=output_csv, mode='a')
    df13 = generate_chr_sequences(13)
    df13.to_csv(path_or_buf=output_csv, mode='a')
    df14 = generate_chr_sequences(14)
    df14.to_csv(path_or_buf=output_csv, mode='a')
    df15 = generate_chr_sequences(15)
    df15.to_csv(path_or_buf=output_csv, mode='a')
    df16 = generate_chr_sequences(16)
    df16.to_csv(path_or_buf=output_csv, mode='a')
    df17 = generate_chr_sequences(17)
    df17.to_csv(path_or_buf=output_csv, mode='a')
    df18 = generate_chr_sequences(18)
    df18.to_csv(path_or_buf=output_csv, mode='a')
    df19 = generate_chr_sequences(19)
    df19.to_csv(path_or_buf=output_csv, mode='a')
    df20 = generate_chr_sequences(20)
    df20.to_csv(path_or_buf=output_csv, mode='a')
    df21 = generate_chr_sequences(21)
    df21.to_csv(path_or_buf=output_csv, mode='a')
    df22 = generate_chr_sequences(22)
    df22.to_csv(path_or_buf=output_csv, mode='a')
    dfX = generate_chr_sequences('X')
    dfX.to_csv(path_or_buf=output_csv, mode='a')
    dfY = generate_chr_sequences('Y')
    dfY.to_csv(path_or_buf=output_csv, mode='a')
    dfM = generate_chr_sequences('M')
    dfM.to_csv(path_or_buf=output_csv, mode='a')

    df = df1.copy()
    df.add(df2)
    df.add(df3)
    df.add(df4)
    df.add(df5)
    df.add(df6)
    df.add(df7)
    df.add(df8)
    df.add(df9)
    df.add(df10)
    df.add(df11)
    df.add(df12)
    df.add(df13)
    df.add(df14)
    df.add(df15)
    df.add(df16)
    df.add(df17)
    df.add(df18)
    df.add(df19)
    df.add(df20)
    df.add(df21)
    df.add(df22)
    print(df)
    df.to_csv(path_or_buf=output_csv, mode='a')
    return df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19,df20,df21,df22,dfX,dfY,dfM,df

# (deprecated)
def count_ref_chr(chrnum):
    if hg_file_generated:
        hg_file_individual = temp_folder+'/>chr'+str(chrnum)+'_hgfile.txt'
    else:
        hg_file_individual = hg_filepath+'/chr'+str(chrnum)+'.fa'
    with open(hg_file_individual, 'r') as hgfile:
        # remove blacklist regions for hg file
        with open(bl_folder+"/chr"+str(chrnum)+"_blacklist.txt", "r") as bl:
            outputf = open(temp_folder+"/chr"+str(chrnum)+"_hgfiltered.txt", "a+")
            f_pos = 0
            position = 0 # position of start of current line in hgfile
            fronttel_over = False
            for linebl in bl:
                start = int(linebl.split()[1])
                end = int(linebl.split()[2])
                for linef in hgfile:
                    if linef.startswith('>chr'):
                        continue
                    if linef in "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n":
                        f_pos += len(linef)
                        position += len(linef)
                        continue
                    if "TTAGGGTTAGGG" in linef.upper(): # end telomere
                        index = linef.upper().find("TTAGGG")
                        outputf.write(linef[:index].upper().replace("TTAGGG",""))
                        outputf.close()
                        return
                    if position+len(linef)-1 < start:
                        # delete telomere sequences, delete N's (only chrM doesn't have telomeres)
                        if not fronttel_over and "CCCTAACCCTAA" in linef.upper(): # front temolere
                            index = linef.upper().find("CCCTAA")
                            outputf.seek(0)
                            outputf.truncate()
                            outputf.write(linef[index:].upper().replace("CCCTAA",""))
                        else:
                            if "NNNNNNNNN" not in linef:
                                fronttel_over = True # once meet a line without telomere sequence, then front telomere is over
                            outputf.write(linef)
                        f_pos += len(linef)
                        position += len(linef)
                    elif position+len(linef)-1 >= start and position+len(linef)-1 <= end:
                        f_pos += len(linef)
                        position += len(linef)
                        if position <= start:
                            linef = linef.strip()[: start-position]
                            if not fronttel_over and "CCCTAACCCTAA" in linef.upper(): # front temolere
                                index = linef.upper().find("CCCTAA")
                                outputf.seek(0)
                                outputf.truncate()
                                outputf.write(linef[index:].upper().replace("CCCTAA",""))
                            else:
                                if "NNNNNNNNN" not in linef:
                                    fronttel_over = True
                                outputf.write(linef)
                        # else delete entire line
                    elif position >= start and position <= end:
                        f_pos += len(linef)
                        position += len(linef)
                        if position+len(linef)-1 >= end:
                            linef = linef.strip()[end-position :]
                            if not fronttel_over and "CCCTAACCCTAA" in linef.upper(): # front temolere
                                index = linef.upper().find("CCCTAA")
                                outputf.seek(0)
                                outputf.truncate()
                                outputf.write(linef[index:].upper().replace("CCCTAA",""))
                            else:
                                if "NNNNNNNNN" not in linef:
                                    fronttel_over = True
                                outputf.write(linef)
                        # else delete entire line
                    elif position > end:
                        # go back one line before break
                        hgfile.seek(f_pos)
                        position += len(linef)
                        break
            # print remaining lines
            for linef in hgfile:
                if linef.startswith('>chr'):
                    continue
                if linef in "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n":
                    f_pos += len(linef)
                    position += len(linef)
                    continue
                if "TTAGGGTTAGGG" in linef.upper(): # end telomere
                    index = linef.upper().find("TTAGGG")
                    outputf.write(linef[:index].upper().replace("TTAGGG",""))
                    outputf.close()
                    return
                else:
                    outputf.write(linef)
            outputf.close()


# (deprecated)
def count_hg_each(chrnum):
    count_ref_chr(chrnum)
    out = pd.DataFrame(0, index=['pos'], columns=['A', 'T', 'C', 'G'])
    with open(temp_folder+"/chr"+str(chrnum)+"_hgfiltered.txt", "r") as hgf:
        for linef in hgf:
            for i in range(0, len(linef)-1):
                if linef[i].upper() != 'N':
                    out.at['pos', linef[i].upper()] += 1
    return out 


def count_ref_chr2(chrnum):
    if hg_file_generated:
        hg_file_individual = temp_folder+'/>chr'+str(chrnum)+'_hgfile.txt'
    else:
        hg_file_individual = hg_filepath+'/chr'+str(chrnum)+'.fa'
    with open(hg_file_individual, 'r') as hgfile:
        # line = hgfile.readline().replace('\n', '')
        # genome = hgfile.read().replace('\n', '')
        # if not line.startswith('>chr'):
        #     line += genome
        #     genome = line

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


def count_hg_old(output_csv_hg):
    df1 = count_hg_each(1)
    df1.to_csv(path_or_buf=output_csv_hg, mode='a')
    df2 = count_hg_each(2)
    df2.to_csv(path_or_buf=output_csv_hg, mode='a')
    df3 = count_hg_each(3)
    df3.to_csv(path_or_buf=output_csv_hg, mode='a')
    df4 = count_hg_each(4)
    df4.to_csv(path_or_buf=output_csv_hg, mode='a')
    df5 = count_hg_each(5)
    df5.to_csv(path_or_buf=output_csv_hg, mode='a')
    df6 = count_hg_each(6)
    df6.to_csv(path_or_buf=output_csv_hg, mode='a')
    df7 = count_hg_each(7)
    df7.to_csv(path_or_buf=output_csv_hg, mode='a')
    df8 = count_hg_each(8)
    df8.to_csv(path_or_buf=output_csv_hg, mode='a')
    df9 = count_hg_each(9)
    df9.to_csv(path_or_buf=output_csv_hg, mode='a')
    df10 = count_hg_each(10)
    df10.to_csv(path_or_buf=output_csv_hg, mode='a')
    df11 = count_hg_each(11)
    df11.to_csv(path_or_buf=output_csv_hg, mode='a')
    df12 = count_hg_each(12)
    df12.to_csv(path_or_buf=output_csv_hg, mode='a')
    df13 = count_hg_each(13)
    df13.to_csv(path_or_buf=output_csv_hg, mode='a')
    df14 = count_hg_each(14)
    df14.to_csv(path_or_buf=output_csv_hg, mode='a')
    df15 = count_hg_each(15)
    df15.to_csv(path_or_buf=output_csv_hg, mode='a')
    df16 = count_hg_each(16)
    df16.to_csv(path_or_buf=output_csv_hg, mode='a')
    df17 = count_hg_each(17)
    df17.to_csv(path_or_buf=output_csv_hg, mode='a')
    df18 = count_hg_each(18)
    df18.to_csv(path_or_buf=output_csv_hg, mode='a')
    df19 = count_hg_each(19)
    df19.to_csv(path_or_buf=output_csv_hg, mode='a')
    df20 = count_hg_each(20)
    df20.to_csv(path_or_buf=output_csv_hg, mode='a')
    df21 = count_hg_each(21)
    df21.to_csv(path_or_buf=output_csv_hg, mode='a')
    df22 = count_hg_each(22)
    df22.to_csv(path_or_buf=output_csv_hg, mode='a')
    dfX = count_hg_each('X')
    dfX.to_csv(path_or_buf=output_csv_hg, mode='a')
    dfY = count_hg_each('Y')
    dfY.to_csv(path_or_buf=output_csv_hg, mode='a')
    dfM = count_hg_each('M')
    dfM.to_csv(path_or_buf=output_csv_hg, mode='a')

    df = df1.copy()
    df.add(df2)
    df.add(df3)
    df.add(df4)
    df.add(df5)
    df.add(df6)
    df.add(df7)
    df.add(df8)
    df.add(df9)
    df.add(df10)
    df.add(df11)
    df.add(df12)
    df.add(df13)
    df.add(df14)
    df.add(df15)
    df.add(df16)
    df.add(df17)
    df.add(df18)
    df.add(df19)
    df.add(df20)
    df.add(df21)
    df.add(df22)
    print(df)
    df.to_csv(path_or_buf=output_csv_hg, mode='a')
    return df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19,df20,df21,df22,dfX,dfY,dfM,df


#print("Usage: python dsbsequence.py temp ../bowtie-files/hg19.fa output blacklist_files $no_control DSB-count-1126 2")
#print("Usage: python dsbsequence.py temp bowtie-files/GRCh38 output blacklist_files $no_control DSB-count-1126 10")
print("Analyzing break sequence bias......")
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

# df_hg = count_hg(output_csv_hg)

for x in xs:
    os.remove(temp_folder+"/chr"+str(x)+"_comparedhits.txt")
    os.remove(temp_folder+"/chr"+str(x)+"_comparedhits_repeats1.txt")
    os.remove(temp_folder+"/chr"+str(x)+"_comparedhits_repeats2.txt")
    os.remove(temp_folder+"/chr"+str(x)+"_comparedhits_repeats_combined.txt")


