import sys
import os
import collections
import pandas as pd

def compare_tc(chrnum):
    # delete sequences slightly divergent owing to sequencing errors
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
                    outputf.write(lt.split()[3]+'\t'+lt.split()[4]+'\t'+lt.split()[2]+'\t'+lt.split()[5]+'\t'+str(len(lt.split()[5]))+'\n')
                    f_pos += len(lt)
                    continue
        # print remaining lines
        for lt in tfile:
            outputf.write(lt.split()[3]+'\t'+lt.split()[4]+'\t'+lt.split()[2]+'\t'+lt.split()[5]+'\t'+str(len(lt.split()[5]))+'\n')
        outputf.close()


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
    # slight_divergent(chrnum)

def slight_divergent(chrnum):
     with open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats2.txt", 'r') as tfile:
        outputf = open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats_combined.txt", "a+")
        prev_line = ''
        temp = []
        n = []
        for current_line in tfile: 
            if prev_line != '':
                if int(prev_line.split()[1])+12 > int(current_line.split()[1]):
                    temp.append(prev_line.split()[1])
                    n.append(int(prev_line.split()[5]))
                else:
                    if n: # if n not empty
                        temp.append(prev_line.split()[1])
                        n.append(int(prev_line.split()[5]))
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
            n.append(int(prev_line.split()[5]))
            index = n.index(max(n))
            outputf.write(prev_line.split()[0]+'\t'+temp[index]+'\t'+str(sum(n))+'\n')
            temp = []
            n = []
        else:
            outputf.write(prev_line.split()[0]+'\t'+str(prev_line.split()[1])+'\t1\n')
        outputf.close()

def generate_hg_files():
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


# output 4x5 array (dataframe): out[0,:] is count of 4 bases (A,T,C,G) for position -2
# out[1,:] is count of 4 bases (A,T,C,G) for position -1
# out[2,:] is count of 4 bases (A,T,C,G) for position 0 (break position)
# out[3,:] is count of 4 bases (A,T,C,G) for position +1
# out[4,:] is count of 4 bases (A,T,C,G) for position +2
# out[:,0] is count of A, out[:,1] is count of T, out[:,2] is count of C, out[:,3] is count of G
def generate_chr_sequences(chrnum):
    with open(temp_folder+"/chr"+str(chrnum)+"_comparedhits_repeats_combined.txt", 'r') as hitsfile:
        with open(temp_folder+'/>chr'+str(chrnum)+'_hgfile.txt', 'r') as hgfile:
            # outputf = open(output_folder+'/chr'+str(chrnum)+'.txt', 'a+')
            pos = 0
            f_pos = 0
            bug_count = 0
            # out = np.zeros((5,4))
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
    print('chr'+str(chrnum)+':')
    print(out)
    print()
    return out


def add_count(fivebases, df):
    if fivebases.find('N') != -1:
        return df
    df.at['pos -2', fivebases[0]] += 1
    df.at['pos -1', fivebases[1]] += 1
    df.at['pos 0', fivebases[2]] += 1
    df.at['pos +1', fivebases[3]] += 1
    df.at['pos +2', fivebases[4]] += 1
    return df


#print("Usage: python dsbsequence.py temp1 ../bowtie-files/hg19.fa output")
print("Counting sequenses near DSB......")
temp_folder = sys.argv[1]
hg_filepath = sys.argv[2]
output_folder = sys.argv[3]

# x=1
# while x<=22:
#     compare_tc(x)
#     x += 1
# compare_tc('X')
# compare_tc('Y')
# compare_tc('M')

# x=1
# while x<=22:
#     repetitive_reads(x)
#     x += 1

# x=1
# while x<=22:
#     slight_divergent(x)
#     x += 1

# generate_hg_files()

print('......generating counts......')
df1 = generate_chr_sequences(1)
df2 = generate_chr_sequences(2)
df3 = generate_chr_sequences(3)
df4 = generate_chr_sequences(4)
df5 = generate_chr_sequences(5)
df6 = generate_chr_sequences(6)
df7 = generate_chr_sequences(7)
df8 = generate_chr_sequences(8)
df9 = generate_chr_sequences(9)
df10 = generate_chr_sequences(10)
df11 = generate_chr_sequences(11)
df12 = generate_chr_sequences(12)
df13 = generate_chr_sequences(13)
df14 = generate_chr_sequences(14)
df15 = generate_chr_sequences(15)
df16 = generate_chr_sequences(16)
df17 = generate_chr_sequences(17)
df18 = generate_chr_sequences(18)
df19 = generate_chr_sequences(19)
df20 = generate_chr_sequences(20)
df21 = generate_chr_sequences(21)
df22 = generate_chr_sequences(22)


df1.to_csv(path_or_buf=output_folder+'/DSB-count-1109-2.csv', mode='a')
df2.to_csv(path_or_buf=output_folder+'/DSB-count-1109-2.csv', mode='a')
df3.to_csv(path_or_buf=output_folder+'/DSB-count-1109-2.csv', mode='a')

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
df.to_csv(path_or_buf=output_folder+'/DSB-count-1109-2.csv', mode='a')




# with open(temp_folder+'/>chr1_hgfile.txt', 'r') as hgfile:
#     for line in hgfile:
#         print(len(line))


# x=1
# while x<=22:
#     generate_chr_sequences(x)
#     x += 1





