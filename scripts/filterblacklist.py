import sys
import os


def remove_blacklist(chrnum, kind):
    with open(bl_folder+"/chr"+str(chrnum)+"_blacklist.txt", "r") as bl:
        with open(chr_folder+"/chr"+str(chrnum)+kind+"_hitssorted.txt", "r") as f:  
            outputf = open(chr_folder+"/chr"+str(chrnum)+kind+"_hitsfiltered.txt", "a+")
            counter = 0
            f_pos = 0
            for linebl in bl:
                start = int(linebl.split()[1])
                end = int(linebl.split()[2])
                #print("start = "+str(start)+", end = "+str(end))
                for linef in f:
                    position = int(linef.split()[4])
                    if position < start:
                        if linef.split()[5] not in "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG" and linef.split()[5] not in "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA" and "TTAGGGTTAGGGTTAGGG" not in linef.split()[5] and "CCCTAACCCTAACCCTAA" not in linef.split()[5]: # filter out telomere sequences
                            outputf.write(linef)
                        f_pos += len(linef)
                    elif position >= start and position <= end:
                        counter += 1
                        f_pos += len(linef)
                        continue
                    elif position > end:
                        # go back one line before break
                        f.seek(f_pos)
                        break
            # print remaining lines
            for linef in f:
                if "TTAGGGTTAGGGTTAGGG" not in linef.split()[5]  and "CCCTAACCCTAACCCTAA" not in linef.split()[5] and linef.split()[5] not in "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG" and linef.split()[5] not in "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA": # filter out telomere sequences
                    outputf.write(linef)
            outputf.close()
            #print("deleted "+str(counter)+" lines")
    os.remove(chr_folder+"/chr"+str(chrnum)+kind+"_hitssorted.txt")



#print("Usage: python filterblacklist.py temp blacklist_files $no_control")
print("\n-> Removing blacklisted alignments......")
chr_folder = sys.argv[1]
bl_folder = sys.argv[2]
no_control = sys.argv[3]

x = 1
while x <= 22:
    remove_blacklist(x, "t")
    if no_control == '0':
        remove_blacklist(x, "c")
    x+=1
    
remove_blacklist('X', "t")
if no_control == '0':
    remove_blacklist('X', "c")
remove_blacklist('Y', "t")
if no_control == '0':
    remove_blacklist('Y', "c")
remove_blacklist('M', "t")
if no_control == '0':
    remove_blacklist('M', "c")





