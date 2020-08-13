import sys
import os


def remove_blacklist(chrnum, kind):
    with open(bl_folder+"/chr"+str(chrnum)+"_blacklist.txt", "r") as bl:
        with open(chr_folder+"/chr"+str(chrnum)+kind+"hitssorted.txt", "r") as f:  
            outputf = open(chr_folder+"/chr"+str(chrnum)+kind+"hitsfiltered.txt", "a+")
            counter = 0
            f_pos = 0
            for linebl in bl:
                start = int(linebl.split()[1])
                end = int(linebl.split()[2])
                #print("start = "+str(start)+", end = "+str(end))
                for linef in f:
                    position = int(linef.split()[4])
                    if position < start:
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
                outputf.write(linef)
            outputf.close()
            #print("deleted "+str(counter)+" lines")
    os.remove(chr_folder+"/chr"+str(chrnum)+kind+"hitssorted.txt")



print("Usage: python filterblacklist.py temp blacklist_files")
chr_folder = sys.argv[1]
bl_folder = sys.argv[2]
x = 1
while x <= 22:
    remove_blacklist(x, "t+")
    remove_blacklist(x, "c+")
    remove_blacklist(x, "t-")
    remove_blacklist(x, "c-")
    x+=1
    
remove_blacklist('X', "t+")
remove_blacklist('X', "c+")
remove_blacklist('X', "t-")
remove_blacklist('X', "c-")
remove_blacklist('Y', "t+")
remove_blacklist('Y', "c+")
remove_blacklist('Y', "t-")
remove_blacklist('Y', "c-")
remove_blacklist('M', "t+")
remove_blacklist('M', "c+")
remove_blacklist('M', "t-")
remove_blacklist('M', "c-")


