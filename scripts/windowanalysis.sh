#!/bin/bash
# Usage: ./windowanalysis.sh test0825/bowtie-outputHAboth.txt test0825/bowtie-outputHC.txt 48000 outputtest08262 hg19-consensusBlacklist.bed

treated=${1}; control=${2}
#wsize=${3}; outputname=${4}
blacklist=${5}

function filterSort() {
    if [ ! -r temp ]; then
        mkdir temp
    fi

    # keep alignments of length >= 23nt
    cat $treated | grep + | awk 'length($6) > 23' > temp/allt+hits.txt
    cat $control | grep + | awk 'length($6) > 23' > temp/allc+hits.txt
    cat $treated | grep - | awk 'length($6) > 23' > temp/allt-hits.txt
    cat $control | grep - | awk 'length($6) > 23' > temp/allc-hits.txt
    # filter mismatch: (not needed here) | awk 'length($9) == 0'

    for((x=1;x<=22;x++)); do
        # filter repetitive reads, and sort in numerical order
        cat temp/allt+hits.txt | grep  -P "chr${x}\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chr${x}t+hitssorted.txt
        cat temp/allc+hits.txt | grep  -P "chr${x}\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chr${x}c+hitssorted.txt
        cat temp/allt-hits.txt | grep  -P "chr${x}\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chr${x}t-hitssorted.txt
        cat temp/allc-hits.txt | grep  -P "chr${x}\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chr${x}c-hitssorted.txt
    done

    cat temp/allt+hits.txt | grep  -P "chrX\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrXt+hitssorted.txt
    cat temp/allc+hits.txt | grep  -P "chrX\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrXc+hitssorted.txt
    cat temp/allt-hits.txt | grep  -P "chrX\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrXt-hitssorted.txt
    cat temp/allc-hits.txt | grep  -P "chrX\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrXc-hitssorted.txt

    cat temp/allt+hits.txt | grep  -P "chrY\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrYt+hitssorted.txt
    cat temp/allc+hits.txt | grep  -P "chrY\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrYc+hitssorted.txt
    cat temp/allt-hits.txt | grep  -P "chrY\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrYt-hitssorted.txt
    cat temp/allc-hits.txt | grep  -P "chrY\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrYc-hitssorted.txt
   
    cat temp/allt+hits.txt | grep  -P "chrM\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrMt+hitssorted.txt
    cat temp/allc+hits.txt | grep  -P "chrM\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrMc+hitssorted.txt
    cat temp/allt-hits.txt | grep  -P "chrM\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrMt-hitssorted.txt
    cat temp/allc-hits.txt | grep  -P "chrM\t" | awk '!seen[$6]++' | sort -k 5 -n > temp/chrMc-hitssorted.txt

    #rm temp/allt+hits.txt temp/allc+hits.txt temp/allt-hits.txt temp/allc-hits.txt
}

function blacklistPrep(){
    if [ ! -r blacklist_files ]; then
        mkdir blacklist_files
    fi

    for((x=1;x<=22;x++)); do
        cat $blacklist | grep -P "chr${x}\t" > blacklist_files/chr${x}_blacklist.txt
    done

    cat $blacklist | grep -P "chrX\t" > blacklist_files/chrX_blacklist.txt
    cat $blacklist | grep -P "chrY\t" > blacklist_files/chrY_blacklist.txt
    cat $blacklist | grep -P "chrM\t" > blacklist_files/chrM_blacklist.txt
}

# generate files for plotting using MATLAB
function forMatlab(){
    for((x=1;x<=22;x++)); do
        awk '{print $2,$6}' output/${1}chr${x}_fwd.txt > output/${1}chr${x}_fwdpval.txt
        awk '{print $2,$7}' output/${1}chr${x}_fwd.txt > output/${1}chr${x}_fwdqval.txt
    done

    awk '{print $2,$6}' output/${1}chrX_fwd.txt > output/${1}chrX_fwdpval.txt
    awk '{print $2,$6}' output/${1}chrY_fwd.txt > output/${1}chrY_fwdpval.txt
    awk '{print $2,$6}' output/${1}chrM_fwd.txt > output/${1}chrM_fwdpval.txt
}



# check “<output-filename>_fwd.txt” and “<output-filename>_rev.txt” does not exist in the current directory

if [ $# -ne 5 ]; then
    echo "Usage: ./windowanalysis.sh <bowtie-output-treated> <bowtie-output-control> <window-size> <output-filename> <blacklist-file>"
    exit 1
fi

echo -e "Filtering and sorting alignments by chromosomes......"
#filterSort

# Remove blacklisted alignments:
#blacklistPrep
#python filterblacklist.py temp blacklist_files

# Calculate p and q values:
if [ ! -r output ]; then
    mkdir output
fi
python windowanalysis.py temp $3 $4 blacklist_files

# generate files for plotting using MATLAB
forMatlab $4


