#!/bin/bash
# Usage: ./windowanalysis.sh bowtie-outputA.txt bowtie-outputC.txt 48000 output.txt hg19-consensusBlacklist.bed

treated=${1}; control=${2}
#wsize=${3}; outputname=${4}
blacklist=${5}

function filterSort() {
    if [ ! -r temp ]; then
        mkdir temp
    fi

    # keep alignments of length >= 23nt, with no mismatch
    cat $treated | grep + | awk 'length($6) > 23' | awk 'length($9) == 0' > temp/allt+hits.txt
    cat $control | grep + | awk 'length($6) > 23' | awk 'length($9) == 0' > temp/allc+hits.txt
    cat $treated | grep - | awk 'length($6) > 23' | awk 'length($9) == 0' > temp/allt-hits.txt
    cat $control | grep - | awk 'length($6) > 23' | awk 'length($9) == 0' > temp/allc-hits.txt

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

    rm temp/allt+hits.txt temp/allc+hits.txt temp/allt-hits.txt temp/allc-hits.txt
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



echo -e "Filtering and sorting alignments by chromosomes......"
filterSort
echo -e "Removing blacklisted alignments......"
#blacklistPrep
#python filterblacklist.py temp blacklist_files
#rm temp/allt+hits.txt temp/allc+hits.txt temp/allt-hits.txt temp/allc-hits.txt
echo -e "Calculating p and q values......"
#python windowanalysis.py temp $3 $4



