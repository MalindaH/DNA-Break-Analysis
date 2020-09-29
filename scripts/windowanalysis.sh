#!/bin/bash
# Usage: ./windowanalysis.sh test0825/bowtie-outputHAboth.txt test0825/bowtie-outputHC.txt 48000 outputtest0827 hg19-consensusBlacklist.bed ../genome-annotation/gencode.v34.annotation.gtf ../genome-annotation/Census_all-Sep17-2020.csv

treated=${1}; control=${2}
#wsize=${3}; 
outputname=${4}; blacklist=${5}; annotationfilegtf=${6}; annotationfilecancer=${7}

function filterSort() {
    echo -e "Filtering and sorting alignments by chromosomes......"
    if [ ! -r temp ]; then
        mkdir temp
    fi

    # filter mismatch: (not needed here) | awk 'length($9) == 0'

    for((x=1;x<=22;x++)); do
        # keep alignments of length >= 23nt, filter repetitive reads, and sort in numerical order
        cat $treated | awk 'length($6) > 23' | grep -P "chr${x}\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chr${x}t_hitssorted.txt
        cat $control | awk 'length($6) > 23' | grep -P "chr${x}\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chr${x}c_hitssorted.txt
    done

    cat $treated | awk 'length($6) > 23' | grep  -P "chrX\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrXt_hitssorted.txt
    cat $treated | awk 'length($6) > 23' | grep  -P "chrX\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrXc_hitssorted.txt

    cat $treated | awk 'length($6) > 23' | grep  -P "chrY\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrYt_hitssorted.txt
    cat $treated | awk 'length($6) > 23' | grep  -P "chrY\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrYc_hitssorted.txt
   
    cat $treated | awk 'length($6) > 23' | grep  -P "chrM\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrMt_hitssorted.txt
    cat $treated | awk 'length($6) > 23' | grep  -P "chrM\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrMc_hitssorted.txt
}

function blacklistPrep() {
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
function forMatlab() {
    for((x=1;x<=22;x++)); do
        awk '{print $2,$6}' output/${1}chr${x}.txt > output/${1}chr${x}_pval.txt
        awk '{print $2,$7}' output/${1}chr${x}.txt > output/${1}chr${x}_qval.txt
    done

    awk '{print $2,$6}' output/${1}chrX.txt > output/${1}chrX_pval.txt
    awk '{print $2,$7}' output/${1}chrX.txt > output/${1}chrX_qval.txt
    awk '{print $2,$6}' output/${1}chrY.txt > output/${1}chrY_pval.txt
    awk '{print $2,$7}' output/${1}chrY.txt > output/${1}chrY_qval.txt
    awk '{print $2,$6}' output/${1}chrM.txt > output/${1}chrM_pval.txt
    awk '{print $2,$7}' output/${1}chrM.txt > output/${1}chrM_qval.txt
}

function forGTF() {
    if [ ! -r annotation_files ]; then
        mkdir annotation_files
    fi

    grep -P "gene\t" $annotationfilegtf > annotation_files/annotation-genes.gtf
    grep -P "protein_coding" $annotationfilegtf > annotation_files/annotation-protein-coding.gtf

    for((x=1;x<=22;x++)); do
        grep -P "chr$x\t" annotation_files/annotation-genes.gtf > annotation_files/chr${x}_annotation-genes.gtf
        grep -P "chr$x\t" annotation_files/annotation-protein-coding.gtf > annotation_files/chr${x}_annotation-protein-coding.gtf
    done

    grep -P "chrX\t" annotation_files/annotation-genes.gtf > annotation_files/chrX_annotation-genes.gtf
    grep -P "chrX\t" annotation_files/annotation-protein-coding.gtf > annotation_files/chrX_annotation-protein-coding.gtf

    grep -P "chrY\t" annotation_files/annotation-genes.gtf > annotation_files/chrY_annotation-genes.gtf
    grep -P "chrY\t" annotation_files/annotation-protein-coding.gtf > annotation_files/chrY_annotation-protein-coding.gtf

    grep -P "chrM\t" annotation_files/annotation-genes.gtf > annotation_files/chrM_annotation-genes.gtf
    grep -P "chrM\t" annotation_files/annotation-protein-coding.gtf > annotation_files/chrM_annotation-protein-coding.gtf
    
}

function filterpeaks() {
    for((x=1;x<=22;x++)); do
        awk '$6<=0.05' output0827/${outputname}chr${x}.txt > output/${outputname}chr${x}_peaks.txt
    done

    awk '$6<=0.05' output0827/${outputname}chrX.txt > output/${outputname}chrX_peaks.txt
    awk '$6<=0.05' output0827/${outputname}chrY.txt > output/${outputname}chrY_peaks.txt
    awk '$6<=0.05' output0827/${outputname}chrM.txt > output/${outputname}chrM_peaks.txt
}




# check “<output-filename>.txt” does not exist in the current directory

if [ $# -ne 7 ]; then
    echo "Usage: ./windowanalysis.sh <bowtie-output-treated> <bowtie-output-control> <window-size> <output-filename> <blacklist-file> <gencode.gtf-annotation-file> <cancer-gene-consensus.csv-annotation-file>"
    exit 10
fi

#filterSort

# Remove blacklisted alignments:
#blacklistPrep
#python filterblacklist.py temp blacklist_files

# Calculate p and q values:
if [ ! -r output ]; then
    mkdir output
fi
#python windowanalysis.py temp $3 $4 blacklist_files
# check error code


# generate files for plotting using MATLAB
#forMatlab $4

# generate peak files
#filterpeaks

# generate gene annotation files for gene analysis
#forGTF

python geneanalysis.py output annotation_files ${4} ${7}



