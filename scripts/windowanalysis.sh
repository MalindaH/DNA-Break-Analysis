#!/bin/bash
# Usage: ./windowanalysis.sh -t bowtie-output/bowtie-outputA.txt -c bowtie-output/bowtie-outputC.txt -s 48000 -o outputtest1125 -B genome-annotation/GRCh38_unified_blacklist.bed -C genome-annotation/Census_all-Sep17-2020.csv -G genome-annotation/gencode.v34.annotation.gtf -hg bowtie-files/GRCh38
# Usage: ./windowanalysis.sh --no-control -t bowtie-outputA.txt -s 4000 -o outputtest1216 -B genome-annotation/GRCh38_unified_blacklist.bed -G genome-annotation/gencode.v34.annotation.gtf -C genome-annotation/Census_all-Sep17-2020.csv -hg bowtie-files/GRCh38

function show_help() {
    echo -e "Usage: ./windowanalysis.sh [options]\n\n[options] include: (all are mandatory)\t
        --no-control\n\t-t | --treated <bowtie-output-treated>\n\t-c | --control <bowtie-output-control>\n\t-s <window-size>\t
        -o | --output <output-filename>\n\t-B <blacklist-file>\n\t-G <gencode.gtf-annotation-file>\t
        -C <cancer-gene-consensus.csv-annotation-file> \n\t-hg <human-reference-genome-file>";
}

if [ $# -ne 15 ] && [ $# -ne 16 ]; then  
    show_help;
    exit 2
fi

no_control=0;

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) show_help; exit 1 ;;
        --no-control) no_control=1 ;;
        -t|--treated) treated="$2"; shift ;;
        -c|--control) control="$2"; shift ;;
        -s) wsize="$2"; shift ;;
        -o|--output) outputname="$2"; shift ;;
        -B) blacklist="$2"; shift ;;
        -G) annotationfilegtf="$2"; shift ;;
        -C) annotationfilecancer="$2"; shift ;;
        -hg) hgfile="$2"; shift ;;
        *) echo "Unknown parameter passed: $1" ;;
    esac
    shift
done

echo "Inputs:"
echo " -- treated (bowtie output): $treated"
if [ $no_control ]; then
    echo " -- no control"
else
    echo " -- control (bowtie output): $control"
fi
echo " -- window size: $wsize"
echo " -- output name: $outputname"
echo " -- blacklist file (.bed): $blacklist"
echo " -- annotation file (gencode, .gtf): $annotationfilegtf"
echo " -- annotation file (cancer census, .csv): $annotationfilecancer"
echo " -- human reference genome file or folder (.fa): $hgfile"

## check files exist


function filterSort() {
    echo -e "Filtering and sorting alignments by chromosomes......"
    if [ ! -r temp ]; then
        mkdir temp
    fi

    # filter mismatch: (not needed here) | awk 'length($9) == 0'
    # keep alignments of length >= 23nt, awk 'length($6) > 23' (not for now)

    for((x=1;x<=22;x++)); do
        # filter repetitive reads, and sort in numerical order
        cat $treated | grep -P "chr${x}\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chr${x}t_hitssorted.txt
        if [ ! $no_control ]; then
            cat $control | grep -P "chr${x}\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chr${x}c_hitssorted.txt
        fi
    done

    cat $treated | grep  -P "chrX\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrXt_hitssorted.txt
    if [ ! $no_control ]; then
        cat $control | grep  -P "chrX\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrXc_hitssorted.txt
    fi

    cat $treated | grep  -P "chrY\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrYt_hitssorted.txt
    if [ ! $no_control ]; then
        cat $control | grep  -P "chrY\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrYc_hitssorted.txt
    fi

    cat $treated | grep  -P "chrM\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrMt_hitssorted.txt
    if [ ! $no_control ]; then
        cat $control | grep  -P "chrM\t" | awk '!seen[$5]++' | sort -k 5 -n > temp/chrMc_hitssorted.txt
    fi
}

function filterSort_KeepRepetitive() {
    echo -e "Filtering and sorting alignments by chromosomes (keep repetitive reads)......"
    if [ ! -r ${1} ]; then
        mkdir ${1}
    fi

    for((x=1;x<=22;x++)); do
        # sort in numerical order
        cat $treated | grep -P "chr${x}\t" | sort -k 5 -n > ${1}/chr${x}t_hitssorted.txt
        if [ ! $no_control ]; then
            cat $control | grep -P "chr${x}\t" | sort -k 5 -n > ${1}/chr${x}c_hitssorted.txt
        fi
    done

    cat $treated | grep  -P "chrX\t" | sort -k 5 -n > ${1}/chrXt_hitssorted.txt
    if [ ! $no_control ]; then
        cat $control | grep  -P "chrX\t" | sort -k 5 -n > ${1}/chrXc_hitssorted.txt
    fi

    cat $treated | grep  -P "chrY\t" | sort -k 5 -n > ${1}/chrYt_hitssorted.txt
    if [ ! $no_control ]; then
        cat $control | grep  -P "chrY\t" | sort -k 5 -n > ${1}/chrYc_hitssorted.txt
    fi
   
    cat $treated | grep  -P "chrM\t" | sort -k 5 -n > ${1}/chrMt_hitssorted.txt
    if [ ! $no_control ]; then
        cat $control | grep  -P "chrM\t" | sort -k 5 -n > ${1}/chrMc_hitssorted.txt
    fi
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
        awk '{print $2,$6}' temp/chr${x}_pval.txt > output/${1}chr${x}_pval.txt
        # awk '{print $2,$7}' output/${1}chr${x}.txt > output/${1}chr${x}_qval.txt
    done

    awk '{print $2,$6}' temp/chrX_pval.txt > output/${1}chrX_pval.txt
    # awk '{print $2,$7}' output/${1}chrX.txt > output/${1}chrX_qval.txt
    awk '{print $2,$6}' temp/chrY_pval.txt > output/${1}chrY_pval.txt
    # awk '{print $2,$7}' output/${1}chrY.txt > output/${1}chrY_qval.txt
    awk '{print $2,$6}' temp/chrM_pval.txt > output/${1}chrM_pval.txt
    # awk '{print $2,$7}' output/${1}chrM.txt > output/${1}chrM_qval.txt
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
        awk '$6<=1e-13' temp/chr${x}_pval.txt > temp/${outputname}chr${x}_peaks.txt
    done

    awk '$6<=1e-13' temp/chrX_pval.txt > temp/${outputname}chrX_peaks.txt
    awk '$6<=1e-13' temp/chrY_pval.txt > temp/${outputname}chrY_peaks.txt
    awk '$6<=1e-13' temp/chrM_pval.txt > temp/${outputname}chrM_peaks.txt
}

function ranksensitivegenes() {
    # cat output/chr1_sensitive-genes.txt output/chr2_sensitive-genes.txt output/chr3_sensitive-genes.txt output/chr4_sensitive-genes.txt output/chr5_sensitive-genes.txt output/chr6_sensitive-genes.txt output/chr7_sensitive-genes.txt output/chr8_sensitive-genes.txt output/chr9_sensitive-genes.txt output/chr10_sensitive-genes.txt output/chr11_sensitive-genes.txt output/chr12_sensitive-genes.txt output/chr13_sensitive-genes.txt output/chr14_sensitive-genes.txt output/chr15_sensitive-genes.txt output/chr16_sensitive-genes.txt output/chr17_sensitive-genes.txt output/chr18_sensitive-genes.txt output/chr19_sensitive-genes.txt output/chr20_sensitive-genes.txt output/chr21_sensitive-genes.txt output/chr22_sensitive-genes.txt output/chrX_sensitive-genes.txt > output/allchr_sensitive-genes.txt
    # sort -k 1 -nr output/allchr_sensitive-genes.txt > output/allchr_sensitive-genes-sorted.txt

    cat output/chr1_sensitive-cancer-genes.txt output/chr2_sensitive-cancer-genes.txt output/chr3_sensitive-cancer-genes.txt output/chr4_sensitive-cancer-genes.txt output/chr5_sensitive-cancer-genes.txt output/chr6_sensitive-cancer-genes.txt output/chr7_sensitive-cancer-genes.txt output/chr8_sensitive-cancer-genes.txt output/chr9_sensitive-cancer-genes.txt output/chr10_sensitive-cancer-genes.txt output/chr11_sensitive-cancer-genes.txt output/chr12_sensitive-cancer-genes.txt output/chr13_sensitive-cancer-genes.txt output/chr14_sensitive-cancer-genes.txt output/chr15_sensitive-cancer-genes.txt output/chr16_sensitive-cancer-genes.txt output/chr17_sensitive-cancer-genes.txt output/chr18_sensitive-cancer-genes.txt output/chr19_sensitive-cancer-genes.txt output/chr20_sensitive-cancer-genes.txt output/chr21_sensitive-cancer-genes.txt output/chr22_sensitive-cancer-genes.txt output/chrX_sensitive-cancer-genes.txt > output/allchr_sensitive-cancer-genes.txt
    sort -k 1 -nr output/allchr_sensitive-cancer-genes.txt > output/allchr_sensitive-cancer-genes-sorted.txt
}

function debug() {
    for((x=1;x<=22;x++)); do
        cat BLESS_output.txt | grep -P "chr${x}\t" | awk '$6<=0.05'  > output/${outputname}chr${x}_peaks.txt
    done

    cat BLESS_output.txt | grep -P "chrX\t" | awk '$6<=0.05' > output/${outputname}chrX_peaks.txt
    cat BLESS_output.txt | grep -P "chrY\t" | awk '$6<=0.05' > output/${outputname}chrY_peaks.txt
    cat BLESS_output.txt | grep -P "chrM\t" | awk '$6<=0.05' > output/${outputname}chrM_peaks.txt
}



# filterSort
# filterSort_KeepRepetitive temp

## Remove blacklisted alignments:
# blacklistPrep
# python filterblacklist.py temp blacklist_files $no_control

# Calculate p and q values:
if [ ! -r output ]; then
    mkdir output
fi
# python windowanalysis.py temp output $wsize $outputname blacklist_files $no_control
## check error code

```
## generate files for plotting using MATLAB
# forMatlab $outputname

## generate peak files
# filterpeaks
```

## generate gene annotation files for gene analysis, only need to run once
# forGTF

# python geneanalysis.py output annotation_files $outputname $annotationfilecancer temp $no_control

ranksensitivegenes

## for debugging:
# debug





## for statistics of DSB sequence bias
# filterSort_KeepRepetitive temp
# blacklistPrep
# python filterblacklist.py temp1 blacklist_files $no_control

# python dsbsequence.py temp ${hgfile} output blacklist_files $no_control /DSB-count-1126.csv /DSB-count-hg-1126.csv