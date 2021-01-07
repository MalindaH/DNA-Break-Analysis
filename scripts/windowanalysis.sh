#!/bin/bash
# Usage: ./windowanalysis.sh -t bowtie-output/bowtie-outputA.txt -c bowtie-output/bowtie-outputC.txt -s 48000 -o outputtest1125 -B genome-annotation/GRCh38_unified_blacklist.bed -C genome-annotation/Census_all-Sep17-2020.csv -G genome-annotation/gencode.v34.annotation.gtf -hg ../bowtie-files/GRCh38
# Usage: ./windowanalysis.sh --no-control -t bowtie-outputA.txt -s 4000 -o outputtest1216 -B ../genome-annotation/GRCh38_unified_blacklist.bed -G ../genome-annotation/gencode.v34.annotation.gtf -R ../genome-annotation/refGene.txt -C ../genome-annotation/Census_all-Sep17-2020.csv -hg ../bowtie-files/GRCh38

function show_help() {
    echo -e "Usage: ./windowanalysis.sh [options]\n\n[options] include:\t
        -t|--treated <bowtie-output-treated> bowtie output file of treated sample
        -c|--control <bowtie-output-control> bowtie output file of control sample
            (provide this or enter --no-control)
        --no-control
        -s|--size <window-size> user-defined window size of peak calling
            (if not provided, find and use window size with highest variance of p-values)
        -o|--output <output-filename>
        -B <blacklist-file>
        -C <cancer-gene-consensus.csv-annotation-file>
        -hg <human-reference-genome-file>
        -R <refGene.txt-annotation-file>
        -G <gencode.gtf-annotation-file>
            (this function is not implemented for use yet)";
}

if [ $# -lt 14 ] && [ $# -gt 16 ]; then  
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
        -s|--size) wsize="$2"; shift ;;
        -o|--output) outputname="$2"; shift ;;
        -B) blacklist="$2"; shift ;;
        -G) annotationfilegtf="$2"; shift ;;
        -R) annotationfilerefgene="$2"; shift ;;
        -C) annotationfilecancer="$2"; shift ;;
        -hg) hgfile="$2"; shift ;;
        *) echo "Unknown parameter passed: $1" ;;
    esac
    shift
done

echo -e "\nInputs:"
echo " -- treated (bowtie output): $treated"
if [[ "$no_control" -eq 1 ]]; then
    echo " -- no control"
else
    echo " -- control (bowtie output): $control"
fi
if [ $wsize ]; then
    echo " -- window size: $wsize"
else
    wsize=-1
    echo " -- window size: use window size with highest variance of p-values"
fi
echo " -- output name: $outputname"
echo " -- human reference genome file or folder (.fa): $hgfile"
echo " -- blacklist file (.bed): $blacklist"
echo " -- annotation file (cancer census, .csv): $annotationfilecancer"
echo " -- annotation file (refGene, .txt): $annotationfilerefgene"
echo " -- annotation file (gencode, .gtf): $annotationfilegtf"



# filter and sort bowtie output (keep repetitive reads)
function filterSort_KeepRepetitive() {
    echo -e "\n-> Filtering and sorting alignments by chromosomes (keep repetitive reads)......"
    if [ ! -r ${1} ]; then
        mkdir ${1}
    fi

    col5=`cat $treated | head -5 | awk '$5 ~ /^chr/'`
    if [ ! -z "$col5" ]; then # if used sratoolkit to get SRRxxx.fastq raw data (contains one extra column 3)
        for((x=1;x<=22;x++)); do
            # sort in numerical order
            cat $treated | grep -P "chr${x}\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chr${x}t_hitssorted.txt
            if [[ "$no_control" -ne 1 ]]; then
                cat $control | grep -P "chr${x}\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chr${x}c_hitssorted.txt
            fi
        done

        cat $treated | grep  -P "chrX\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrXt_hitssorted.txt
        if [[ "$no_control" -ne 1 ]]; then
            cat $control | grep  -P "chrX\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrXc_hitssorted.txt
        fi

        cat $treated | grep  -P "chrY\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrYt_hitssorted.txt
        if [[ "$no_control" -ne 1 ]]; then
            cat $control | grep  -P "chrY\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrYc_hitssorted.txt
        fi
    
        cat $treated | grep  -P "chrM\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrMt_hitssorted.txt
        if [[ "$no_control" -ne 1 ]]; then
            cat $control | grep  -P "chrM\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrMc_hitssorted.txt
        fi
    else    
        for((x=1;x<=22;x++)); do
            # sort in numerical order
            cat $treated | grep -P "chr${x}\t" | sort -k 5 -n > ${1}/chr${x}t_hitssorted.txt
            if [[ "$no_control" -ne 1 ]]; then
                cat $control | grep -P "chr${x}\t" | sort -k 5 -n > ${1}/chr${x}c_hitssorted.txt
            fi
        done

        cat $treated | grep  -P "chrX\t" | sort -k 5 -n > ${1}/chrXt_hitssorted.txt
        if [[ "$no_control" -ne 1 ]]; then
            cat $control | grep  -P "chrX\t" | sort -k 5 -n > ${1}/chrXc_hitssorted.txt
        fi

        cat $treated | grep  -P "chrY\t" | sort -k 5 -n > ${1}/chrYt_hitssorted.txt
        if [[ "$no_control" -ne 1 ]]; then
            cat $control | grep  -P "chrY\t" | sort -k 5 -n > ${1}/chrYc_hitssorted.txt
        fi
    
        cat $treated | grep  -P "chrM\t" | sort -k 5 -n > ${1}/chrMt_hitssorted.txt
        if [[ "$no_control" -ne 1 ]]; then
            cat $control | grep  -P "chrM\t" | sort -k 5 -n > ${1}/chrMc_hitssorted.txt
        fi
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




## --- prep steps --- ##

# filterSort_KeepRepetitive temp

## Remove blacklisted alignments: (only need to run blacklistPrep once)
# blacklistPrep
# python filterblacklist.py temp blacklist_files $no_control

if [ ! -r output ]; then
    mkdir output
fi

## --- for statistics of break sequence motif --- ##
# python dsbsequence.py temp ${hgfile} output blacklist_files $no_control DSB-count-1218 10

if [ ! -r annotation_files ]; then
    mkdir annotation_files
fi

## --- Calculate p values for window of cencer genes, and break density wrt TSS and TTS--- ##
python geneanalysis.py output annotation_files $outputname $annotationfilecancer temp $no_control $annotationfilerefgene


## --- Calculate p (and q) values for window of user-defined size --- ##
# python windowanalysis.py temp output $wsize $outputname blacklist_files $no_control



echo "Done :)"



