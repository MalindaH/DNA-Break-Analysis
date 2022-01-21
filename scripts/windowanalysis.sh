#!/bin/bash
# Usage: ./windowanalysis.sh -t bowtie-output/bowtie-outputA.txt -c bowtie-output/bowtie-outputC.txt -s 48000 -o outputtest1125 -B genome-annotation/GRCh38_unified_blacklist.bed -C genome-annotation/Census_all-Sep17-2020.csv -G genome-annotation/gencode.v34.annotation.gtf -hg ../bowtie-files/GRCh38
# Usage: ./windowanalysis.sh --no-control -t bowtie-outputA.txt -s 4000 -o outputtest1216 -B ../genome-annotation/GRCh38_unified_blacklist.bed -G ../genome-annotation/gencode.v34.annotation.gtf -R ../genome-annotation/refGene.txt -C ../genome-annotation/Census_all-Sep17-2020.csv -hg ../bowtie-files/GRCh38
# ./windowanalysis.sh -t ../work_dir_eachadapter57/bowtie-output-A2.txt -c ../work_dir_eachadapter57/bowtie-output-A5.txt -o outputtest0520 -B ../genome-annotation/GRCh38_unified_blacklist.bed -G ../genome-annotation/gencode.v34.annotation.gtf -R ../genome-annotation/refGene.txt -C ../genome-annotation/Census_all-Sep17-2020.csv -hg ../bowtie-files/GRCh38 -s 4000
# ./windowanalysis.sh --no-control -t bowtie-output-1.txt -o outputtest0701 -B ../genome-annotation/GRCh38_unified_blacklist.bed -G ../genome-annotation/gencode.v34.annotation.gtf -C ../genome-annotation/Census_all-Sep17-2020.csv -hg ../bowtie-files/GRCh38 -R ../genome-annotation/refGene.txt -r 1 (-s 10000)
# ./windowanalysis.sh -t ../work_dir_eto11-UMI/bowtie-output-1.txt -c ../work_dir_dmso11-UMI/bowtie-output-1.txt -o outputtest0701 -B ../genome-annotation/GRCh38_unified_blacklist.bed -G ../genome-annotation/gencode.v34.annotation.gtf -C ../genome-annotation/Census_all-Sep17-2020.csv -hg ../bowtie-files/GRCh38 -R ../genome-annotation/refGene.txt -r 1 -s 8000
# ./windowanalysis.sh -t bowtie-output-1.txt -c ../work_dir_56-dmso/bowtie-output-1.txt -o outputtest0701 -B ../genome-annotation/GRCh38_unified_blacklist.bed -G ../genome-annotation/GRCh38_latest_genomic.gff -C ../genome-annotation/Census_all-Sep17-2020.csv -hg ../bowtie-files/GRCh38 -R ../genome-annotation/refGene.txt -r 1 -E ../genome-annotation/U2OS_gene_expression.txt -s 10000


function show_help() {
    echo -e "Usage: ./windowanalysis.sh [options]\n\n[options] include:\t
        -t|--treated <bowtie-output-treated> bowtie output file of treated sample
        -c|--control <bowtie-output-control> bowtie output file of control sample
            (provide this or enter --no-control)
        --no-control
        -r|--read-name <number of columns taken by read names in the input files>
        -s|--size <window-size> user-defined window size of peak calling
            (if not provided, find and use window size with highest variance of p-values)
        -o|--output <output-filename>
        -B <blacklist-file>
        -C <cancer-gene-consensus.csv-annotation-file>
        -hg <human-reference-genome-file>
        -R <refGene.txt-annotation-file>
        -E <U2OS_gene_expression.txt-gene-expression-file>
        -G <gencode.gtf-annotation-file>
            (this function is not implemented for use yet)";
}

if [ $# -lt 16 ] || [ $# -gt 22 ]; then  
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
        -r|--read-name) readnamelength="$2"; shift;;
        -s|--size) wsize="$2"; shift ;;
        -o|--output) outputname="$2"; shift ;;
        -B) blacklist="$2"; shift ;;
        -G) annotationfilegtf="$2"; shift ;;
        -R) annotationfilerefgene="$2"; shift ;;
        -C) annotationfilecancer="$2"; shift ;;
        -E) annotationfilegeneexp="$2"; shift ;;
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
echo " -- input file(s) read name length: $readnamelength"
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
echo " -- annotation file (gene expression, .txt): $annotationfilegeneexp"
echo " -- annotation file (gencode, .gtf): $annotationfilegtf"



# filter and sort bowtie output (keep repetitive reads)
function filterSort_KeepRepetitive() {
    echo -e "\n-> Filtering and sorting alignments by chromosomes (keep repetitive reads)......"
    if [ ! -r ${1} ]; then
        mkdir ${1}
    fi

    # the filtered file should have 9 columns; the 5th column is the position number
    if [[ "$readnamelength" -eq 2 ]]; then
        for((x=1;x<=22;x++)); do
            # sort in numerical order
            cat $treated | grep "chr${x}\t" | sort -k 5 -n > ${1}/chr${x}t_hitssorted.txt
            if [[ "$no_control" -eq 0 ]]; then
                cat $control | grep "chr${x}\t" | sort -k 5 -n > ${1}/chr${x}c_hitssorted.txt
            fi
        done

        cat $treated | grep  "chrX\t" | sort -k 5 -n > ${1}/chrXt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep  "chrX\t" | sort -k 5 -n > ${1}/chrXc_hitssorted.txt
        fi

        cat $treated | grep "chrY\t" | sort -k 5 -n > ${1}/chrYt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep "chrY\t" | sort -k 5 -n > ${1}/chrYc_hitssorted.txt
        fi
    
        cat $treated | grep "chrM\t" | sort -k 5 -n > ${1}/chrMt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep "chrM\t" | sort -k 5 -n > ${1}/chrMc_hitssorted.txt
        fi
    elif [[ "$readnamelength" -eq 1 ]]; then
        for((x=1;x<=22;x++)); do
            cat $treated | grep "chr${x}\t" | awk '{ print $1,$1,$2,$3,$4,$5,$6,$7 }' | sort -k 5 -n > ${1}/chr${x}t_hitssorted.txt
            if [[ "$no_control" -eq 0 ]]; then
                cat $control | grep "chr${x}\t" | awk '{ print $1,$1,$2,$3,$4,$5,$6,$7 }' | sort -k 5 -n > ${1}/chr${x}c_hitssorted.txt
            fi
        done

        cat $treated | grep "chrX\t" | awk '{ print $1,$1,$2,$3,$4,$5,$6,$7 }' | sort -k 5 -n > ${1}/chrXt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep "chrX\t" | awk '{ print $1,$1,$2,$3,$4,$5,$6,$7 }' | sort -k 5 -n > ${1}/chrXc_hitssorted.txt
        fi

        cat $treated | grep "chrY\t" | awk '{ print $1,$1,$2,$3,$4,$5,$6,$7 }' | sort -k 5 -n > ${1}/chrYt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep "chrY\t" | awk '{ print $1,$1,$2,$3,$4,$5,$6,$7 }' | sort -k 5 -n > ${1}/chrYc_hitssorted.txt
        fi
    
        cat $treated | grep "chrM\t" | awk '{ print $1,$1,$2,$3,$4,$5,$6,$7 }' | sort -k 5 -n > ${1}/chrMt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep "chrM\t" | awk '{ print $1,$1,$2,$3,$4,$5,$6,$7 }' | sort -k 5 -n > ${1}/chrMc_hitssorted.txt
        fi
    elif [[ "$readnamelength" -eq 3 ]]; then # if used sratoolkit to get SRRxxx.fastq raw data (contains one extra column 3)
        for((x=1;x<=22;x++)); do
            cat $treated | grep "chr${x}\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chr${x}t_hitssorted.txt
            if [[ "$no_control" -eq 0 ]]; then
                cat $control | grep "chr${x}\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chr${x}c_hitssorted.txt
            fi
        done

        cat $treated | grep  "chrX\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrXt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep  "chrX\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrXc_hitssorted.txt
        fi

        cat $treated | grep  "chrY\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrYt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep  "chrY\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrYc_hitssorted.txt
        fi
    
        cat $treated | grep  "chrM\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrMt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep  "chrM\t" | awk '{ print $1,$2,$4,$5,$6,$7,$8,$9 }' | sort -k 5 -n > ${1}/chrMc_hitssorted.txt
        fi
    elif [[ "$readnamelength" -eq 4 ]]; then # BLISS output after fastp removing UMI
        for((x=1;x<=22;x++)); do
            cat $treated | grep "chr${x}\t" | awk '{ print $1,$4,$5,$6,$7,$8,$9,$10 }' | sort -k 5 -n > ${1}/chr${x}t_hitssorted.txt
            if [[ "$no_control" -eq 0 ]]; then
                cat $control | grep "chr${x}\t" | awk '{ print $1,$4,$5,$6,$7,$8,$9,$10 }' | sort -k 5 -n > ${1}/chr${x}c_hitssorted.txt
            fi
        done

        cat $treated | grep  "chrX\t" | awk '{ print $1,$4,$5,$6,$7,$8,$9,$10 }' | sort -k 5 -n > ${1}/chrXt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep  "chrX\t" | awk '{ print $1,$4,$5,$6,$7,$8,$9,$10 }' | sort -k 5 -n > ${1}/chrXc_hitssorted.txt
        fi

        cat $treated | grep  "chrY\t" | awk '{ print $1,$4,$5,$6,$7,$8,$9,$10 }' | sort -k 5 -n > ${1}/chrYt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep  "chrY\t" | awk '{ print $1,$4,$5,$6,$7,$8,$9,$10 }' | sort -k 5 -n > ${1}/chrYc_hitssorted.txt
        fi
    
        cat $treated | grep  "chrM\t" | awk '{ print $1,$4,$5,$6,$7,$8,$9,$10 }' | sort -k 5 -n > ${1}/chrMt_hitssorted.txt
        if [[ "$no_control" -eq 0 ]]; then
            cat $control | grep  "chrM\t" | awk '{ print $1,$4,$5,$6,$7,$8,$9,$10 }' | sort -k 5 -n > ${1}/chrMc_hitssorted.txt
        fi
    else
        echo "read name length invalid, please edit windowanalysis.sh"
        ecit 3 
    fi

    # num_cols=`awk '{print NF}' $treated | head -n 1`
    # # col5=`cat $treated | head -5 | awk '$5 ~ /^chr/'`
    # # if [ ! -z "$col5" ]; then # if used sratoolkit to get SRRxxx.fastq raw data (contains one extra column 3)
    # if [[ "$num_cols" -eq 10 ]]; then # BLISS output after fastp removing UMI
}

function blacklistPrep() {
    if [ ! -r blacklist_files ]; then
        mkdir blacklist_files
    fi

    for((x=1;x<=22;x++)); do
        cat $blacklist | grep "chr${x}\t" > blacklist_files/chr${x}_blacklist.txt
    done

    cat $blacklist | grep "chrX\t" > blacklist_files/chrX_blacklist.txt
    cat $blacklist | grep "chrY\t" > blacklist_files/chrY_blacklist.txt
    cat $blacklist | grep "chrM\t" > blacklist_files/chrM_blacklist.txt
}

# for gencode.v34.annotation.gtf 
function forGTF() {
    if [ ! -r annotation_files ]; then
        mkdir annotation_files
    fi

    grep "gene\t" $annotationfilegtf > annotation_files/annotation-genes.txt
    grep "protein_coding" annotation_files/annotation-genes.txt > annotation_files/annotation-protein-coding.txt

    for((x=1;x<=22;x++)); do
        grep "chr$x\t" annotation_files/annotation-genes.txt > annotation_files/chr${x}_annotation-genes.txt
        grep "chr$x\t" annotation_files/annotation-protein-coding.txt > annotation_files/chr${x}_annotation-protein-coding.txt
    done

    grep "chrX\t" annotation_files/annotation-genes.txt > annotation_files/chrX_annotation-genes.txt
    grep "chrX\t" annotation_files/annotation-protein-coding.txt > annotation_files/chrX_annotation-protein-coding.txt

    grep "chrY\t" annotation_files/annotation-genes.txt > annotation_files/chrY_annotation-genes.txt
    grep "chrY\t" annotation_files/annotation-protein-coding.txt > annotation_files/chrY_annotation-protein-coding.txt

    grep "chrM\t" annotation_files/annotation-genes.txt > annotation_files/chrM_annotation-genes.txt
    grep "chrM\t" annotation_files/annotation-protein-coding.txt > annotation_files/chrM_annotation-protein-coding.txt
}

# for GRCh38_latest_genomic.gff
function forGTF-gff() {
    if [ ! -r annotation_files ]; then
        mkdir annotation_files
    fi

    grep "\tgene\t" $annotationfilegtf > annotation_files/annotation-genes.txt

    # grep "NC_000001.11" annotation_files/annotation-genes.txt > annotation_files/chr1_annotation-genes.txt
    # grep "NC_000002.12" annotation_files/annotation-genes.txt > annotation_files/chr2_annotation-genes.txt
    # grep "NC_000003.12" annotation_files/annotation-genes.txt > annotation_files/chr3_annotation-genes.txt
    # grep "NC_000004.12" annotation_files/annotation-genes.txt > annotation_files/chr4_annotation-genes.txt
    # grep "NC_000005.10" annotation_files/annotation-genes.txt > annotation_files/chr5_annotation-genes.txt
    # grep "NC_000006.12" annotation_files/annotation-genes.txt > annotation_files/chr6_annotation-genes.txt
    # grep "NC_000007.14" annotation_files/annotation-genes.txt > annotation_files/chr7_annotation-genes.txt
    # grep "NC_000008.11" annotation_files/annotation-genes.txt > annotation_files/chr8_annotation-genes.txt
    # grep "NC_000009.12" annotation_files/annotation-genes.txt > annotation_files/chr9_annotation-genes.txt
    # grep "NC_000010.11" annotation_files/annotation-genes.txt > annotation_files/chr10_annotation-genes.txt
    # grep "NC_000011.10" annotation_files/annotation-genes.txt > annotation_files/chr11_annotation-genes.txt
    # grep "NC_000012.12" annotation_files/annotation-genes.txt > annotation_files/chr12_annotation-genes.txt
    # grep "NC_000013.11" annotation_files/annotation-genes.txt > annotation_files/chr13_annotation-genes.txt
    # grep "NC_000014.9" annotation_files/annotation-genes.txt > annotation_files/chr14_annotation-genes.txt
    # grep "NC_000015.10" annotation_files/annotation-genes.txt > annotation_files/chr15_annotation-genes.txt
    # grep "NC_000016.10" annotation_files/annotation-genes.txt > annotation_files/chr16_annotation-genes.txt
    # grep "NC_000017.11" annotation_files/annotation-genes.txt > annotation_files/chr17_annotation-genes.txt
    # grep "NC_000018.10" annotation_files/annotation-genes.txt > annotation_files/chr18_annotation-genes.txt
    # grep "NC_000019.10" annotation_files/annotation-genes.txt > annotation_files/chr19_annotation-genes.txt
    # grep "NC_000020.11" annotation_files/annotation-genes.txt > annotation_files/chr20_annotation-genes.txt
    # grep "NC_000021.9" annotation_files/annotation-genes.txt > annotation_files/chr21_annotation-genes.txt
    # grep "NC_000022.11" annotation_files/annotation-genes.txt > annotation_files/chr22_annotation-genes.txt
    # grep "NC_000023.11" annotation_files/annotation-genes.txt > annotation_files/chrX_annotation-genes.txt
    # grep "NC_000024.10" annotation_files/annotation-genes.txt > annotation_files/chrY_annotation-genes.txt
    # grep "NC_012920.1" annotation_files/annotation-genes.txt > annotation_files/chrM_annotation-genes.txt    
}

function ranksensitivegenes() {
    # cat output/chr1_sensitive-genes.txt output/chr2_sensitive-genes.txt output/chr3_sensitive-genes.txt output/chr4_sensitive-genes.txt output/chr5_sensitive-genes.txt output/chr6_sensitive-genes.txt output/chr7_sensitive-genes.txt output/chr8_sensitive-genes.txt output/chr9_sensitive-genes.txt output/chr10_sensitive-genes.txt output/chr11_sensitive-genes.txt output/chr12_sensitive-genes.txt output/chr13_sensitive-genes.txt output/chr14_sensitive-genes.txt output/chr15_sensitive-genes.txt output/chr16_sensitive-genes.txt output/chr17_sensitive-genes.txt output/chr18_sensitive-genes.txt output/chr19_sensitive-genes.txt output/chr20_sensitive-genes.txt output/chr21_sensitive-genes.txt output/chr22_sensitive-genes.txt output/chrX_sensitive-genes.txt > output/allchr_sensitive-genes.txt
    # sort -k 1 -gr output/allchr_sensitive-genes.txt > output/allchr_sensitive-genes-sorted.txt
    # rm output/chr1_sensitive-genes.txt output/chr2_sensitive-genes.txt output/chr3_sensitive-genes.txt output/chr4_sensitive-genes.txt output/chr5_sensitive-genes.txt output/chr6_sensitive-genes.txt output/chr7_sensitive-genes.txt output/chr8_sensitive-genes.txt output/chr9_sensitive-genes.txt output/chr10_sensitive-genes.txt output/chr11_sensitive-genes.txt output/chr12_sensitive-genes.txt output/chr13_sensitive-genes.txt output/chr14_sensitive-genes.txt output/chr15_sensitive-genes.txt output/chr16_sensitive-genes.txt output/chr17_sensitive-genes.txt output/chr18_sensitive-genes.txt output/chr19_sensitive-genes.txt output/chr20_sensitive-genes.txt output/chr21_sensitive-genes.txt output/chr22_sensitive-genes.txt output/chrX_sensitive-genes.txt

    cat output/chr1_sensitive-cancer-genes.txt output/chr2_sensitive-cancer-genes.txt output/chr3_sensitive-cancer-genes.txt output/chr4_sensitive-cancer-genes.txt output/chr5_sensitive-cancer-genes.txt output/chr6_sensitive-cancer-genes.txt output/chr7_sensitive-cancer-genes.txt output/chr8_sensitive-cancer-genes.txt output/chr9_sensitive-cancer-genes.txt output/chr10_sensitive-cancer-genes.txt output/chr11_sensitive-cancer-genes.txt output/chr12_sensitive-cancer-genes.txt output/chr13_sensitive-cancer-genes.txt output/chr14_sensitive-cancer-genes.txt output/chr15_sensitive-cancer-genes.txt output/chr16_sensitive-cancer-genes.txt output/chr17_sensitive-cancer-genes.txt output/chr18_sensitive-cancer-genes.txt output/chr19_sensitive-cancer-genes.txt output/chr20_sensitive-cancer-genes.txt output/chr21_sensitive-cancer-genes.txt output/chr22_sensitive-cancer-genes.txt output/chrX_sensitive-cancer-genes.txt > output/allchr_sensitive-cancer-genes.txt
    sort -k 1 -gr output/allchr_sensitive-cancer-genes.txt > output/allchr_sensitive-cancer-genes-sorted.txt
    rm output/chr1_sensitive-cancer-genes.txt output/chr2_sensitive-cancer-genes.txt output/chr3_sensitive-cancer-genes.txt output/chr4_sensitive-cancer-genes.txt output/chr5_sensitive-cancer-genes.txt output/chr6_sensitive-cancer-genes.txt output/chr7_sensitive-cancer-genes.txt output/chr8_sensitive-cancer-genes.txt output/chr9_sensitive-cancer-genes.txt output/chr10_sensitive-cancer-genes.txt output/chr11_sensitive-cancer-genes.txt output/chr12_sensitive-cancer-genes.txt output/chr13_sensitive-cancer-genes.txt output/chr14_sensitive-cancer-genes.txt output/chr15_sensitive-cancer-genes.txt output/chr16_sensitive-cancer-genes.txt output/chr17_sensitive-cancer-genes.txt output/chr18_sensitive-cancer-genes.txt output/chr19_sensitive-cancer-genes.txt output/chr20_sensitive-cancer-genes.txt output/chr21_sensitive-cancer-genes.txt output/chr22_sensitive-cancer-genes.txt output/chrX_sensitive-cancer-genes.txt

    cat output/chr*_refgene_counts.txt > output/chrAll_refGene_counts.txt
    rm output/chr*_refgene_counts.txt
}



## --- prep steps --- ##

filterSort_KeepRepetitive temp

## Remove blacklisted alignments: (only need to run blacklistPrep once)
blacklistPrep
python filterblacklist.py temp blacklist_files $no_control

if [ ! -r output ]; then
    mkdir output
fi

## --- for statistics of sequence bias --- ##
python dsbsequence.py temp ${hgfile} output blacklist_files $no_control DSB-motif-output 10

if [ ! -r annotation_files ]; then
    mkdir annotation_files
fi

# ## --- Calculate p values for window of cencer genes, and break density wrt TSS and TTS --- ##
# python geneanalysis.py output annotation_files $outputname $annotationfilecancer temp $no_control $annotationfilerefgene
if [[ "$no_control" -eq 1 ]]; then # no control
    tail -r output/allchr_sensitive-cancer-genes-sorted.txt > output/allchr_sensitive-cancer-genes-sorted-bigtosmall.txt
else # with control
    cat output/allchr_sensitive-cancer-genes-sorted.txt > output/allchr_sensitive-cancer-genes-sorted-bigtosmall.txt
fi
# tail -r output/allchr_sensitive-cancer-genes-sorted.txt > output/allchr_sensitive-cancer-genes-sorted-bigtosmall.txt

## ranksensitivegenes # don't need to run ranksensitivegenes, included in ths python script above

## --- Calculate break density wrt TSS, TTS, wrt gene expression level --- ##
forGTF-gff
python geneanalysis2.py output annotation_files $outputname temp $no_control $annotationfilegeneexp

## --- Calculate p (and q) values for window of user-defined size (or auto-select size) --- ##
python windowanalysis.py temp output $wsize $outputname blacklist_files $no_control



echo "Done :)"



