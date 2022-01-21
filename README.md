# DSB-Analysis
Scripts I wrote for DNA double-strand break analysis.
Referred to methods of:
- iSeq Hygestat_BLESS (http://breakome.utmb.edu/software.html)
- MACS peak calling (https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html)

**System:**
 - OS: Windows subsystem for Linux, Ubuntu
 - Bowtie version 1.2.3
 - Python version 3.8.2
 - Google Colab / Jupyter Notebook
 - (MATLAB R2019a)

### Step 0.0: Download programs
Bowtie: http://bowtie-bio.sourceforge.net/manual.shtml


### Step 0.1: Download and process preparatory files
1. Download and unzip human reference genome file hg38 (GRCh38) .fasta files (chr1.fa, chr2.fa, ..., chrX.fa, chrY.fa, chrM.fa, from January 2014): http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/, (http://hgdownload.soe.ucsc.edu/downloads.html#human)

2. `bowtie-build` to build the reference genome files: `bowtie-build chr1.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chr10.fa,chr11.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17.fa,chr18.fa,chr19.fa,chr20.fa,chr21.fa,chr22.fa,chrX.fa,chrY.fa,chrM.fa ../hg38.bowtie`

3. Download raw sequencing (.fastq) files: (https://www.ebi.ac.uk/ena) using accession numbers
- Here I use files from DSBCapture as an example:
> Lensing, S., Marsico, G., Hänsel-Hertsch, R. et al. DSBCapture: in situ capture and sequencing of DNA breaks. Nat Methods 13, 855–857 (2016). https://doi.org/10.1038/nmeth.3960

4. Download blacklist file of HUMAN (hg38/GRCh38): https://sites.google.com/site/anshulkundaje/projects/blacklists
- Blacklisted Regions: set of regions in the human genome that have anomalous, unstructured, high signal/read counts in next gen sequencing experiments independent of cell line and type of experiment, which are troublesome for high throughput sequencing aligners.

### Step 1: Remove Illumina adaptors (AGATCGGAAGAGC) using cutadapt (https://cutadapt.readthedocs.io/)
- for paired-end illumina sequencing:
```
cutadapt -m 10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --discard-untrimmed -o break_out1.fastq -p break_out2.fastq SRR3182242_1.fastq SRR3182242_2.fastq > catadapt-output.txt
cat break_out1.fastq break_out2.fastq > break_out.fastq
```
- for single-end illumina sequencing:
```
cutadapt -m 10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --discard-untrimmed -o break_out.fastq SRR10173728.fastq > catadapt-output.txt
```

### Step 2: Align to reference genome using Bowtie
```
bowtie -q ../bowtie-files/hg38.bowtie -v 0 -M 1 --best ../NHEK_break_seq_rep1/break_out.fastq bowtie-outputA.txt
```
`-q` means input file is in .fastq format, `-v 0` allows 0 mismatch, `-M 1 --best` reports the best read if a read has more than 1 reportable alignments

### Step 3: Filter alignments and analyze alignment
Script file: `scripts/windowanalysis.sh`
```
Usage: ./windowanalysis.sh [options]

[options] include:
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
        -G <gencode.gtf-annotation-file>
            (this function is not implemented yet)
```
Example: `./windowanalysis.sh --no-control -t bowtie-outputA.txt -s 5000 -o outputtest -B ../genome-annotation/GRCh38_unified_blacklist.bed -G ../genome-annotation/gencode.v34.annotation.gtf -C ../genome-annotation/Census_all-Sep17-2020.csv -hg ../bowtie-files/GRCh38`

Detailed steps in `windowanalysis.sh`:
1. Filter bowtie outputs: delete alignments in blacklisted regions (using `scripts/filterblacklist.py`)

2. Analyze sequence bias of the break: can modify window size, by default analyze ±10bp of the break (using `scripts/dsbsequence.py`)


2. Analyze each chromosome by windows of user-defined size, and calculate p-values (using `scripts/windowanalysis.py`)
    1) Scan each chromosome by windows of user-defined size
&nbsp;

    2) hypergeometric test for p-values
        * For p-values:
          * N = total number in population = number of reads in the given chromosome for both treated and non-treated sample
          * k = total number with condition in population = number of reads in the given chromosome for the treated sample
          * m = number in subset = number of reads in the given window for treated and non-treated samples
          * x = number with condition in subset = number of reads in a given window for the treated sample
        * Document: https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html
        * Online p-value calculator: https://stattrek.com/online-calculator/hypergeometric.aspx

    3) Benjamini-Hochberg correction to produce q-values (q-values = corrected p-values for multiple hypothesis testing)
        * The corrected P value for a test is either the raw P value times m/i or the adjusted P value for the next higher raw P value, whichever is smaller (m = number of tests, i = rank of each test, with 1 the rank of the smallest P value)
        * Document: http://www.biostathandbook.com/multiplecomparisons.html
        * Online B-H correction calculator: https://www.sdmproject.com/utilities/?show=FDR

**Output:**
- Column 1: chromosome number
- Column 2: window start position
- Column 3: window end position
- Column 4: number of alignments within the window in treated sample
- Column 5: number of alignments within the window in control sample
- Column 6: hypergeometric p-value of the window
- Column 7: Benjamini-Hochberg correction q-value of the window


### Step 3: Filter significant intervals (q-value <= 0.05)
```
awk ‘$7 < 0.05’ output/<output-filename>chr<chrnum>_fwd.txt
awk ‘$7 < 0.05’ output/<output-filename>chr<chrnum>_rev.txt
```

### Step 4: Visualize significant intervals using MATLAB
`chrplot_pval.m`


