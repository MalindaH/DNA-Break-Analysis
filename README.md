# DSB-Analysis
A repository to track scripts I wrote for DSB analysis. <br>
Scripts follow steps of iSeq Hygestat_BLESS (http://breakome.utmb.edu/software.html). <br>

**System:**
 - OS: Windows subsystem for Linux, Ubuntu
 - Bowtie version 1.2.3
 - Python version 3.8.2
 - MATLAB R2019a

### Step 0.0: Download programs
Bowtie: http://bowtie-bio.sourceforge.net/manual.shtml


### Step 0.1: Download preparatory files
1. Download human reference genome file hg19.fa: https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies#grch37

2. Build the reference genome files: `bowtie-build hg19.fasta bowtie-files/hg19.bowtie`

3. Download BLESS raw sequencing files: (https://www.ebi.ac.uk/ena) SRR codes: SRR695402, SRR695426 (aphidicolin-treated), SRR695321 (control)
 ```
 cat SRR695402.fastq, SRR695426.fastq > DSB-quantification/treated.fastq
 mv SRR695321.fastq DSB-quantification/control.fastq
 ```
4. Download blacklist file of HUMAN (hg19/GRCh37): https://sites.google.com/site/anshulkundaje/projects/blacklists
- Blacklisted Regions: set of regions in the human genome that have anomalous, unstructured, high signal/read counts in next gen sequencing experiments independent of cell line and type of experiment, which are troublesome for high throughput sequencing aligners.


### Step 1: Align to reference genome using Bowtie
```
bowtie -q bowtie-files/hg19.bowtie -v 0 -M 1 --best DSB-quantification/treated.fastq bowtie-outputT.txt
bowtie -q bowtie-files/hg19.bowtie -v 0 -M 1 --best DSB-quantification/control.fastq bowtie-outputC.txt
```
`-q` means input file is in .fastq format, `-v 0` allows 0 mismatch, `-M 1 --best` reports the best read if a read has more than 1 reportable alignments

### Step 2: Filter alignments and analyze alignment
Usage: `./windowanalysis.sh <bowtie-output-treated> <bowtie-output-control> <window-size> <output-filename> <blacklist-file>`

Detailed steps in `windowanalysis.sh`:
1. Filter bowtie outputs: keep alignments of length >= 23 nt, with no mismatch, and filter repetitive reads; delete alignments in blacklisted regions using `filterblacklist.py`

2. Analyze each chromosome by windows of user-defined size, calculate hypergeometric p-values, and do Benjamini-Hochberg correction to produce q-values using `windowanalysis.py`
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
&nbsp;

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


