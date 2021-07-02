# ./umiprocess.sh ../fastq_eto-rep1/SRR11119501.fastq NNNNNNNNXXXXXXXX ^GTCGTCGC ../bowtie-files/hg38.bowtie

echo "Processing UMI..."

if [ ! -r UMI-process-files ]; then
    mkdir UMI-process-files
fi

umi_tools extract --stdin=$1 --bc-pattern=$2 --log=UMI-process-files/processed.log --stdout UMI-process-files/processed.fastq

cutadapt -q 10 -g $3 -o UMI-process-files/cutadapt_out.fastq UMI-process-files/processed.fastq -m 10 --untrimmed-output UMI-process-files/cutadapt_out_untrimmed.fastq &> UMI-process-files/catadapt.log

bowtie -q $4 -v 0 -M 1 --threads 4 --best UMI-process-files/cutadapt_out.fastq --sam UMI-process-files/mapped.sam &> UMI-process-files/bowtie-log-1.txt

samtools view -bS UMI-process-files/mapped.sam > UMI-process-files/mapped.bam

samtools sort UMI-process-files/mapped.bam -o UMI-process-files/example.bam

samtools index UMI-process-files/example.bam

umi_tools dedup --output-stats=deduplicated -I UMI-process-files/example.bam -S UMI-process-files/deduplicated.bam > UMI-process-files/deduplicated.log

samtools bam2fq UMI-process-files/deduplicated.bam > UMI-process-files/deduplicated.fastq

bowtie -q $4 -v 0 -M 1 --best UMI-process-files/deduplicated.fastq bowtie-output-1.txt &> UMI-process-files/bowtie-log-deduplicated.txt
