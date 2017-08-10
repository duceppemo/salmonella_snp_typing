## Description

Salmonella typing based on SNP at 60 sites. Clades are set by comparing the 60 SNP sequence to a reference database.

## Publication

https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-713

## Dependencies

* Java (https://java.com/en/download/manual.jsp)

* BBDuk (http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/)

## Usage

A typical command line using paired-end fastq files:
```
bash ~/scripts/salmonella\_SNP\_PCR.sh \
    -q -m \
    -p SE_clades_primers_21-mer.fasta \
    -o salmonellaPCR/ \
    mysample_R1.fastq.gz
    mysample_R2.fastq.gz
```
