# 16SSEARCHER

Sequencing runs ideally only contain the intended DNA. Often however, they
are contaminated with unintended bacteria.

16ssearcher matches FASTQ files to databases of 16s rRNA, and generates
a list of bacteria (or archea or fungi) present in your reads.

## Databases
For [Green Genes](http://greengenes.secondgenome.com/downloads), 
download gg_xx_y.fasta.gz and gg_xx_y_accessions.txt.gz.

For the [Ribosomal Datbase Project](http://rdp.cme.msu.edu/), download the
unaligned FASTA.gz.

There is no need to decompress the files, leave them as is.

## Syntax

The --help output is instructive, but in short:

$ 16ssearcher --mode gg gg_13_5.fasta.gz gg_13_5_accessions.txt.gz *.fastq

or

$ 16ssearcher --mode rdp release11_1_Bacteria_unaligned.fa.gz *.fastq

A running output is emitted to stderr, a final summary to stdout.

