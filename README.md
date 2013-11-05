ANTONIE
=======
Antonie is an integrated, robust, reliable and fast processor of DNA reads,
mostly from Next Generation Sequencing platforms. 

Antonie is open source software.

Initial focus is on automatically & quickly producing the most useful
results on prokaryotic sized genomes. A second goal is to make the program
robust against bad input: out of the box it will refuse to draw conclusions
based on low quality or unnaturally distributed data.

Antonie is named after [Antonie van
Leeuwenhoek](http://en.wikipedia.org/wiki/Antonie_van_Leeuwenhoek), the
Delft inventor of microscopes and the discoverer of bacteria.

AUTHORS
=======
Antonie is developed at the Bertus Beaumontlab of the Bionanoscience
Department of Delft University of Technology.  The lead author is bert
hubert.

For more information, please see:

 * [Bertus Beaumont Lab](http://bertusbeaumontlab.tudelft.nl/)
 * [bert hubert](http://ds9a.nl/)

Please contact us at a.w.r.hubert at tudelft dot nl, or report issues
through our Github page on <https://github.com/beaumontlab/antonie>.

CAPABILITIES
============
Currently, Antonie can map the FASTQ output of sequencers to a FASTA
reference genome.  In addition, it can also exclude known contaminants, like
for example PhiX. Finally, if annotation of the reference genome is available, 
features found by Antonie will be annotated.

The output of Antonie is:

 * A JSON-compatible file with analysis, graphs, data, log, annotations
 * A pretty webpage displaying the JSON data
 * SAM file, mapping the reads to the reference genome

The analysis includes calls for:

 * SNPs
 * Indels
 * Metagenomically variable loci

In addition, there are graphs of:

 * Distribution of reported Phred scores (global, per read position)
 * Distribution of actually measured Phred scores 
 * K-mer variability per read position
 * GC-content per read position

So as a formula:

FASTQ + FASTA + ANNOTATIONS -> JSON + SAM

![First graphs](http://ds9a.nl/antonie/antonie1.png)

![Second graphs](http://ds9a.nl/antonie/antonie2.png)

![Third graphs](http://ds9a.nl/antonie/antonie3.png)

CODE
====
Antonie is written in C++ and has no external dependencies. It can be distributed as 
a single file for Mac, Linux and Windows platforms. 

To compile, get a recent C++ compiler, and run 'make'. 

LIMITATIONS
===========
The current algorithm is fast on common hardware, but needs around 200MB of
memory for a typical prokaryote.  It also assumes it is aligning against a
single chromosome.  Combined, this means that right now, eukaryotic
processing hard to do using Antonie.

SAMPLE USE
==========

> $ antonie -f tot.fastq -r PfSBW25\_genome\_FASTA.fa -b 15 -x phi-x174.fasta -a NC_012660.csv -u > report

This will align the reads from 'tot.fastq' against the Pseudomonas SBW25
reference genome, while stripping out any PhiX reads. Annotations will be read from
'NC_012660.csv'. A human readable, but large, text based report will be written to 'report'.

The mapping will be saved as 'data.sam', and can for example be viewed in
[Tablet](http://bioinf.scri.ac.uk/tablet/) from the James Hutton Institute,
or post-processed using [samtools](http://samtools.sourceforge.net/).

Meanwhile, unmatched reads from 'tot.fastq' will be written to 'unfound.fastq', and could
for example be reprocessed against another reference file to see what is in
there. Alternatively, paste output from 'unfound.fastq' into BLAST.

Finally, in 'data.js', all interesting features found are encoded in JSON format. To view this,
point your browser at 'report.html', and it will source 'data.js' and print pretty graphs.

Sample output:

	Done reading 12125 annotations
	Snipping 15 from beginning of reads, 0 from end of reads
	Reading FASTA reference genome of 'ENA|AM181176|AM181176.4 Pseudomonas fluorescens SBW25 complete genome'
	Indexing 6722540 nucleotides for length 135.. done
	Sorting hashes.. done
	Average hash fill: 1.00445
	Reading FASTA reference genome of 'gi|216019|gb|J02482.1|PX1CG Coliphage phi-X174, complete genome'
	Loading positive control filter genome(s)
	Indexing 5387 nucleotides for length 135.. done
	Sorting hashes.. done
	Average hash fill: 1.0259
	Performing exact matches of reads to reference genome
	Total reads:                                5167210 (0.70 gigabps)
	Excluded control reads:                 -     27917
	Quality excluded:                       -         0
	Ignored reads with N:                   -      5439
	Full matches:                           -   3720893 (72.01%)
	Not fully matched:                      =   1412961 (27.34%)
	Mean Q:                                          34.68 +- 35.41
	Average depth:                                   69.60
	Undercovered nucleotides:                      1393 (0.02%), 73 ranges
	Indexing 6722540 nucleotides for length 11.. done
	Sorting hashes.. done
	Average hash fill: 2.9893
	Indexing 5387 nucleotides for length 11.. done
	Sorting hashes.. done
	Average hash fill: 1.0041
	Performing sliding window partial matches
	Performing sliding window partial matches for control/exclude
	Fuzzy found:                            -   1106081
	Fuzzy found in excluded set:            -     13714
	Unmatchable reads:                      =    298605 (5.78%)
	Average depth:                                   84.84
	Undercovered nucleotides:                       613 (0.01%), 37 ranges
	Found 147036 varying loci
	Found 25 seriously variable loci
	Found 872 loci with at least one insert in a read
	Found 3 serious inserts


