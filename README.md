ANTONIE
=======
Antonie is an integrated, robust, reliable and fast processor of DNA reads,
mostly from Next Generation Sequencing platforms (typically Illumina, but we
strive to be multiplatform).  It is currently focussed on prokaryotic and
other small genomes.

Antonie is free open source software, and we welcome contributions!

Initial focus is on automatically & quickly producing the most useful
results on prokaryotic sized genomes. A second goal is to make the program
robust against bad input: out of the box it should refuse to draw conclusions
based on low quality or unnaturally distributed data. 

Antonie is named after [Antonie van
Leeuwenhoek](http://en.wikipedia.org/wiki/Antonie_van_Leeuwenhoek), the
Delft inventor of microscopes and the discoverer of bacteria.

AUTHORS
=======
Antonie is developed at the Bertus Beaumontlab of the Bionanoscience
Department of Delft University of Technology.  The lead author is Bert
Hubert.

For more information, please see:

 * [Bionanoscience department](http://www.tnw.tudelft.nl/en/about-faculty/departments/bionanoscience/)
 * [Bertus Beaumont Lab](http://bertusbeaumontlab.tudelft.nl/)
 * [bert hubert](http://ds9a.nl/)

Please contact us at a.w.r.hubert at tudelft dot nl, or report issues
through our Github page on <https://github.com/beaumontlab/antonie>.

CAPABILITIES
============
Currently, Antonie can map the FASTQ output of sequencers to a FASTA
reference genome.  In addition, it can also exclude known contaminants, like
for example PhiX.  Finally, if GFF3 annotation of the reference genome is
available, features found by Antonie will be annotated.

Antonie performs similar functions as for example
[bowtie](http://bowtie-bio.sourceforge.net/index.shtml), except somewhat
faster for small genomes, while also performing some of the analysis usually
performed further downstream, for example by
[fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

The output of Antonie is:

 * A JSON-compatible file with analysis, graphs, data, log, annotations
 * A pretty webpage displaying the JSON data ([sample](http://ds9a.nl/antonie/SRR956947/report.html))
 * SAM file, mapping the reads to the reference genome

The analysis includes calls for:

 * SNPs ('undercovered regions')
 * Indels
 * Metagenomically variable loci

In addition, there are graphs of:

 * Distribution of reported Phred scores (global, per read position)
 * Distribution of actually measured Phred scores 
 * Q-Q plot of empirical versus reported Phred scores
 * K-mer variability per read position
 * GC-content per read position
 * GC-content distribution of reads (versus genome wide)
 * Duplication count of reads

So as a formula:

FASTQ + FASTA + GFF3 -> JSON + SAM -> PRETTY HTML

![First graphs](http://ds9a.nl/antonie/antonie1.png)

![Second graphs](http://ds9a.nl/antonie/antonie2.png)

![Third graphs](http://ds9a.nl/antonie/antonie3.png)

CODE
====
Antonie is written in C++ and has no external dependencies. It can be distributed as 
a single file for Mac, Linux and Windows platforms. 

To compile, get a recent C++ compiler, and run 'make':

	$ git clone git@github.com:beaumontlab/antonie.git
	$ cd antonie
	$ make

Antonie depends on Boost being installed at compile time, but not at runtime.

LIMITATIONS
===========

 * The current algorithm is fast on common hardware, but needs around 200MB of
   memory for a typical prokaryote.  It also assumes it is aligning against a
   single chromosome.  Combined, this means that right now, eukaryotic
   processing is hard to do using Antonie.
 * Additionally, Antonie does not yet benefit from the position information that
   can be inferred from paired end reads, and in fact you'll have to concatenate both
   paired reads first.
 * Antonie can't yet deal with reads of varying lengths.
 * We only do indels of 1 nucleotide as of now
 * Our HTML is known not to be displayed correctly on some versions of Internet Explorer.
 * Finally, Antonie may be of limited use for reads shorter than 75
   nucleotides, as it hasn't been tried in that domain a lot yet.

SAMPLE USE
==========

> $ antonie -f tot.fastq -r PfSBW25\_genome\_FASTA.fna -x phi-x174.fasta -a NC_012660.gff -u > report

This will align the reads from 'tot.fastq' against the Pseudomonas SBW25
reference genome, while stripping out any PhiX reads. Annotations will be read from
'NC_012660.gff'. A human readable, but large, text based report will be written to 'report'.

The mapping will be saved as 'data.sam', and can for example be viewed in
[Tablet](http://bioinf.scri.ac.uk/tablet/) from the James Hutton Institute,
or post-processed using [samtools](http://samtools.sourceforge.net/).

Meanwhile, because we passed -u, unmatched reads from 'tot.fastq' will be
written to 'unfound.fastq', and could for example be reprocessed against
another reference file to see what is in there.  Alternatively, paste output
from 'unfound.fastq' into BLAST.

Finally, in 'data.js', all interesting features found are encoded in JSON format. To view this,
point your browser at 'report.html', and it will source 'data.js' and print pretty graphs.

Try 'antonie --help' for a full listing of options.

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

GETTING SAMPLE DATA
===================
Reference materials can be found from <ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/>  
The .fna file is FASTA, and corresponds to our '-f' and '-x' fields.  
The .gff file is GFF3, and contains annotations understood by our '-a' field.

As an example, for "Escherichia coli str. K-12 substr. MG1655", head to
<ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/>

To get sample FASTQ, read the 'species' line from the .gff file, and visit 
the URL found there.  This may present you with a series of substrains, or
directly give you a listing of Entrez records.

From there, find SRA or Sequence Read Archive files. SRA is a compressed representation 
of FASTQ, to convert SRA to FASTQ, use 'fastq_dump' the SRA toolkit which can be found on 
<http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software>.

For K12 MG1655, this may work: <http://www.ncbi.nlm.nih.gov/sra/SRX339396>, which will lead you to:
<ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR956/SRR956947>

Once the SRA is converted to FASTQ, run Antonie like this:

> $ antonie -f SRR956947.fastq -a Escherichia\_coli\_K\_12\_substr\_\_MG1655\_uid57779/\*.gff -r Escherichia\_coli\_K\_12\_substr\_\_MG1655\_uid57779/\*.fna > report

And the resulting report should look like this <http://ds9a.nl/antonie/SRR956947/report.html>

FUTURE DIRECTION
================

 * We are aiming for compatibility with the [Galaxy
   Project](http://galaxyproject.org/).
 * The program needs to automate its detection of quality, and not draw conclusions on bad data, but 
   suggest filtering instead: "Antonie is unhappy with the Phred scores at positions < 12".
 * BAM file output (smaller, compressed version of SAM)
 * Automated post-analysis of unmatched reads against Golomb compressed set of common contaminants
 * Utilize the extra information in paired-end reads
 * Detect indels of >1 nucleotide
