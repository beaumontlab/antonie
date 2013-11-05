ANTONIE
=======
Antonie is an integrated, robust, reliable and fast processor of DNA reads,
mostly from Next Generation Sequencing platforms. 

Antonie is open source software.

Initial focus is on automatically & quickly producing the most useful
results on prokaryotic sized genomes.

Antonie is named after Antonie van Leeuwenhoek, the Delft inventor of
microscopes and the discoverer of bacteria.

AUTHORS
=======
Antonie is developed at the Bertus Beaumontlab of the Bionanoscience
Department of Delft University of Technology.  The lead author is bert
hubert.

http://bertusbeaumontlab.tudelft.nl/
http://ds9a.nl/

CAPABILITIES
============
Currently, Antonie can map the FASTQ output of sequencers to a FASTA
reference genome.  In addition, it can also exclude known contaminants, like
for example PhiX. Finally, if annotation of the reference genome is available, 
features found by Antonie will be annotated.

The output of Antonie is:
 * A JSON-compatible file with analysis, graphs, data, log, annotations
 * A pretty webpage displaying the JSON data
 * SAM file

So as a formula:

FASTQ + FASTA + ANNOTATIONS -> JSON + SAM

CODE
====
Antonie is written in C++ and has no external dependencies. It can be distributed as 
a single file for Mac, Linux and Windows platforms. 

LIMITATIONS
===========
The current algorithm is fast on common hardware, but needs around 200MB of
memory for a typical prokaryote.  It also assumes it is aligning against a
single chromosome.  Combined, this means that eukaryotic processing is
currently hard to do using Antonie.

