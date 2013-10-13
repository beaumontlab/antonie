CXXFLAGS=-Wall -I.  -MMD -O3
CFLAGS=-I. -O3 -MMD

all: powerdna 16ssearcher

-include *.d

powerdna: powerdna.o hash.o geneannotated.o misc.o fastq.o
	g++ $(LDFLAGS) powerdna.o hash.o geneannotated.o misc.o fastq.o -o $@

16ssearcher: 16ssearcher.o hash.o misc.o fastq.o
	g++ $(LDFLAGS) 16ssearcher.o hash.o misc.o fastq.o -lz -o $@

clean:
	rm -f *~ *.o *.d powerdna

report.html: report-header.html data report-footer.html
	cat report-header.html data report-footer.html > $@	