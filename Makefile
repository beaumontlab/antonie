-include sysdeps/$(shell uname).inc

CXXFLAGS=-Wall -I. -Iext/libmba -MMD -O3 $(CXX2011FLAGS)
CFLAGS=-I. -Iext/libmba -O3 -MMD 
LDFLAGS=$(CXX2011FLAGS) 
CHEAT_ARG := $(shell ./update-git-hash-if-necessary)

PROGRAMS=antonie 16ssearcher

all: ext/libmba/libdiff.a $(PROGRAMS)

-include *.d

ext/libmba/libdiff.a:
	cd ext/libmba/; make

ANTONIE_OBJECTS = antonie.o hash.o geneannotated.o misc.o fastq.o saminfra.o dnamisc.o ext/libmba/libdiff.a githash.o

strdiff: strdiff.o ext/libmba/libdiff.a
	$(CC) strdiff.o ext/libmba/libdiff.a -o $@

antonie: $(ANTONIE_OBJECTS)
	$(CXX) $(ANTONIE_OBJECTS) $(LDFLAGS) -o $@

SEARCHER_OBJECTS=16ssearcher.o hash.o misc.o fastq.o githash.o

16ssearcher: $(SEARCHER_OBJECTS)
	$(CXX) $(LDFLAGS) $(SEARCHER_OBJECTS) -lz -o $@

clean:
	rm -f *~ *.o *.d $(PROGRAMS) githash.h remake-please
	cd ext/libmba/;	make clean


