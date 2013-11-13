-include sysdeps/$(shell uname).inc

CXXFLAGS=-Wall -I. -Iext/libmba -MMD -O3 $(CXX2011FLAGS)
CFLAGS=-I. -Iext/libmba -O3 -MMD 
LDFLAGS=$(CXX2011FLAGS) 

PROGRAMS=antonie 16ssearcher

all: gitversion ext/libmba/libdiff.a $(PROGRAMS)

-include *.d
-include remake-please

.PHONY: gitversion

gitversion: 
	@./update-git-version-if-necessary

gitversion.h:
	./update-git-version-if-necessary

ext/libmba/libdiff.a:
	cd ext/libmba/; make

ANTONIE_OBJECTS = antonie.o hash.o geneannotated.o misc.o fastq.o saminfra.o dnamisc.o ext/libmba/libdiff.a 

strdiff: strdiff.o ext/libmba/libdiff.a
	$(CC) strdiff.o ext/libmba/libdiff.a -o $@

antonie: $(ANTONIE_OBJECTS)
	$(CXX) $(ANTONIE_OBJECTS) $(LDFLAGS) -o $@

16ssearcher: 16ssearcher.o hash.o misc.o fastq.o
	$(CXX) $(LDFLAGS) 16ssearcher.o hash.o misc.o fastq.o -lz -o $@

clean:
	rm -f *~ *.o *.d $(PROGRAMS)
	cd ext/libmba/;	make clean


