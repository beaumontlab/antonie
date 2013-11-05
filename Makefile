-include sysdeps/$(shell uname).inc

CXXFLAGS=-Wall -I. -Iext/libmba  -MMD -O3 $(CXX2011FLAGS) 
CFLAGS=-I. -Iext/libmba -O3 -MMD 
LDFLAGS=$(CXX2011FLAGS) 

all: antonie 16ssearcher

-include *.d

ANTONIE_OBJECTS = antonie.o hash.o geneannotated.o misc.o fastq.o saminfra.o ext/libmba/libdiff.a

antonie: $(ANTONIE_OBJECTS)
	$(CXX) $(ANTONIE_OBJECTS) $(LDFLAGS) -o $@

16ssearcher: 16ssearcher.o hash.o misc.o fastq.o
	$(CXX) $(LDFLAGS) 16ssearcher.o hash.o misc.o fastq.o -lz -o $@

clean:
	rm -f *~ *.o *.d powerdna

