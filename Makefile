-include sysdeps/$(shell uname).inc

CXXFLAGS=-Wall -I. -Iext/libmba -MMD -O3 $(CXX2011FLAGS) -Wno-unused-local-typedefs 
CFLAGS=-I. -Iext/libmba -O3 -MMD  
LDFLAGS=$(CXX2011FLAGS)  
CHEAT_ARG := $(shell ./update-git-hash-if-necessary)

PROGRAMS=antonie 16ssearcher

all: $(PROGRAMS)

-include *.d

MBA_OBJECTS = ext/libmba/allocator.o  ext/libmba/diff.o  ext/libmba/msgno.o  ext/libmba/suba.o  ext/libmba/varray.o 
ANTONIE_OBJECTS = antonie.o hash.o geneannotated.o misc.o fastq.o saminfra.o dnamisc.o githash.o $(MBA_OBJECTS)
STATICFLAGS=-Wl,-Bstatic -lstdc++ $(LUALIBS) -lgcc -Wl,-Bdynamic -static-libgcc -lm -lc

# STATICFLAGS= -static -static-libgcc -static-libstdc++

strdiff: strdiff.o $(MBA_OBJECTS)
	$(CC) strdiff.o $(MBA_OBJECTS) -o $@

antonie: $(ANTONIE_OBJECTS)
	$(CXX) $(ANTONIE_OBJECTS) $(LDFLAGS) $(STATICFLAGS) -o $@

SEARCHER_OBJECTS=16ssearcher.o hash.o misc.o fastq.o githash.o

16ssearcher: $(SEARCHER_OBJECTS)
	$(CXX) $(LDFLAGS) $(SEARCHER_OBJECTS) -lz -o $@

clean:
	rm -f *~ *.o $(MBA_OBJECTS) *.d $(PROGRAMS) githash.h 


