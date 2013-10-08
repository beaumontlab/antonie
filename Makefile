CXXFLAGS=-Wall -I. -O3 -MMD
CFLAGS=-I. -O3 -MMD

all: powerdna

-include *.d

powerdna: powerdna.o hash.o
	g++ $(LDFLAGS) powerdna.o hash.o -o $@
	
clean:
	rm -f *~ *.o *.d powerdna
	