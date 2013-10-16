CFLAGS=-Wall -O3 -I. -MMD

all: libdiff.a

-include *.d

clean:
	rm -f *.o *.a *.d *~

libdiff.a: suba.o diff.o allocator.o msgno.o varray.o
	ar rc libdiff.a suba.o diff.o allocator.o msgno.o varray.o