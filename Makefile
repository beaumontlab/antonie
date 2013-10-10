CXXFLAGS=-Wall -I.  -MMD -O3
CFLAGS=-I. -O3 -MMD

all: powerdna

-include *.d

powerdna: powerdna.o hash.o
	g++ $(LDFLAGS) powerdna.o hash.o -o $@
	
clean:
	rm -f *~ *.o *.d powerdna

report.html: report-header.html data report-footer.html
	cat report-header.html data report-footer.html > $@	