CC=g++
CFLAGS=-I.
DEPS = trio_binning.h

all: classify_longreads classify_pairedend

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

classify_longreads: classify_longreads.o kmer.o seq.o
	$(CC) -o $@ $^ $(CFLAGS)

classify_pairedend: classify_pairedend.o kmer.o seq.o
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o *.gch classify_longreads
