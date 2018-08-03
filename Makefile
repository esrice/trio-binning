CC=g++
CFLAGS=-I.
DEPS = trio_binning.h

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

classify_longreads: classify_longreads.o kmer.o seq.o
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm *.o classify_longreads
