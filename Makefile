CXX=g++
CXXFLAGS=-I.
LDLIBS = -lgzstream -lz
DEPS = trio_binning.h

all: classify_longreads classify_pairedend

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(LDLIBS)

classify_longreads: classify_longreads.o kmer.o seq.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDLIBS)

classify_pairedend: classify_pairedend.o kmer.o seq.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDLIBS)

.PHONY: clean

clean:
	rm -f *.o *.gch classify_longreads
