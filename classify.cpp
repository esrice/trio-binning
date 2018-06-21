#include <set>
#include <string>
#include <iostream>
#include <stdint.h>
#include <fstream>
#include <stdlib.h>

struct haplotype_counts_t {
    int hapA_count, hapB_count;
};

/**
 * Convert a kmer in string form to a bunch of bits, where
 * A=0, C=1, G=2, T=3.
 */
uint64_t kmer_to_bits(const std::string& kmer_string)
{
    unsigned char i;
    uint64_t two_bits;
    uint64_t bit_repr = 0;
    for (i = 0; i < kmer_string.length(); i++)
    {
        switch(kmer_string[i]) {
            case 'A': two_bits = 0; break;
            case 'C': two_bits = 1; break;
            case 'G': two_bits = 2; break;
            case 'T': two_bits = 3; break;
        }

        bit_repr |= two_bits << (2*i);
    }

    return bit_repr;
}

std::string bits_to_kmer(uint64_t bit_repr, size_t k)
{
    uint64_t two_bits;
    char base;
    std::string kmer (k, 'N');
    uint64_t mask = 3;

    int i;
    for (i=0; i<k; i++) {
        switch((bit_repr & (mask<<(2*i))) >>(2*i) )
        {
            case 0: base = 'A'; break;
            case 1: base = 'C'; break;
            case 2: base = 'G'; break;
            case 3: base = 'T'; break;
        }

        kmer[i] = base;
    }

    return kmer;
}

size_t get_kmer_size(char* file_path)
{
    std::ifstream infile;
    std::string kmer;

    infile.open(file_path);
    if (!infile) {
        std::cerr << "Can't open " << file_path;
        exit(1);
    }

    infile >> kmer;
    infile.close();
    return kmer.length();
}

haplotype_counts_t count_kmers_in_read(std::string read,
        std::set<uint64_t>& hapA_kmers, std::set<uint64_t>& hapB_kmers,
        size_t k)
{
    int i;
    std::string kmer;
    haplotype_counts_t counts;
    counts.hapA_count = 0; counts.hapB_count = 0;

    for (i=0; i < read.length()-k+1; i++) {
        kmer = read.substr(i, k);
        if (hapA_kmers.find(kmer_to_bits(kmer)) != hapA_kmers.end()) {
            counts.hapA_count++;
        }

        if (hapB_kmers.find(kmer_to_bits(kmer)) != hapB_kmers.end()) {
            counts.hapB_count++;
        }
    }

    return counts;
}

std::set<uint64_t> read_kmers_into_set(char* file_path)
{
    std::set<uint64_t> kmers;
    std::string kmer;
    std::ifstream infile;

    infile.open(file_path);
    if (!infile) {
        std::cerr << "Can't open " << file_path;
        exit(1);
    }

    // TODO figure out how to use the insert function that gives a "hint" to
    // optimize stuff, however that works
    while (infile >> kmer) {
        kmers.insert(kmer_to_bits(kmer));
    }
    infile.close();

    return kmers;
}

int main(int argc, char** argv)
{
    std::set<uint64_t> hapA_kmers, hapB_kmers;
    std::set<uint64_t>::iterator it;

    size_t k = get_kmer_size(argv[1]);
    uint8_t *seq;

    int i;

    // read k-kmers into sets
    hapA_kmers = read_kmers_into_set(argv[1]);
    hapB_kmers = read_kmers_into_set(argv[2]);

    // TODO go through reads and classify. perhaps require fastq?

    // this is test code
    haplotype_counts_t counts = count_kmers_in_read(
            "AAAAAAAAAAAAAAAAAAAAATTTATATATTATTTTATTATT", hapA_kmers,
            hapB_kmers, k);
    std::cout << counts.hapA_count << ' ' << counts.hapB_count << '\n';
    counts = count_kmers_in_read("CCACCATTACATTGAGGAGTCGCTACTGATCAGGGGACTGCA",
            hapA_kmers, hapB_kmers, k);
    std::cout << counts.hapA_count << ' ' << counts.hapB_count << '\n';
    counts = count_kmers_in_read("TACAGCGGAGCGCCCAAATTT",
            hapA_kmers, hapB_kmers, k);
    std::cout << counts.hapA_count << ' ' << counts.hapB_count << '\n';
}
