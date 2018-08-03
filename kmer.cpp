/**
 * This file contains functions for manipulating k-mers.
 * See classify.cpp's main() for usage examples, and
 * trio_binning.h for documentation.
 */

#include "trio_binning.h"

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

std::string reverse_complement(std::string sequence) {
    std::string revcomp_seq  (sequence.length(), '*');
    int i;
    char base;

    for (i = 0; i < sequence.length(); i++) {
        switch (sequence[i]) {
            case 'A': base = 'T'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
            case 'T': base = 'A'; break;
        }
        revcomp_seq[sequence.length() - 1 - i] = base;
    }

    return revcomp_seq;
}

std::string get_canonical_representation(std::string sequence) {
    std::string revcomp_seq = reverse_complement(sequence);
    if (sequence.compare(revcomp_seq) < 0) {
        return sequence;
    } else {
        return revcomp_seq;
    }
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
        kmer = get_canonical_representation(kmer);
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

    while (infile >> kmer) {
        kmers.insert(kmer_to_bits(kmer));
    }
    infile.close();

    return kmers;
}
