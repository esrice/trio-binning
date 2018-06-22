#include <set>
#include <stdint.h>
#include <cstdlib>
#include "trio_binning.h"

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
    unsigned int num_hapA_kmers, num_hapB_kmers, max_num_kmers;
    double scaling_factor_A, scaling_factor_B;
    float hapA_score, hapB_score;
    char best_haplotype;

    size_t k = get_kmer_size(argv[1]);
    uint8_t *seq;

    int i;

    // read k-kmers into sets
    hapA_kmers = read_kmers_into_set(argv[1]);
    hapB_kmers = read_kmers_into_set(argv[2]);

    // look at the sizes of the k-mer sets and use these to calculating scaling
    // factors. The original program divides read haplotype counts by the size
    // of the k-mer set for that haplotype, resulting in a really tiny number.
    // In order to make it more pleasant for humans to look at, we multiply
    // both scaling factors by the size of the larger k-mer set so that the
    // scaling factors are close to 1.
    num_hapA_kmers = hapA_kmers.size();
    num_hapB_kmers = hapB_kmers.size();
    max_num_kmers = std::max(num_hapA_kmers, num_hapB_kmers);
    scaling_factor_A = (double) max_num_kmers / (double) num_hapA_kmers;
    scaling_factor_B = (double) max_num_kmers / (double) num_hapB_kmers;

    // go through reads
    // TODO add gzip support and automatic fasta/fastq detection
    FastqParser reads_parser(argv[3]);
    seq_entry_t entry; // TODO convert this to a class
    haplotype_counts_t counts;

    // TODO actually output reads into files in addition to printing scores to
    // stdout
    while (!reads_parser.done) {
        entry = reads_parser.next_sequence();
        counts = count_kmers_in_read(entry.read, hapA_kmers, hapB_kmers, k);

        hapA_score = counts.hapA_count * scaling_factor_A;
        hapB_score = counts.hapB_count * scaling_factor_B;

        if (hapA_score > hapB_score) {
            best_haplotype = 'A';
        } else if (hapB_score > hapA_score) {
            best_haplotype = 'B';
        } else {
            best_haplotype = '-';
        }
        printf("%s\t%c\t%.2f\t%.2f\n", entry.id.c_str(), best_haplotype,
                hapA_score, hapB_score);
    }

    // this is test code
    FastaParser faparser("test.fa");
    while (!faparser.done) {
        entry = faparser.next_sequence();
        std::cout << entry.id << '\t' << entry.read << '\n';
    }
}
