#include <set>
#include <stdint.h>
#include <getopt.h>
#include <cstdlib>
#include "trio_binning.h"

void print_usage_message()
{
    std::cerr <<
        "--hapA-kmers <a_kmers>:     file of hapA-unique kmers\n"
        "--hapB-kmers <b_kmers>:     file of hapB-unique kmers\n"
        "--input-reads <reads_file>: fasta/q containing reads to classify\n"
        "--hapA-out <a_out>:         output file for hapA reads\n"
        "--hapB-out <b_out>:         output file for hapB reads\n"
        "--hapU-out <u_out>:         output file for unknown hap reads\n"
        "--help:                     print this message\n";
    exit(1);
}

struct options {
    options() : hapA_kmer_filepath((char *) ""),
                hapB_kmer_filepath((char *) ""),
                reads_inpath((char *) ""), hapA_reads_outpath((char *) ""),
                hapB_reads_outpath((char *) ""),
                hapU_reads_outpath((char *) "") {}

    char* hapA_kmer_filepath;
    char* hapB_kmer_filepath;
    char* reads_inpath;
    char* hapA_reads_outpath;
    char* hapB_reads_outpath;
    char* hapU_reads_outpath;
};

options process_args(int argc, char** argv)
{
    options opts;
    int opt;

    const char* const short_opts = "a:b:i:s:t:u:h";
    const option long_opts[] = {
        {"hapA-kmers", required_argument, 0, 'a'},
        {"hapB-kmers", required_argument, 0, 'b'},
        {"input-reads", required_argument, 0, 'i'},
        {"hapA-out", required_argument, 0, 's'},
        {"hapB-out", required_argument, 0, 't'},
        {"hapU-out", required_argument, 0, 'u'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1)
    {
        switch (opt)
        {
            case 'a':
                opts.hapA_kmer_filepath = optarg;
                break;
            case 'b':
                opts.hapB_kmer_filepath = optarg;
                break;
            case 'i':
                opts.reads_inpath = optarg;
                break;
            case 's':
                opts.hapA_reads_outpath = optarg;
                break;
            case 't':
                opts.hapB_reads_outpath = optarg;
                break;
            case 'u':
                opts.hapU_reads_outpath = optarg;
                break;
            case 'h':
                print_usage_message();
                break;
            case ':':
                std::cerr << "option " << optopt << " requires argument.\n";
            case '?':
            default:
                std::cerr << "option " << optopt << " is invalid. Ignoring.\n";
        }
    }

    // make sure all required options have been filled
    if (opts.hapA_kmer_filepath[0] == '\0')
    {
        std::cerr << "Must specify -a option.\n";
        print_usage_message();
    }

    else if (opts.hapB_kmer_filepath[0] == '\0')
    {
        std::cerr << "Must specify -b option.\n";
        print_usage_message();
    }

    else if (opts.reads_inpath[0] == '\0')
    {
        std::cerr << "Must specify -i option.\n";
        print_usage_message();
    }

    else if (opts.hapA_reads_outpath[0] == '\0')
    {
        std::cerr << "Must specify -s option.\n";
        print_usage_message();
    }

    else if (opts.hapB_reads_outpath[0] == '\0')
    {
        std::cerr << "Must specify -t option.\n";
        print_usage_message();
    }

    else if (opts.hapU_reads_outpath[0] == '\0')
    {
        std::cerr << "Must specify -u option.\n";
        print_usage_message();
    }

    return opts;
}

int main(int argc, char** argv)
{
    options opts; // struct to hold command-line arguments

    // sets of hapA or hapB specific k-mers, encoded as 64-bit ints
    std::set<uint64_t> hapA_kmers, hapB_kmers;

    // an iterator for above sets
    std::set<uint64_t>::iterator it;

    // output streams for haplotype-specific reads (A and B) and reads which
    // cannot be assigned to a haplotype (U)
    std::ofstream hapA_reads_out, hapB_reads_out, hapU_reads_out;

    // counts of number of k-mers unique to each haplotype
    unsigned int num_hapA_kmers, num_hapB_kmers, max_num_kmers;

    // scaling factors for scores, to account for the difference in size
    // between the two sets of k-mers
    double scaling_factor_A, scaling_factor_B;

    // scores for both haplotypes for a read
    float hapA_score, hapB_score;
    char best_haplotype;

    // get arguments
    opts = process_args(argc, argv);

    size_t k = get_kmer_size(opts.hapA_kmer_filepath);
    uint8_t *seq;

    int i;

    // read k-kmers into sets
    hapA_kmers = read_kmers_into_set(opts.hapA_kmer_filepath);
    hapB_kmers = read_kmers_into_set(opts.hapB_kmer_filepath);

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

    // set up some output streams for haplotype reads
    hapA_reads_out.open(opts.hapA_reads_outpath, std::ofstream::out);
    hapB_reads_out.open(opts.hapB_reads_outpath, std::ofstream::out);
    hapU_reads_out.open(opts.hapU_reads_outpath, std::ofstream::out);

    // go through reads
    // TODO add gzip support and automatic fasta/fastq detection
    FastqParser reads_parser(opts.reads_inpath);
    SeqEntry entry;
    haplotype_counts_t counts;

    while (!reads_parser.done) {
        entry = reads_parser.next_sequence();
        counts = count_kmers_in_read(entry.read, hapA_kmers, hapB_kmers, k);

        hapA_score = counts.hapA_count * scaling_factor_A;
        hapB_score = counts.hapB_count * scaling_factor_B;

        if (hapA_score > hapB_score) {
            best_haplotype = 'A';
            hapA_reads_out << entry;
        } else if (hapB_score > hapA_score) {
            best_haplotype = 'B';
            hapB_reads_out << entry;
        } else {
            best_haplotype = '-';
            hapU_reads_out << entry;
        }
        printf("%s\t%c\t%.2f\t%.2f\n", entry.id.c_str(), best_haplotype,
                hapA_score, hapB_score);
    }
}
