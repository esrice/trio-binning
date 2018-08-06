#include <set>
#include <stdint.h>
#include <getopt.h>
#include <cstdlib>
#include "trio_binning.h"

void print_usage_message()
{
    std::cerr <<
        "--hapA-kmers <a_kmers>:  file of hapA-unique kmers\n"
        "--hapB-kmers <b_kmers>:  file of hapB-unique kmers\n"
        "--input-r1 <reads_file>: fastq containing forward reads to classify\n"
        "--input-r2 <reads_file>: fastq containing reverse reads to classify\n"
        "--hapA-out-r1 <a_out>:   output fastq for forward hapA reads\n"
        "--hapA-out-r2 <a_out>:   output fastq for reverse hapA reads\n"
        "--hapB-out-r1 <a_out>:   output fastq for forward hapB reads\n"
        "--hapB-out-r2 <a_out>:   output fastq for reverse hapB reads\n"
        "--hapU-out-r1 <b_out>:   output fastq for forward hapU reads\n"
        "--hapU-out-r2 <b_out>:   output fastq for reverse hapU reads\n"
        "--help:                  print this message\n";
    exit(1);
}

struct options {
    options() : hapA_kmer_filepath((char *) ""),
                hapB_kmer_filepath((char *) ""),
                forward_inpath((char *) ""),
                reverse_inpath((char *) ""),
                hapA_r1_outpath((char *) ""),
                hapA_r2_outpath((char *) ""),
                hapB_r1_outpath((char *) ""),
                hapB_r2_outpath((char *) ""),
                hapU_r1_outpath((char *) ""),
                hapU_r2_outpath((char *) "") {}

    char* hapA_kmer_filepath;
    char* hapB_kmer_filepath;
    char *forward_inpath, *reverse_inpath;
    char *hapA_r1_outpath, *hapA_r2_outpath;
    char *hapB_r1_outpath, *hapB_r2_outpath;
    char *hapU_r1_outpath, *hapU_r2_outpath;
};

options process_args(int argc, char** argv)
{
    options opts;
    int opt;

    const char* const short_opts = "a:b:i:j:s:t:u:v:w:x:h";
    const option long_opts[] = {
        {"hapA-kmers", required_argument, 0, 'a'},
        {"hapB-kmers", required_argument, 0, 'b'},
        {"input-r1", required_argument, 0, 'i'},
        {"input-r2", required_argument, 0, 'j'},
        {"hapA-out-r1", required_argument, 0, 's'},
        {"hapA-out-r2", required_argument, 0, 't'},
        {"hapB-out-r1", required_argument, 0, 'u'},
        {"hapB-out-r2", required_argument, 0, 'v'},
        {"hapU-out-r1", required_argument, 0, 'w'},
        {"hapU-out-r2", required_argument, 0, 'x'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1)
    {
        switch (opt)
        {
            case 'a': opts.hapA_kmer_filepath = optarg; break;
            case 'b': opts.hapB_kmer_filepath = optarg; break;
            case 'i': opts.forward_inpath = optarg; break;
            case 'j': opts.reverse_inpath = optarg; break;
            case 's': opts.hapA_r1_outpath = optarg; break;
            case 't': opts.hapA_r2_outpath = optarg; break;
            case 'u': opts.hapB_r1_outpath = optarg; break;
            case 'v': opts.hapB_r2_outpath = optarg; break;
            case 'w': opts.hapU_r1_outpath = optarg; break;
            case 'x': opts.hapU_r2_outpath = optarg; break;
            case 'h': print_usage_message(); break;
            case ':':
                std::cerr << "option " << optopt << " requires argument.\n";
                break;
            case '?':
            default:
                std::cerr << "option " << optopt << " is invalid. Ignoring.\n";
        }
    }

    // make sure all required options have been filled
    if (opts.hapA_kmer_filepath[0] == '\0')
    {
        std::cerr << "Must specify --hapA-kmers option.\n";
        print_usage_message();
    }

    else if (opts.hapB_kmer_filepath[0] == '\0')
    {
        std::cerr << "Must specify --hapB-kmers option.\n";
        print_usage_message();
    }

    else if (opts.forward_inpath[0] == '\0')
    {
        std::cerr << "Must specify --input-r1 option.\n";
        print_usage_message();
    }

    else if (opts.reverse_inpath[0] == '\0')
    {
        std::cerr << "Must specify --input-r2 option.\n";
        print_usage_message();
    }

    else if (opts.hapA_r1_outpath[0] == '\0')
    {
        std::cerr << "Must specify --hapA-out-r1 option.\n";
        print_usage_message();
    }

    else if (opts.hapA_r2_outpath[0] == '\0')
    {
        std::cerr << "Must specify --hapA-out-r2 option.\n";
        print_usage_message();
    }

    else if (opts.hapB_r1_outpath[0] == '\0')
    {
        std::cerr << "Must specify --hapB-out-r1 option.\n";
        print_usage_message();
    }

    else if (opts.hapB_r2_outpath[0] == '\0')
    {
        std::cerr << "Must specify --hapB-out-r2 option.\n";
        print_usage_message();
    }

    else if (opts.hapU_r1_outpath[0] == '\0')
    {
        std::cerr << "Must specify --hapU-out-r1 option.\n";
        print_usage_message();
    }

    else if (opts.hapU_r2_outpath[0] == '\0')
    {
        std::cerr << "Must specify --hapU-out-r2 option.\n";
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
    std::ofstream hapA_r1_out, hapA_r2_out, hapB_r1_out, hapB_r2_out;
    std::ofstream hapU_r1_out, hapU_r2_out;

    // counts of number of k-mers unique to each haplotype
    unsigned int num_hapA_kmers, num_hapB_kmers, max_num_kmers;

    // scaling factors for scores, to account for the difference in size
    // between the two sets of k-mers
    double scaling_factor_A, scaling_factor_B;

    // scores for both haplotypes for a read
    int hapA_sum, hapB_sum;
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
    hapA_r1_out.open(opts.hapA_r1_outpath, std::ofstream::out);
    hapA_r2_out.open(opts.hapA_r2_outpath, std::ofstream::out);
    hapB_r1_out.open(opts.hapB_r1_outpath, std::ofstream::out);
    hapB_r2_out.open(opts.hapB_r2_outpath, std::ofstream::out);
    hapU_r1_out.open(opts.hapU_r1_outpath, std::ofstream::out);
    hapU_r2_out.open(opts.hapU_r2_outpath, std::ofstream::out);

    // go through reads
    // TODO add gzip support
    FastqParser r1_parser(opts.forward_inpath), r2_parser(opts.forward_inpath);
    SeqEntry r1_entry, r2_entry;
    haplotype_counts_t r1_counts, r2_counts;

    while (!r1_parser.done) {
        // TODO catch errors caused by r1 and r2 having different # of reads
        r1_entry = r1_parser.next_sequence();
        r2_entry = r2_parser.next_sequence();

        r1_counts = count_kmers_in_read(r1_entry.read,
                hapA_kmers, hapB_kmers, k);
        r2_counts = count_kmers_in_read(r2_entry.read,
                hapA_kmers, hapB_kmers, k);

        hapA_sum = r1_counts.hapA_count + r2_counts.hapA_count;
        hapA_score = hapA_sum * scaling_factor_A;
        hapB_sum = r1_counts.hapB_count + r2_counts.hapB_count;
        hapB_score = hapB_sum * scaling_factor_B;

        if (hapA_score > hapB_score) {
            best_haplotype = 'A';
            hapA_r1_out << r1_entry;
            hapA_r2_out << r2_entry;
        } else if (hapB_score > hapA_score) {
            best_haplotype = 'B';
            hapB_r1_out << r1_entry;
            hapB_r2_out << r2_entry;
        } else {
            best_haplotype = '-';
            hapU_r1_out << r1_entry;
            hapU_r2_out << r2_entry;
        }

        printf("%s\t%c\t%.2f\t%.2f\n", r1_entry.id.c_str(), best_haplotype,
                hapA_score, hapB_score);
    }
}
