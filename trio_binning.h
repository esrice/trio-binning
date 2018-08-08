#ifndef __TRIO_BINNING_H_INCLUDED__
#define __TRIO_BINNING_H_INCLUDED__

#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <set>
#include <cstdlib>
#include <memory>
#include <zlib.h>
#include <gzstream.h>

//---------------------------------------------//
// definitions for things defined in kmer.cpp. //
//---------------------------------------------//

struct haplotype_counts_t {
    int hapA_count, hapB_count;
};

/**
 * Convert a kmer in string form to a bunch of bits, where
 * A=0, C=1, G=2, T=3. For example, ACCGT = 11 10 01 01 00 (little-endian)
 */
uint64_t kmer_to_bits(const std::string& kmer_string);

/**
 * Given a bit representation of a string created by the
 * kmer_to_bits function, convert it back to a string.
 *
 * Arguments:
 * uint64_t bit_repr -- the bit representation of a k-mer
 * size_t k -- the length of the k-mer
 *
 * Returns: std::string of length k containing the k-mer
 */
std::string bits_to_kmer(uint64_t bit_repr, size_t k);

/**
 * Given a DNA sequence, returns the reverse complement.
 */
std::string reverse_complement(std::string sequence);

/**
 * Given a k-mer string, finds the reverse complement, and
 * then returns whichever appears first in alphabetical
 * order, the k-mer or its reverse complement.
 */
std::string get_canonical_representation(std::string sequence);

/**
 * Figures out what size of k-mer is being used by taking
 * a peek at the first line of a file containing a list of
 * unique k-mers.
 */
size_t get_kmer_size(char* file_path);

/**
 * Given a read sequence and two sets of k-mers, counts the
 * number of k-mers from each set that appears in the read
 * and returns a haplotype_counts_t containing these counts.
 *
 * Arguments:
 * std::string read -- the sequence of a read
 * std::set<uint64_t>& hapA_kmers
 *     -- a set of bit representations of k-mers in haplotype A
 * std::set<uint64_t>& hapB_kmers -- same, but for haplotype B
 * size_t k -- k-mer size
 *
 * Returns: haplotype_counts_t instance where
 *          counts.hapA_count is the number of k-mers from
 *          haplotype A that appear in the read and
 *          counts.hapB_count is the number of k-mers from
 *          haplotype B that appear in the read.
 */
haplotype_counts_t count_kmers_in_read(std::string read,
        std::set<uint64_t>& hapA_kmers, std::set<uint64_t>& hapB_kmers,
        size_t k);

/**
 * Reads a file containing a list of k-mers, one per line,
 * into a set in which the k-mers are stored in bit
 * representations.
 *
 * Arguments:
 * char* file_path
 *     -- path to a file containing a list of k-mers, one
 *        per line
 *
 * Returns: std::set<uint64_t>, a set of 64-bit
 *          representations of the k-mers in the given file
 *
 */
std::set<uint64_t> read_kmers_into_set(char* file_path);

//-----------------------------------//
// definitions for things in seq.cpp //
//-----------------------------------//

/*
 * Opens an output stream that is either gzip-compressed or
 * uncompressed, depending on the extension.
 *
 * Argument: char* filepath -- path to file to write
 * Returns: std::unique_ptr<std::ostream>, a smart pointer
 *          to an ostream into that file. To write to it,
 *          just dereference, e.g.
 *              *outstream_p << "text" '
 *          Because it's a smart pointer, it will close
 *          automatically once the program is done with it.
 */
std::unique_ptr<std::ostream> ostream_gz_or_uncompressed(char* filepath);

/**
 * Struct for a single entry of a sequence file, containing
 * an entry id and a read.
 */ // TODO do this with polymorphism instead of overloading
class SeqEntry {
    public:
        std::string id, read, qual;
        SeqEntry (): id(""), read(""), qual("") {};
        SeqEntry (std::string id, std::string read):
            id(id), read(read), qual("") {};
        SeqEntry (std::string id, std::string read, std::string qual):
            id(id), read(read), qual(qual) {};
        friend std::ostream& operator<< (std::ostream& stream,
                const SeqEntry& entry);
};

/*
 * Given a sequence header string, return the sequence id,
 * i.e., everything before the first space.
 */
std::string get_id(std::string header);

/*
 * Abstract class for parsing a sequence file. Contains
 * code for opening and reading the first line of a file,
 * but the next_sequence() method must be implemented by
 * child class.
 */
class SeqParser {
    public:
        bool done; // have we reached EOF?
        virtual SeqEntry next_sequence() = 0;
    protected:
        /* this is a pointer because our input stream could either be an
           ifstream or an igzstream */
        std::istream* infile_p;

        // we will use *either* ifs or igz depending on whether file is gzipped
        std::ifstream ifs;
        igzstream igz;

        std::string line, header, id, read;
        SeqParser (const char*);
};

/*
 * Implementation of SeqParser for parsing fasta files.
 */
class FastaParser : public SeqParser {
    public:
        FastaParser(const char* infile_path) : SeqParser(infile_path) {};
        SeqEntry next_sequence();
};

/*
 * Implementation of SeqParser for parsing fastq files.
 */
class FastqParser : public SeqParser {
    public:
        FastqParser(const char* infile_path) : SeqParser(infile_path) {};
        SeqEntry next_sequence();
    private:
        std::string qual;
};
#endif
