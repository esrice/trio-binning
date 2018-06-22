#ifndef __TRIO_BINNING_H_INCLUDED__
#define __TRIO_BINNING_H_INCLUDED__

#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>

/**
 * Struct for a single entry of a sequence file, containing
 * an entry id and a read.
 */
struct seq_entry_t { std::string id, read; };

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
        virtual seq_entry_t next_sequence() = 0;
    protected:
        std::ifstream infile;
        std::string line, header, id, read;
        SeqParser (const char*);
};

/*
 * Implementation of SeqParser for parsing fasta files.
 */
class FastaParser : public SeqParser {
    public:
        FastaParser(const char* infile_path) : SeqParser(infile_path) {};
        seq_entry_t next_sequence();
};

/*
 * Implementation of SeqParser for parsing fastq files.
 */
class FastqParser : public SeqParser {
    public:
        FastqParser(const char* infile_path) : SeqParser(infile_path) {};
        seq_entry_t next_sequence();
    private:
        std::string qual;
};
#endif
