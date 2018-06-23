#ifndef __TRIO_BINNING_H_INCLUDED__
#define __TRIO_BINNING_H_INCLUDED__

#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>

/**
 * Struct for a single entry of a sequence file, containing
 * an entry id and a read.
 */ // TODO figure out how to do this with polymorphism instead of overloading
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
