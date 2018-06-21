#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iterator>

/**
 * This file contains functions for reading sequence files
 * like fasta and fastq. See main() for usage examples.
 */

/**
 * Struct for a single entry of a sequence file, containing
 * an entry id and a read.
 */
struct seq_entry_t { std::string id, read; };

/*
 * Given a sequence header string, return the sequence id,
 * i.e., everything before the first space.
 */
std::string get_id(std::string header) {
    size_t first_space = header.find(' ');
    if (first_space == std::string::npos) {
        return header;
    } else {
        return header.substr(0, first_space);
    }
}

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
 * Create a new SeqParser by opening the given file,
 * reading the first line, and making sure it's not
 * an empty file.
 */
SeqParser::SeqParser (const char* infile_path) : done(false) {
    infile.open(infile_path);
    if (!infile) {
        std::cerr << "Can't open " << infile_path << "\n";
        exit(1);
    }
    if(!getline(infile, line)) {
        std::cerr << "File " << infile_path << " is empty.\n";
    }
}

/*
 * Implementation of SeqParser for parsing fasta files.
 */
class FastaParser : public SeqParser {
    public:
        FastaParser(const char* infile_path) : SeqParser(infile_path) {};
        seq_entry_t next_sequence();
};

seq_entry_t FastaParser::next_sequence() {
    seq_entry_t entry;

    header = line.substr(1, line.length() - 1);
    id = get_id(header);
    read = "";

    for (std::getline(infile, line); line[0] != '>' && !infile.eof();
            std::getline(infile, line)) {
        read += line;
    }

    if (infile.eof()) {
        done = true;
    }
    entry.id = id;
    entry.read = read;
    return entry;
}

int main() {
    seq_entry_t entry;
    FastaParser parser("test.fa");
    while (!parser.done) {
        entry = parser.next_sequence();
        std::cout << entry.id << '\t' << entry.read << '\n';
    }
}
