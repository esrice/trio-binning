#include "trio_binning.h"
/**
 * This file contains functions for reading sequence files
 * like fasta and fastq. See main() for usage examples.
 */

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

std::ostream& operator<< (std::ostream& stream, const SeqEntry& entry) {
    if (entry.qual == "") {
        stream << ">" << entry.id << "\n" << entry.read << "\n";
    } else {
        stream << "@" << entry.id << "\n" << entry.read << "\n+\n"
            << entry.qual << "\n";
    }
}

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

SeqEntry FastaParser::next_sequence() {
    SeqEntry entry;

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

SeqEntry FastqParser::next_sequence() {
    SeqEntry entry;

    header = line.substr(1, line.length() - 1);
    id = get_id(header);

    std::getline(infile, read); // get the read sequence
    std::getline(infile, line); // get the quality header line
    std::getline(infile, qual); // get the quality string

    std::getline(infile, line); // advance to the next entry
    if (infile.eof()) {
        done = true;
    }

    entry.id = id;
    entry.read = read;
    entry.qual = qual;
    return entry;
}
