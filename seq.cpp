#include "trio_binning.h"

/**
 * This file contains functions for reading sequence files
 * like fasta and fastq. See main() for usage examples.
 */

/*
 *
 */
std::unique_ptr<std::ostream> ostream_gz_or_uncompressed(char* filepath)
{
    /* outstream_p points to either a std::ofstream or an ogzstream, depending
       on the extension of filepath. */
    std::unique_ptr<std::ostream> outstream_p;
   
    std::string filepath_str(filepath);
    if (filepath_str.substr(filepath_str.length() - 3, 3) == ".gz") {
        outstream_p = std::unique_ptr<std::ostream>(new ogzstream(filepath));
    } else {
        outstream_p = std::unique_ptr<std::ostream>
            (new std::ofstream(filepath));
    }

    return outstream_p;
}

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
    std::string infile_str(infile_path);
    if (infile_str.substr(infile_str.length() - 3, 3) == ".gz") {
        std::cerr << "opening gz file..\n";
        igz.open(infile_path);
        infile_p = &igz;
    } else {
        std::cerr << "opening uncompressed file..\n";
        ifs.open(infile_path);
        infile_p = &ifs;
    }

    if (! *infile_p) {
        std::cerr << "Can't open " << infile_path << "\n";
        exit(1);
    }

    std::getline(*infile_p, line);
    if(infile_p->eof()) {
        std::cerr << "File " << infile_path << " is empty.\n";
        exit(1);
    }
}

SeqEntry FastaParser::next_sequence() {
    SeqEntry entry;

    header = line.substr(1, line.length() - 1);
    id = get_id(header);
    read = "";

    for (std::getline(*infile_p, line); line[0] != '>' && !infile_p->eof();
            std::getline(*infile_p, line)) {
        read += line;
    }

    if (infile_p->eof()) {
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

    std::getline(*infile_p, read); // get the read sequence
    std::getline(*infile_p, line); // get the quality header line
    std::getline(*infile_p, qual); // get the quality string

    std::getline(*infile_p, line); // advance to the next entry
    if (infile_p->eof()) {
        done = true;
    }

    entry.id = id;
    entry.read = read;
    entry.qual = qual;
    return entry;
}
