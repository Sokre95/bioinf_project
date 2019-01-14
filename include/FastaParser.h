#ifndef BIOINF_FASTAPARSER_H
#define BIOINF_FASTAPARSER_H

#include "Sequence.h"

class FastaParser{

    const std::string filepath;

public:
    /*
     * Creates an instance of fasta parser.
     *
     */
    FastaParser(const std::string f);

    /*
     *Parses the input .fasta file and returns the sequences.
     */
    std::vector<Sequence*> parse();
};

#endif //BIOINF_FASTAPARSER_H
