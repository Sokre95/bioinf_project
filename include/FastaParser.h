#ifndef BIOINF_FASTAPARSER_H
#define BIOINF_FASTAPARSER_H

#include "Sequence.h"

class FastaParser{

    const std::string filepath;

public:
    FastaParser(const std::string f);

    std::vector<Sequence*> parse();
};

#endif //BIOINF_FASTAPARSER_H
