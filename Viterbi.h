//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_VITERBI_H
#define BIOINF_VITERBI_H


#include <string>
#include <map>
#include "Sequence.h"

class Viterbi {
    std::map<char, int> _lookupTable;
    double *transition_probabilities;
    double **emission_probabilities;

public:
    std::string alignSequences(Sequence &first, Sequence &second);

    Viterbi(double *transition_probabilities, double **emission_probabilities);
};


#endif //BIOINF_VITERBI_H
