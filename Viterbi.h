//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_VITERBI_H
#define BIOINF_VITERBI_H


#include <string>
#include <map>
#include "Sequence.h"
#include "IViterbi.h"

class Viterbi : public IViterbi {

public:

    void alignSequences(Sequence &first, Sequence &second) override;

    Viterbi(const float *transition_probabilities, const float **emission_probabilities);
};


#endif //BIOINF_VITERBI_H
