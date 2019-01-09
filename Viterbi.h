#ifndef BIOINF_VITERBI_H
#define BIOINF_VITERBI_H


#include <string>
#include <map>
#include "Sequence.h"
#include "IViterbi.h"

class Viterbi : public IViterbi {

public:
    void alignSequences(Sequence *first, Sequence *second) override;

    Viterbi(const float *transition_probabilities,
            float **emission_probabilities,
            std::map<char, int> &lookup);
};


#endif //BIOINF_VITERBI_H
