#ifndef BIOINF_VITERBI_H
#define BIOINF_VITERBI_H


#include <string>
#include <map>
#include "Sequence.h"
#include "IViterbi.h"

class Viterbi : public IViterbi {

public:
    void alignSequences(Sequence *first, Sequence *second) override;
    Viterbi(const float (&emission_probabilities)[5][5], const float (&trans_prob)[3][3]);
};


#endif //BIOINF_VITERBI_H
