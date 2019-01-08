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

    std::string alignSequences(Sequence &first, Sequence &second) override;

    Viterbi(const double *transition_probabilities,const double **emission_probabilities);

private:
    float max(float v1, float v2);
    float max(float v1, float v2, float v3);
};


#endif //BIOINF_VITERBI_H
