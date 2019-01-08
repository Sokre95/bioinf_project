//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_VITERBILOGODDS_H
#define BIOINF_VITERBILOGODDS_H


#include <map>
#include "IViterbi.h"

class ViterbiLogOdds : public IViterbi {
    double eta;

public:

    std::string alignSequences(Sequence &first, Sequence &second) override;

    ViterbiLogOdds(const double *transition_probabilities,const double **emission_probabilities, const double eta);
};



#endif //BIOINF_VITERBILOGODDS_H
