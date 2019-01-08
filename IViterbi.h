//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_IVITERBI_H
#define BIOINF_IVITERBI_H


#include <string>
#include "Sequence.h"

class IViterbi {

protected:
    const double *transition_probabilities;
    const double **emission_probabilities;
    std::map<char, int> _lookupTable;

public:
    // pure virtual function providing interface framework.
    virtual std::string alignSequences(Sequence &first, Sequence &second) = 0;

    explicit
    IViterbi(const double *transition_probabilities,const  double **emission_probabilities) : transition_probabilities(
            transition_probabilities), emission_probabilities(emission_probabilities) {

        this->_lookupTable.insert(std::pair<char, int>('A', 0));
        this->_lookupTable.insert(std::pair<char, int>('C', 1));
        this->_lookupTable.insert(std::pair<char, int>('G', 2));
        this->_lookupTable.insert(std::pair<char, int>('T', 3));
        this->_lookupTable.insert(std::pair<char, int>('-', 4));


        this->_lookupTable.insert(std::pair<char, int>('d', 0));
        this->_lookupTable.insert(std::pair<char, int>('t', 1));
        this->_lookupTable.insert(std::pair<char, int>('e', 2));
    }

};


#endif //BIOINF_IVITERBI_H
