//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_MLEESTIMATOR_H
#define BIOINF_MLEESTIMATOR_H


#include <string>
#include <map>

class MleEstimator {
private:
    float _transition_probabilities[5][5];
    float _emission_probabilities[5][5];
    char *_directoryPath;
    std::map<char, int> _lookupTable;

    void increaseFrequency(std::map<std::pair<char, char>, unsigned long> &dictionary, std::pair<char, char> &pair);

    void
    setProbabilities(std::map<std::pair<char, char>, unsigned long> &dictionary, unsigned long number_of_pairs,
                     bool emission);

public:

    MleEstimator(char *directory_path);

    void estimate();

    float **getTransitionProbabilities();

    float **getEmissionProbabilities();
};


#endif //BIOINF_MLEESTIMATOR_H
