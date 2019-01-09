//
// Created by tomo on 1/8/19.
//

#ifndef BIOINF_MLEESTIMATOR_H
#define BIOINF_MLEESTIMATOR_H


#include <string>
#include <map>
#include <vector>

class MleEstimator {
private:
    float _transition_probabilities[5][5];
    float _emission_probabilities[5][5];
    char *_directoryPath;
    std::map<char, int> _lookupTable;
    std::vector<std::pair<char, char>> statesTransitionInPairHmm;

    void increaseFrequency(std::map<std::pair<char, char>, unsigned long> &dictionary, std::pair<char, char> &pair);

    void
    setProbabilities(std::map<std::pair<char, char>, unsigned long> &dictionary, unsigned long number_of_pairs,
                     bool emission);

public:

    MleEstimator(char *directory_path);

    void estimate();

    float **getTransitionProbabilities();

    float **getEmissionProbabilities();

    float *getAveragedTransitionProbabilities();
};


#endif //BIOINF_MLEESTIMATOR_H
