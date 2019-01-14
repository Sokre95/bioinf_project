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
    float _transition_probabilities[4][5];
    float _emission_probabilities[5][5];
    char *_directoryPath;
    std::vector<std::pair<char, char>> statesTransitionInPairHmm;
    void increaseFrequency(std::map<std::pair<char, char>, unsigned long> &dictionary, std::pair<char, char> &pair);

    void
    setProbabilities(std::map<std::pair<char, char>, unsigned long> &dictionary, unsigned long number_of_pairs,
                     bool emission);

public:
    /*
     *Creates an instance of Maximum Likelihood Estimator with smoothing parameter.
     */
    MleEstimator(char *directory_path);

    static const std::map<char, int> lookupTable;

    /*
     * Estimates parameters (delta, tau , epsilon) for pairwise HMM.
     */
    void estimate();

    /*
     * Returns the transition probabilities between HMM states.
     */
    float **getTransitionProbabilities();

    /*
     * Returns the emission probabilities.
     */
    float **getEmissionProbabilities();
    /*
     * Returns the averaged transition probabilities, ie the method returns the probabilities used in HMM.
     */
    float *getAveragedTransitionProbabilities();
    /*
     * Returns the lookup table for easier index operation.
     */
    const std::map<char, int> &getLookupTable();
};


#endif //BIOINF_MLEESTIMATOR_H
