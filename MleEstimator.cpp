#include "MleEstimator.h"
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <cstring>
#include<cmath>

const std::map<char, int> MleEstimator::lookupTable = {
        {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'-', 4},
        {'B', 0}, {'M', 1}, {'X', 2}, {'Y', 3}, {'E', 4},
        {'d', 0}, {'t', 1}, {'e', 2}
};

MleEstimator::MleEstimator(char *directory_path) : _directoryPath(directory_path) {
    this->statesTransitionInPairHmm.emplace_back('B', 'M');
    this->statesTransitionInPairHmm.emplace_back('B', 'X');
    this->statesTransitionInPairHmm.emplace_back('B', 'Y');
    this->statesTransitionInPairHmm.emplace_back('B', 'E');
    this->statesTransitionInPairHmm.emplace_back('M', 'M');
    this->statesTransitionInPairHmm.emplace_back('M', 'X');
    this->statesTransitionInPairHmm.emplace_back('M', 'Y');
    this->statesTransitionInPairHmm.emplace_back('M', 'E');
    this->statesTransitionInPairHmm.emplace_back('X', 'X');
    this->statesTransitionInPairHmm.emplace_back('X', 'M');
    this->statesTransitionInPairHmm.emplace_back('X', 'E');
    this->statesTransitionInPairHmm.emplace_back('Y', 'Y');
    this->statesTransitionInPairHmm.emplace_back('Y', 'M');
    this->statesTransitionInPairHmm.emplace_back('Y', 'E');

    memset(_emission_probabilities, 0, sizeof(_emission_probabilities[0][0]) * 5 * 5);
    memset(MleEstimator::_transition_probabilities, 0, sizeof(_transition_probabilities[0][0]) * 4 * 5);
}

float **MleEstimator::getEmissionProbabilities() {
    auto **emission = new float *[5];
    for (size_t i = 0; i < 5; ++i) {
        emission[i] = _emission_probabilities[i];
    }
    return emission;
}

float **MleEstimator::getTransitionProbabilities() {
    auto **emission = new float *[4];
    for (size_t i = 0; i < 4; ++i) {
        emission[i] = _transition_probabilities[i];
    }
    return emission;
}

float *MleEstimator::getAveragedTransitionProbabilities() {
    auto *averagedTransitionProbabilities = new float[3];

    float accumulator = 0.0f;
    accumulator += _transition_probabilities[lookupTable.at('B')][lookupTable.at('X')];
    accumulator += _transition_probabilities[lookupTable.at('B')][lookupTable.at('Y')];
    accumulator += _transition_probabilities[lookupTable.at('M')][lookupTable.at('X')];
    accumulator += _transition_probabilities[lookupTable.at('M')][lookupTable.at('Y')];

    averagedTransitionProbabilities[lookupTable.at('d')] = accumulator / 4.0f;

    accumulator = 0.0f;
    accumulator += _transition_probabilities[lookupTable.at('X')][lookupTable.at('X')];
    accumulator += _transition_probabilities[lookupTable.at('Y')][lookupTable.at('Y')];

    averagedTransitionProbabilities[lookupTable.at('e')] = accumulator / 2.0f;

    accumulator = 0.0f;
    accumulator += _transition_probabilities[lookupTable.at('B')][lookupTable.at('E')];
    accumulator += _transition_probabilities[lookupTable.at('Y')][lookupTable.at('E')];
    accumulator += _transition_probabilities[lookupTable.at('X')][lookupTable.at('E')];
    accumulator += _transition_probabilities[lookupTable.at('M')][lookupTable.at('E')];

    averagedTransitionProbabilities[lookupTable.at('t')] = accumulator / 4.0f;
    return averagedTransitionProbabilities;
}

const std::map<char, int> &MleEstimator::getLookupTable() {
    return MleEstimator::lookupTable;
}

void MleEstimator::estimate() {
    DIR *dir;
    struct dirent *ent;
    struct stat s{};

    if ((dir = opendir(_directoryPath)) != nullptr) {
        /* print all the files and directories within directory */


        unsigned long  long number_of_transitions = 0;
        unsigned long  long number_of_pairs = 0;


        std::map<std::pair<char, char>, unsigned long> match_frequency;
        std::map<std::pair<char, char>, unsigned long> transition_frequency;
        std::map<std::pair<char, char>, unsigned long> emission_x_frequency;
        std::map<std::pair<char, char>, unsigned long> emission_y_frequency;

        auto *states = new std::vector<char>;


        while ((ent = readdir(dir)) != nullptr) {

            auto *file_path = new char[1024];

            strcpy(file_path, _directoryPath);
            strcat(file_path, ent->d_name);
            stat(file_path, &s);

            //it is a file
            if (s.st_mode & S_IFREG) {
                std::ifstream file(file_path);
                std::string id_a;
                std::string sequence_a;
                std::string id_b;
                std::string sequence_b;


                std::getline(file, id_a);
                std::getline(file, sequence_a);
                std::getline(file, id_b);
                std::getline(file, sequence_b);


                states->push_back('B');

                const unsigned long sizeOfSequence = sequence_a.size();

                for (int i = 0; i < sizeOfSequence; i++) {
                    char &a = sequence_a[i];
                    char &b = sequence_b[i];

                    number_of_pairs++;

                    if (a != '-' && b != '-') {
                        states->push_back('M');
                        std::pair<char, char> pair(a, b);
                        increaseFrequency(match_frequency, pair);
                    }

                    if (a != '-' && b == '-') {
                        states->push_back('X');
                        std::pair<char, char> pair(a, b);
                        increaseFrequency(emission_x_frequency, pair);
                    }

                    if (a == '-' && b != '-') {
                        states->push_back('Y');
                        std::pair<char, char> pair(a, b);
                        increaseFrequency(emission_y_frequency, pair);
                    }
                }

                states->push_back('E');
                long sizeOfStates = states->size();
                for (unsigned long i = 0; i < sizeOfStates - 1; i++) {
                    std::pair<char, char> pair(states->at(i), states->at(i + 1));
                    number_of_transitions++;
                    increaseFrequency(transition_frequency, pair);
                }

                states->clear();
            }


        }

        setProbabilities(match_frequency, number_of_pairs, true);
        setProbabilities(emission_x_frequency, number_of_pairs, true);
        setProbabilities(emission_y_frequency, number_of_pairs, true);
        setProbabilities(transition_frequency, number_of_transitions, false);

        delete (states);
        closedir(dir);
    } else {
        perror("Could not open directory");
    }
}

void MleEstimator::increaseFrequency(std::map<std::pair<char, char>, unsigned long> &dictionary,
                                     std::pair<char, char> &pair) {

    if (dictionary.count(pair) == 0) {
        dictionary[{pair.first, pair.second}] = 0;
    }
    dictionary[{pair.first, pair.second}] = dictionary[{pair.first, pair.second}] + 1;
}


void MleEstimator::setProbabilities(std::map<std::pair<char, char>, unsigned long> &dictionary,
                                    unsigned long number_of_pairs, bool emission) {

    for (auto const &item : dictionary) {
        std::pair<char, char> emission_x = item.first;
        unsigned long emission_frequency = item.second;

        char &x = emission_x.first;
        char &y = emission_x.second;
        if (emission) {
            _emission_probabilities[lookupTable.at(x)][lookupTable.at(y)] =
                    (emission_frequency + 1.0f) / (number_of_pairs + 2.0f);
        } else {
            _transition_probabilities[lookupTable.at(x)][lookupTable.at(y)] =
                    (emission_frequency + 1.0f) / (number_of_pairs + 2.0f);
        }
    }

    if (!emission) {
        const double epsilon = 1e-6 /* some small number such as 1e-5 */;
        //za sve prijelaze koji nisu u skupu za ucenje, daj neku malu vjerojatnost
        for (auto &stateTransition :statesTransitionInPairHmm) {
            char &x = stateTransition.first;
            char &y = stateTransition.second;
            //ako je vjerojatnost 0.0
            if (std::abs(_transition_probabilities[lookupTable.at(x)][lookupTable.at(y)] - 0.0f) <=
                epsilon * std::abs(_transition_probabilities[lookupTable.at(x)][lookupTable.at(y)])) {
                _transition_probabilities[lookupTable.at(x)][lookupTable.at(y)] = (1.0f / (number_of_pairs + 2.0f));
            }
        }
    }
}



