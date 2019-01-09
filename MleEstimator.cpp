//
// Created by tomo on 1/8/19.
//

#include "MleEstimator.h"
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <cstring>
#include <vector>


MleEstimator::MleEstimator(char *directory_path) : _directoryPath(directory_path) {

    this->_lookupTable.insert(std::pair<char, int>('A', 0));
    this->_lookupTable.insert(std::pair<char, int>('C', 1));
    this->_lookupTable.insert(std::pair<char, int>('G', 2));
    this->_lookupTable.insert(std::pair<char, int>('T', 3));
    this->_lookupTable.insert(std::pair<char, int>('-', 4));


    this->_lookupTable.insert(std::pair<char, int>('B', 0));
    this->_lookupTable.insert(std::pair<char, int>('M', 1));
    this->_lookupTable.insert(std::pair<char, int>('X', 3));
    this->_lookupTable.insert(std::pair<char, int>('Y', 4));
    this->_lookupTable.insert(std::pair<char, int>('E', 5));

    memset(_emission_probabilities, 0, sizeof(_emission_probabilities[0][0]) * 5 * 5);
    memset(_transition_probabilities, 0, sizeof(_transition_probabilities[0][0]) * 5 * 5);


}

float **MleEstimator::getEmissionProbabilities() {
    auto **emission = new float *[5];
    for (size_t i = 0; i < 5; ++i) {
        emission[i] = _emission_probabilities[i];
    }
    return emission;
}

float **MleEstimator::getTransitionProbabilities() {
    auto **emission = new float *[5];
    for (size_t i = 0; i < 5; ++i) {
        emission[i] = _transition_probabilities[i];
    }
    return emission;
}

void MleEstimator::estimate() {
    DIR *dir;
    struct dirent *ent;
    struct stat s{};

    if ((dir = opendir(_directoryPath)) != nullptr) {
        /* print all the files and directories within directory */


        unsigned long number_of_transitions = 0;
        unsigned long number_of_pairs = 0;


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

        delete(states);
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
            _emission_probabilities[_lookupTable.at(x)][_lookupTable.at(y)] = float(
                    (emission_frequency + 1.0) / (float) (number_of_pairs + 2.0));
        } else {
            _transition_probabilities[_lookupTable.at(x)][_lookupTable.at(y)] = (float) (
                    (emission_frequency + 1.0) / (float) (number_of_pairs + 2.0));
        }
    }
}



