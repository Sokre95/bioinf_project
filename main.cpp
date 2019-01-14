#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <bits/stdc++.h>

#include "FastaParser.h"
#include "Viterbi.h"
#include "ViterbiLogOdds.h"
#include "MleEstimator.h"

#define DEBUG 1

void run_viterbi(std::string filename);
void run_estimator(std::string path);
void load_params(float* averaged_transition_probabilities, float** emission_probabilities);

void print_sequences(const std::vector<Sequence *> sequences);
void print_options();

void print(long duration);

int main(int arg, char* argv[]) {

    if(arg < 2 || arg > 3){
        std::cerr << "Wrong number of arguments";
        print_options();
        std::exit(EXIT_FAILURE);
    }

    std::string option(argv[1]);

    if( option == "-v" || option == "--viterbi"){
        std::string pair_file_path(argv[2]);
        run_viterbi(pair_file_path);
    }
    else if (option == "-e" || option == "--estimate"){
        std::string path(argv[2]);
        run_estimator(path);
    }
    else if (option == "-h" || option == "--help"){
        print_options();
    }
    else {
        std::cerr << "Wrong option of argument:" << std::endl;
        print_options();
        std::exit(EXIT_FAILURE);
    }
    return 0;
}

void run_estimator(std::string path) {
    std::cout << "Running estimator, please wait" << std::endl;
    auto *mleEstimator = new MleEstimator( (char*)path.c_str());
    mleEstimator->estimate();

    std::cout << "Writing to params.txt file:" << std::endl;
    std::ofstream params("./params.txt");

    params << "transitions:" << std::endl;
    for(int i = 0; i < 3; i++){
        params << mleEstimator->getAveragedTransitionProbabilities()[i] << "  " << std::flush;
    }
    params << std::endl;

    params << "emissions:" << std::endl;
    for(int i = 0; i < 5; i++){
        for(int j= 0; j < 5; j++){
            params << mleEstimator->getEmissionProbabilities()[i][j] << "  " << std::flush;
        }
        params << std::endl;
    }

    params.close();
    std::cout << "End" << std::endl;
}

void run_viterbi(std::string file_path) {
    std::cout << "Parsing input file..." << std::endl;
    FastaParser parser(file_path);
    const std::vector<Sequence *> sequences = parser.parse();

#ifdef DEBUG
    print_sequences(sequences);
#endif

    float* averaged_transition_probabilities = new float[3];
    float** emission_probabilities = new float*[5] ;
    std::cout << "Loading params from params.txt ..." << std::endl;

    load_params(averaged_transition_probabilities, emission_probabilities);
    std::map <char, int> lookup_copy(MleEstimator::lookupTable);

    std::cout << "Running Viterbi algorithm..." << std::endl;
    auto *logOdds = new ViterbiLogOdds(averaged_transition_probabilities,
                                       emission_probabilities,
                                       lookup_copy, 0.01);

    std::vector<char> top;
    std::vector<char> bottom;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    logOdds->alignSequences(sequences.at(0), sequences.at(1), &top, &bottom);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    long duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();

    // ispis...

    /* for (auto it = top.rbegin(); it != top.rend(); ++it) {
         std::cout << *it << std::flush;
     }

     std::cout << std::endl;

     for (auto it = bottom.rbegin(); it != bottom.rend(); ++it) {
         std::cout << *it << std::flush;
     } */

    const int line_width = 80;
    ulong curr = top.size() - 1;

    bool isEnd;

    while (curr > 0) {
        int i;
        for (i = 0; i < line_width; i++) {
            if (curr == 0) {
                isEnd = true;
                break;
            }
            curr--;
            std::cout << top.at(curr) << std::flush;
        }

        std::cout << std::endl;

        if (isEnd) {
            curr += i;
        } else {
            curr += line_width - 1;
        }

        for (i = 0; i < line_width; i++) {
            if (curr == 0) break;
            curr--;
            std::cout << bottom.at(curr) << std::flush;
        }

        std::cout << std::endl;
        std::cout << std::endl;
    }

    print(duration);
}

void load_params(float* averaged_transition_probabilities, float** emission_probabilities) {
    std::ifstream params("../params.txt");
    std::string line;
    std::string token;
    int i = 0;

    if (!params.good()){
        std::cerr << "Error occurred while opening params file" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    //check header
    std::getline(params, line);
    if (line.compare("transitions:") != 0){
        std::cerr << "Invalid params file" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Start reading averaged transition probabilities
    std::getline(params, line);
    std::istringstream ss(line);

    for(int i = 0; i < 3; i++) {
        ss >> token;
        averaged_transition_probabilities[i] = std::stof(token);
    };
    // End reading averaged transition probabilities

    // check header
    std::getline(params, line);
    if (line.compare("emissions:") != 0){
        std::cerr << "Invalid params file" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    //Start reading emission probabilities
    i = 0;
    while(std::getline(params, line)){
        std::istringstream ss(line);
        float * arr = new float[5];
        for(int j = 0; j < 5; j++) {
            ss >> token;
            arr[j] = std::stof(token);
        };
        emission_probabilities[i] = arr;
        i++;
    }
    //End reading emission probabilities
}

void print_sequences(const std::vector<Sequence *> sequences) {
    for (int i = 0; i < sequences.size(); ++i) {
        Sequence *s = sequences.at(i);

        std::cout << s->getDescription() << std::endl;

        for (int j = 0; j < s->getSequence().size(); ++j) {
            std::cout << s->getSequence().at(j) << std::flush;
        }
        std::cout << std::endl << std::endl;
    }
}

void print_options() {
    std::cout << "Usage:" << std::endl;
    std::cout << "\tbioinf [options] PATH_TO_FILE_OR_FOLDER:" << std::endl << std::endl;
    std::cout << "Options" << std::endl;
    std::cout << "\t-v [--viterbi] \t\t# Run sequence align algorithm on given pair in fasta file specified by PATH_TO_FILE_OR_FOLDER" << std::endl;
    std::cout << "\t-e [--estimate]\t\t# Run hmm params estimator with path to learning database specified by PATH_TO_FILE_OR_FOLDER" << std::endl;
    std::cout << "\t-h [--help]    \t\t# Show help" << std::endl;
}

void print(long duration) {
    std::cout << "Sequence alignment finished." << std::endl;
    std::cout << "Duration: " << duration << std::endl;
}
