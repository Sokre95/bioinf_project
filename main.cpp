#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

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
        std::cerr << "\tOption can be either -v (--viterbi), -e (--estimate) or -h (--help)" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return 0;
}

void run_estimator(std::string path) {
    auto *mleEstimator = new MleEstimator("../database/outputs_mafft/upcase/");
    mleEstimator->estimate();

    std::ofstream params("../params.txt");

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
}

void run_viterbi(std::string filename) {
    FastaParser parser("../database/pairs/HIV1_REF_2010/" + filename);
    const std::vector<Sequence *> sequences = parser.parse();

#ifdef DEBUG
    print_sequences(sequences);
#endif

    float* averaged_transition_probabilities = new float[3];
    float** emission_probabilities = new float*[5] ;
    load_params(averaged_transition_probabilities, emission_probabilities);
    std::map <char, int> lookup_copy(MleEstimator::lookupTable);

    auto *logOdds = new ViterbiLogOdds(averaged_transition_probabilities,
                                       emission_probabilities,
                                       lookup_copy, 0.01);

    logOdds->alignSequences(sequences.at(0), sequences.at(1));
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
    std::cout << line << std::endl;
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
    std::cout << "Treba opcije nadopisat" << std::endl;
}

