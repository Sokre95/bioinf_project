#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <bits/stdc++.h>

#include "../include/FastaParser.h"
#include "../include/Viterbi.h"
#include "../include/ViterbiLogOdds.h"
#include "../include/ViterbiLogOddsOptimal.h"
#include "../include/MleEstimator.h"

//https://github.com/jarro2783/cxxopts
#include "../include/cxxopts.hpp"


void run_viterbi(std::string filename);
void run_estimator(std::string path);
void load_params(float* averaged_transition_probabilities, float** emission_probabilities);

void print_aligned_sequences_to_console(std::vector<char> top, std::vector<char> bottom, int line_width);
void write_aligned_sequences_to_file(std::string file_path, std::vector<char> top, std::vector<char> bottom, bool single_line, int line_width);
void write_sequence(std::ostream &os, std::vector<char> sequence, bool multiline, int line_width);
void print_sequences(const std::vector<Sequence *> sequences);
void print_duration(long durationS, long durationMs);

cxxopts::Options options("bioinf", "Run either with -v [--viterbi] or -e [--estimate] option. Both options can't be used at same time");

bool print_to_console = false;
bool write_to_file = true;
bool multiline = false;
int line_width = 80;

int main(int argc, char* argv[]) {

    options.add_options("OPTIONS")
            ("v,viterbi", "# Run sequence alignment algorithm on sequence pair given in fasta file specified by <arg> path", cxxopts::value<std::string>())
            ("e,estimate", "# Run HMM parameters estimator with path to directory holding learning database specified by <arg> path. Calculated parameters are stored in ./params.txt file", cxxopts::value<std::string>())
            ("o,out", "# [Use only with -v option] Write aligned sequences to ./aligned/{pair_file_name}.fasta", cxxopts::value<std::string>()->default_value("true"))
            ("c,console" ,"# [Use only with -v option] Print aligned sequences to console")
            ("m,multiline", "# [Use only with -v option] Write/Print aligned sequences in multiple lines, each line containg N chars", cxxopts::value<int>()->implicit_value("100"), "N")
            ("h,help", "# Show help");

    auto result = options.parse(argc, argv);

    if(argc < 1){
        std::cerr << "Wrong number of arguments" << std::endl;
        std::cout << options.help() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (result.count("viterbi") > 0){
        std::string file_path = result["viterbi"].as<std::string>();
        if(result.count("out") > 0){
           if (result["out"].as<std::string>() == "false"){
               write_to_file = false;
           }
        }
        if(result.count("console") > 0){
            print_to_console = true;
        }
        if(result.count("multiline") > 0){
            line_width = result["multiline"].as<int>();
            multiline = true;
        }
        run_viterbi(file_path);
    }
    else if(result.count("estimate") > 0){
        std::string directory_path = result["estimate"].as<std::string>();
        run_estimator(directory_path);
    }
    else{
        std::cout << options.help() << std::endl;
    }
    return 0;
}

void run_estimator(std::string path) {
    std::cout << "Running estimator, please wait..." << std::endl;
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
    FastaParser parser(file_path);
    const std::vector<Sequence *> sequences = parser.parse();

#ifdef DEBUG
    print_sequences(sequences);
#endif

    float* averaged_transition_probabilities = new float[3];
    float** emission_probabilities = new float*[5] ;
    std::cout << "Loading params from params.txt..." << std::endl;

    load_params(averaged_transition_probabilities, emission_probabilities);
    std::map <char, int> lookup_copy(MleEstimator::lookupTable);

    std::cout << "Running Viterbi algorithm. Please wait..." << std::endl;
    auto *logOdds = new ViterbiLogOddsOptimal(averaged_transition_probabilities,
                                       emission_probabilities,
                                       lookup_copy, 0.01, true);

    std::vector<char> top;
    std::vector<char> bottom;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    logOdds->alignSequences(sequences.at(0), sequences.at(1), &top, &bottom);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    long durationS = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    long durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

    if(write_to_file){
        write_aligned_sequences_to_file(file_path, top, bottom, multiline, line_width);
    }
    if(print_to_console){
        print_aligned_sequences_to_console(top, bottom, line_width);
    }

    print_duration(durationS, durationMs);
}

void load_params(float* averaged_transition_probabilities, float** emission_probabilities) {
    std::ifstream params("./params.txt");
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

void print_duration(long durationS, long durationMs) {
    std::cout << "Sequence alignment finished." << std::endl;
    std::cout << "Duration: " << durationS << " s" << std::endl;
    std::cout << "Duration: " << durationMs << " ms" << std::endl;
}


void write_aligned_sequences_to_file(std::string file_path, std::vector<char> top, std::vector<char> bottom, bool multiline, int line_width){

    std::size_t last_slash_index = file_path.find_last_of("/");
    std::string file_name = file_path.substr(last_slash_index +1);
    std::string directory_path = file_path.substr(0, last_slash_index) + "/aligned/";

    system( (std::string("mkdir -p ") + directory_path).c_str() );

    std::ofstream aligned(directory_path + file_name, std::ofstream::out);

    if (!aligned.good()){
        std::cerr << "Error occurred while opening output file: " + directory_path + file_name << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::cout << "Writing aligned sequences to file: " + directory_path + file_name << std::endl;

    aligned << ">" << file_name + "_top" << std::endl;
    write_sequence(aligned, top, multiline, line_width);
    aligned << std::endl;

    aligned << ">" << file_name + "_bottom" << std::endl;
    write_sequence(aligned, bottom, multiline, line_width);
    aligned.close();
}

void write_sequence(std::ostream &os, std::vector<char> sequence, bool multiline, int line_width){
    int curret_line_width = 0;
    for(auto const& c: sequence) {
        os << c << std::flush;
        if (multiline == true){
            curret_line_width++;
            if (curret_line_width == line_width){
                os << std::endl;
                curret_line_width = 0;
            }
        }
    }
}

void print_aligned_sequences_to_console(
        std::vector<char> top,
        std::vector<char> bottom,
        int line_width = 80)
{
    std::cout << "Aligned sequences:" << std::endl;
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
}
