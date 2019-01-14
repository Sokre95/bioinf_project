#include <fstream>
#include <iostream>

#include "../include/FastaParser.h"

FastaParser::FastaParser(const std::string f) : filepath(f){
    std::cout << "Filepath: " << filepath << std::endl;
}

std::vector<Sequence*> FastaParser::parse() {

    std::ifstream file(filepath);
    std::string line;

    if (!file.good()){
        std::cout << "Error occurred while opening input file" << std::endl;
        std::exit(1);
    }

    Sequence *s = new Sequence();

    std::vector<char> sequence;
    std::vector<Sequence*> res;

    while(std::getline(file, line)){
        if (line.empty()){
            continue;
        }
        if(line[0] == '>'){
            if(!sequence.empty()){
                s->setSequence(sequence);
                res.push_back(s);
                sequence.clear();
                s = new Sequence();
            }
            s->setDescription(line.substr(1));
        }
        else{
            for (int i = 0; i < line.size(); ++i) {
                sequence.push_back(line[i]);
            }
        }
    }

    s->setSequence(sequence);
    res.push_back(s);

    return res;
}