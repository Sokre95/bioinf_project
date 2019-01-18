# Pairwise sequence alignment using Hidden Markov Models

Project repository for [Bioinformatics course](http://www.fer.unizg.hr/predmet/bio) at Faculty of Electrical Engineering and Computing, University of Zagreb, 2018./2019.

[![License](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://img.shields.io/packagist/l/doctrine/orm.svg)

## Description
Main goal of this project is implementation of an algorithm for pairwise sequence alignment using HMM.

* [Project documentation](https://github.com/Sokre95/bioinf_project/blob/master/documentation/dokumentacija.pdf)
* [Presentation](https://github.com/Sokre95/bioinf_project/blob/master/documentation/Presentation.pdf)

## Authors
- [Tomislav Božurić](https://github.com/tbozuric)
- [Martin Pisačić ](https://github.com/mpisacic)
- [Krešimir Topolovec](https://github.com/Sokre95)

## Installation 
  Clone this repo and go to `./build` directory
  ```
  git clone git@github.com:Sokre95/bioinf_project.git
  cd bioinf_project/build
  ```
  Build `bioinf` executable
  ```
  cmake ../
  make
  ```
## Usage
Execute `./bioinf --help` to print command info

```
Run either with -v [--viterbi] or -e [--estimate] option. Both options can't be used at the  same time
Usage:
  bioinf [OPTION...]

 OPTIONS options:
  -v, --viterbi arg           # Run sequence alignment algorithm on sequence pair given in fasta file 
                                specified by <arg> path
  -e, --estimate arg          # Run HMM parameters estimator with path to directory holding learning database
                                specified by  <arg> path. Calculated parameters are stored in ./params.txt file
  -o, --out arg               # [Use only with -v option] Write aligned sequences to ./aligned/{pair_file_name}.fasta
                                (default: true)
  -c, --console               # [Use only with -v option] Print aligned sequences to console
  -m, --multiline [=N(=100)]  # [Use only with -v option] Write/Print aligned sequences in multiple lines, each line
                                containg N chars
  -p, --progress              # Show progress while running algorithm
  -t, --mem-optimized         # Run memory optimized version of Viterbi algorithm
  -h, --help                  # Show help
```
## Input
All input files must be in FASTA format and can be placed in arbitrary location. You always have to specify path to the file or folder you are using.

## Examples
Run sequence alignment for `p2.fasta` file writing 200 chars in each line of output file and print output to console also
```
./bioinf --viterbi ../database/pairs/hepatitis/p2.fasta --console --multiline=200 --progress

```
Estimate HMM parameters (probabilities) using `outputs_mafft/upcase` as learning database
```
./bioinf --estimate ../database/outputs_mafft/upcase

```
## Project structure
All c++ files are located in `/src` and `/include` folder. In `/helpers` folder are small python and ruby scripts used as a help in development of this project. 
 
License
---------
MIT License
---------
2018 Tomislav Božurić, Krešimir Topolovec & Martin Pisačić
