#include "data.hpp"
#include "jmlib.hpp"
#include "types.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <armadillo>

rentix::Data::Data(std::string inputfile){
	std::ifstream input(inputfile.c_str());
	std::string all_names;
	std::getline(input,all_names);
	std::istringstream names_input(all_names);

	// read the first line of the file and extract
	// the variable names it contains
	int index = 0;
	while(names_input){
		std::string s;
		names_input >> s;
		if(s.compare("") != 0){
			namesRaw[s] = index++;
			indexRaw.push_back(s);
			}
		}

	// get sizes
	nobs = jmlib::countLines(inputfile) -1; // -1 for first line with var names
	nvarRaw = namesRaw.size();

	// read the data set
	data = arma::mat(nobs,nvarRaw);
	for(int i = 0; i < nobs; ++i)
		for(int j = 0; j < nvarRaw; ++j)
			input >> data(i,j);
	
	input.close();
	}

void rentix::Data::print(int from, int to) const {;
	for(int i = from; i < to; ++i)
		std::cout << data.row(i);
	return;
	}

arma::rowvec rentix::Data::getRow(int i) const {
	return data.row(i);
	}

arma::colvec rentix::Data::getVarRaw(int i) const {
	return data.col(i);
	}

arma::colvec rentix::Data::getVar(std::string varName) const {
	// namesRaw.find(varNames) returns a pair ("VAR_NAME",ID)
	// we need the ID, accessed with second
	return data.col(namesRaw.find(varName)->second);
	}

double rentix::Data::getValRaw(int rowNum, int colNum) const {
	return data(rowNum,colNum);
	}

double rentix::Data::operator()(int rowNum, const std::string& varName) const {
	return data(rowNum,namesRaw.find(varName)->second);
	}
/*
double& rentix::Data::operator()(int rowNum, const std::string& varName){
	return data(rowNum,namesRaw.find(varName)->second);
	}
*/
std::string rentix::Data::getVarNameRaw(int index) const {
	return indexRaw.at(index);
	}

int rentix::Data::getVarIndexRaw(std::string varName) const {
	return namesRaw.find(varName)->second;
	}

int rentix::Data::getNobs() const {
	return nobs;
	}

int rentix::Data::getNvarRaw() const {
	return nvarRaw;
	}
