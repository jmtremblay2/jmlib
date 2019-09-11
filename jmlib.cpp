#include "jmlib.hpp"
#include "types.hpp"

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>


int jmlib::countLines(std::string inputfile){
	int nlines = 0;
	std::ifstream input(inputfile.c_str());
	std::string line;

	while(input){
		std::getline(input,line);
		if(line.compare("") != 0) // make sure we read something
			nlines++;
		}
	return nlines;
	}

int jmlib::countColumns(std::string inputfile){
	std::ifstream input(inputfile.c_str());
	std::string line;
	while(input){
		std::getline(input, line);
		if(line.compare("") != 0){
			int numCol = 0;
			std::stringstream ss(line);
			std::string tok;
			while(ss){
				ss >> tok;
				numCol++;
			}
			return numCol -1;
		}
	}
	return -1;
}

arma::mat jmlib::readMat(std::string inputfile){
	int nrow = countLines(inputfile);
	int ncol = countColumns(inputfile);
	arma::mat M(nrow, ncol);

	std::ifstream input(inputfile.c_str());
	for(int row = 0; row < nrow; ++row)
		for(int col = 0; col < ncol; ++col)
			input >> M(row,col);
	return M;
}

jmlib::vvstring jmlib::getTokens(std::string startsWith, std::string file){
	jmlib::vvstring v;

	std::ifstream input(file.c_str());
	std::string line;
	while(input){
		// Look for a line that start with startsWith
		line = "";
		while(input && line.compare("") == 0)
			std::getline(input,line);
		std::istringstream linestream(line);
		std::string token("");
		while(linestream && token.compare("") == 0)
			linestream >> token;

		// if not, we found a line that starts with something else

		if(token.compare(startsWith) == 0){
			jmlib::vstring c;
			while(linestream){
				std::string var;
				linestream >> var;
				if(var.compare("")!=0)
					c.push_back(var);
				}
			v.push_back(c);
			}
		}
	input.close();
	return v;
	}

/*
void jmlib::removeOptions(std::string model, std::string model_noopt){
	vstring remove;
	remove.push_back("rel_tol");
	remove.push_back("abs_tol");
	remove.push_back("verbose");
	remove.push_back("max_time");
	remove.push_back("max_iter");
	remove.push_back("output");
	remove.push_back("delta");

	std::ifstream in(model.c_str());
	std::ofstream out(model_noopt.c_str());

	std::string line;
	while( ! in.eof()){
	std::getline(in,line);
	bool keep = true;
	for(int i = 0; i < remove.size(); ++i)
		if(line.find(remove.at(i)) != std::string::npos)
			keep = false;

		if(keep)
			out << line << "\n";
		}
	out << "se 0";    
	in.close();
	out.close();
	return;
	}
*/

