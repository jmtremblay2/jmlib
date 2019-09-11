#ifndef JMLIB
#define JMLIB

#include <string>
#include <armadillo>
#include <iostream>

#include "types.hpp"

#include <armadillo>

namespace jmlib{

/**
 * counts the number of non-empty lines in a file
 *
 * @param inputfile the name and path of the file
 *
 * @return the number of lines that are not empty
 */
int countLines(std::string inputfile);

int countColumns(std::string inputfile);

arma::mat readMat(std::string inputfile);
/**
 * get a the list of words on lines that start with startsWith
 *
 * parses the file and get the words on lines that start with startsWith
 * the words are stored in a vector<vector<string> > each element of the
 * first vector is a line, and each element of the second vector
 * is a word
 *
 * @param startsWith the token that lines start with
 * @param file the name of the file
 *
 * @return a vector that contains all the words in one line
 */
vvstring getTokens(std::string startsWith, std::string file);

/**
 * get a value from a file
 *
 * getVal will look for a line that will look like :
 * startsWith value
 * and will return the value. If the value is not
 * specified the def value will be returned.
 *
 * @param startsWith the identifier of the variable
 * @param file the name of the file to look into
 * @param def a default value to return
 *
 * @return the value as read in the data file (or def)
 */
template<class T>
T getVal(std::string startswith, std::string file, T def){
    vvstring fetch = getTokens(startswith,file);
    if(0 == fetch.size())
        return def;
    std::stringstream ss(fetch.at(0).at(0));
    T val;
    ss >> val;
    return val;
    }
} // namespace jmlib

#endif
