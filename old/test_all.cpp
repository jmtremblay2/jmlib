#include "test_all.hpp"

#include "data.hpp"
#include "datareg.hpp"

#include "model.hpp"
#include "oprobit.hpp"
#include "logit.hpp"
#include "probit.hpp"
#include "regsim.hpp"
#include "discrete_cont1.hpp"


#include "jmlib.hpp"
#include "jmlib_types.hpp"


#include <armadillo>
#include <iostream>
#include <cmath>

bool rentix::test_data() {
    double tol = 0.0000000001;
    bool success = true;
    rentix::Data D("test_files/data_test.txt");

    // test getRow(int)
    if( ! jmlib::matEq(D.getRow(0), arma::rowvec("1 2 3 4")) ||
            ! jmlib::matEq(D.getRow(1), arma::rowvec("4 3 2 1"))) {
        std::cerr << "problem with rentix::Data::getRow(int)\n";
        success = false;
        }

    // test getVarRaw(int)
    if( ! jmlib::matEq(D.getVarRaw(0), arma::colvec("1 4")) ||
            ! jmlib::matEq(D.getVarRaw(1), arma::colvec("2 3")) ||
            ! jmlib::matEq(D.getVarRaw(2), arma::colvec("3 2")) ||
            ! jmlib::matEq(D.getVarRaw(3), arma::colvec("4 1"))) {
        std::cerr << "problem with rentix::Data::getVarRaw(int)\n";
        success = false;
        }

    // test getVar(string)
    if( ! jmlib::matEq(D.getVar("a"), arma::colvec("1 4")) ||
            ! jmlib::matEq(D.getVar("b"), arma::colvec("2 3")) ||
            ! jmlib::matEq(D.getVar("c"), arma::colvec("3 2")) ||
            ! jmlib::matEq(D.getVar("d"), arma::colvec("4 1"))) {
        std::cerr << "problem with rentix::Data::getVar(string)\n";
        success = false;
        }


    double realVals[][4] = {{1,2,3,4},{4,3,2,1}};
    jmlib::vstring realNames;
    realNames.push_back("a");
    realNames.push_back("b");
    realNames.push_back("c");
    realNames.push_back("d");

    // test getValRaw(int,int)
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 4; ++j)
            if(std::fabs(realVals[i][j] - D.getValRaw(i,j)) > tol) {
                std::cerr << "problem with rentix::Data::getValRaw(int,int)\n";
                success = false;
                }


    // test operator()
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 4; ++j)
            if(std::fabs(realVals[i][j] - D(i,realNames.at(j))) > tol) {
                std::cerr << "problem with rentix::Data::getValRaw(int,int)\n";
                success = false;
                }

    // test getVarNameRaw(int)
    for(int i = 0; i < 4; ++i)
        if(D.getVarNameRaw(i).compare(realNames.at(i)) != 0) {
            std::cerr << "problem with rentix::Data::getVarNameRaw(int)\n";
            success = false;
            }

    // test getVarIndexRaw
    for(int i = 0; i < 4; ++i)
        if(D.getVarIndexRaw(realNames.at(i)) != i) {
            std::cerr << "problem with rentix::Data::getVarIndexRaw(string)\n";
            success = false;
            }

    // test getNobs
    if( 2 != D.getNobs()) {
        std::cerr << "problem with rentix::Data::getNobs()\n";
        success = false;
        }

    // test getNvarRaw
    if ( 4 != D.getNvarRaw()) {
        std::cerr << "problem with rentix::Data::getNvarRaw()\n";
        success = false;
        }

    if(success)
        std::cout << "class Data succesfully checked\n";
    return success;
    }

bool rentix::test_dataReg() {
    DataReg D("test_files/data_test.txt", "test_files/datareg_model.txt");
    bool success = true;

    // test getObs;
    if(! jmlib::matEq(arma::rowvec("1 2"), D.getObs(0)) ||
            ! jmlib::matEq(arma::rowvec("4 3"), D.getObs(1))) {
        std::cerr << "problem with rentix::DataReg::getObs(int)\n";
        success = false;
        }

    // test getXBeta
    arma::colvec beta1("1 1");
    jmlib::vdouble beta2(2);
    beta2.at(0) = beta2.at(1) = 1;

    if( (! jmlib::matEq(arma::colvec("3 7"), D.getXBeta(beta1))) ||
            (! jmlib::matEq(arma::colvec("3 7"), D.getXBeta(beta2)))) {
        std::cerr << "problem with rentix::DataReg::getXBeta(vdouble or colvec)\n";
        success = false;
        }


    double tol = 0.00000001;
    // test getY
    if( ! jmlib::matEq(arma::colvec("4 1"), D.getY())) {
        std::cerr << "problem with rentix::DataReg::getY()\n";
        success = false;
        }
    if( fabs( D.getY(0) - 4) > tol  || fabs( D.getY(1) - 1) > tol ) {
        std::cerr << "problem with rentix::DataReg::getY(int)\n";
        success = false;
        }

    // test getVarName
    if( D.getVarName(0).compare("a") != 0 || D.getVarName(1).compare("b") != 0) {
        std::cerr << "problem with rentix::DataReg::getVarName(int)\n";
        success = false;
        }

    // test getVarIndex
    if( D.getVarIndex("a") != 0 || D.getVarIndex("b") != 1) {
        std::cerr << "problem with rentix::DataReg::getVarIndex(string)\n";
        success = false;
        }

    // test getNvar()
    if( D.getNvar() != 2 ) {
        std::cerr << "problem with rentix::DataReg::getNvar()\n";
        success = false;
        }

    if(success)
        std::cout << "class DataReg succesfully checked\n";
    return success;
    }

bool rentix::test_oprobit() {
    Oprobit op("/home/jm/Dropbox/CPP/rentix/test_files/data_oprobit","/home/jm/Dropbox/CPP/rentix/test_files/model_oprobit");
    jmlib::vvdouble x;
    x.push_back(jmlib::vdouble(16,0));
    x.push_back(jmlib::vdouble(1,1));

    jmlib::vvdouble g;
    g.push_back(jmlib::vdouble(16,1));
    g.push_back(jmlib::vdouble(1,1));

    //op.LL(x,g);
    op.solve();
    if(op.getMaxLL() > -10448 && op.getMaxLL() < -10447){
        std::cout << "Oprobit succesfully checked\n";
        return true;
        }
    else{
        std::cerr << "Problem with Oprobit*****************\n";
        return false;
        }
    }

bool rentix::test_logit() {
	Logit l("/home/jm/Dropbox/CPP/rentix/test_files/logit_example_data","/home/jm/Dropbox/CPP/rentix/test_files/logit_example_model");
	l.solve();
	if(l.getMaxLL() > -7694 && l.getMaxLL() < -7693){
        std::cout << "Logit succesfully checked\n";
        return true;
        }
    else{
        std::cerr << "Problem with Logit*****************\n";
        return false;
        }
    return true;
    }

bool rentix::test_probit(){
    Probit p("/home/jm/Dropbox/CPP/rentix/test_files/probit_example_data","/home/jm/Dropbox/CPP/rentix/test_files/probit_example_model");
    p.solve(1);
    return true;
    }

bool rentix::test_regSim(){
    RegSim rs("/home/jm/Dropbox/CPP/rentix/test_files/regsim_data", "/home/jm/Dropbox/CPP/rentix/test_files/regsim_model");
    //rs.solve(1);
    return true;
    }


bool rentix::test_DC1(){
    DC1 dc("/home/jm/Dropbox/CPP/rentix/test_files/dc1_data","/home/jm/Dropbox/CPP/rentix/test_files/dc1_model");
    dc.solve();
    }


bool rentix::test_rentix() {

    return test_data() && test_dataReg() && test_logit();
    }
