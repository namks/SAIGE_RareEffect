// [[Rcpp::depends(BH)]]

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <sstream>
#include <time.h>
#include <Rcpp.h>
#include <stdint.h>
#include <boost/math/distributions/normal.hpp>
#include "SPA_binary.hpp"
#include "SPA_survival.hpp"
#include "UTIL.hpp"

// [[Rcpp::export]]
void SPA(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, double & pval, bool & isSPAConverge){
        using namespace Rcpp;
        List result;
        double p1, p2;
        bool Isconverge = true;
	Rcpp::List outuni1;
	Rcpp::List outuni2;
        if( traitType == "binary"){
          outuni1 = getroot_K1_Binom(0, mu, g, q, tol);
          outuni2 = getroot_K1_Binom(0, mu, g, qinv, tol);
        }else if(traitType == "survival"){
          outuni1 = getroot_K1_Poi(0, mu, g, q, tol);
          outuni2 = getroot_K1_Poi(0, mu, g, qinv, tol);
        }

/*
        double outuni1root = outuni1["root"];
        double outuni2root = outuni2["root"];
        bool Isconverge1 = outuni1["Isconverge"];
        bool Isconverge2 = outuni2["Isconverge"];

        std::cout << "outuni1root" << outuni1root << std::endl;
        std::cout << "outuni2root" << outuni2root << std::endl;
        std::cout << "Isconverge1" << Isconverge1 << std::endl;
        std::cout << "Isconverge2" << Isconverge2 << std::endl;
*/

        Rcpp::List getSaddle;
        Rcpp::List getSaddle2;
        if(outuni1["Isconverge"]  && outuni2["Isconverge"])
        {
                if( traitType == "binary"){
                  getSaddle = Get_Saddle_Prob_Binom(outuni1["root"], mu, g, q, logp);
                  getSaddle2 = Get_Saddle_Prob_Binom(outuni2["root"], mu, g, qinv, logp);
                }else if(traitType == "survival"){
		/*
		std::cout << "q " << q << std::endl;
		std::cout << "qinv " << qinv << std::endl;
		mu.print("mu");
		g.print("g");
		*/
                  getSaddle = Get_Saddle_Prob_Poi(outuni1["root"], mu, g, q, logp);
                  getSaddle2 = Get_Saddle_Prob_Poi(outuni2["root"], mu, g, qinv, logp);
		/*bool isSaddleprint = getSaddle["isSaddle"];
		bool isSaddle2print = getSaddle2["isSaddle"];
		std::cout << "isSaddle " << isSaddleprint << std::endl;
		std::cout << "isSaddle2 " << isSaddle2print << std::endl;
		*/
                }


                if(getSaddle["isSaddle"]){
                        p1 = getSaddle["pval"];
        		//std::cout << "p1 " << p1 << std::endl;

                }else{
			Isconverge = false;
                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }
                if(getSaddle2["isSaddle"]){
                        p2 = getSaddle2["pval"];
        		//std::cout << "p2 " << p2 << std::endl;
                }else{
			Isconverge = false;
                        if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
                }

                if(logp)
                {
                        pval = add_logp(p1,p2);
                } else {
                        pval = std::abs(p1)+std::abs(p2);
                }
                //Isconverge=true;
        }else {
                        //std::cout << "Error_Converge" << std::endl;
                        pval = pval_noadj;
                        Isconverge=false;
                }
        isSPAConverge = Isconverge;
	//std::cout << "pval " << pval << std::endl;
	//result["pvalue"] = pval;
        //result["Isconverge"] = Isconverge;
        //return(result);
}


// [[Rcpp::export]]
void SPA_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol, std::string traitType, double & pval, bool & isSPAConverge){

        using namespace Rcpp;
        List result;
        double p1, p2;
        bool Isconverge = true;
	Rcpp::List outuni1;
        Rcpp::List outuni2;
/*	std::cout << "mu.n_elem " << mu.n_elem << std::endl;
	std::cout << "g.n_elem " << g.n_elem << std::endl;
	std::cout << "gNA.n_elem " << gNA.n_elem << std::endl;
	std::cout << "gNB.n_elem " << gNB.n_elem << std::endl;
	std::cout << "muNA.n_elem " << muNA.n_elem << std::endl;
	std::cout << "muNB.n_elem " << muNB.n_elem << std::endl;
*/
        if( traitType == "binary"){
	    //   arma::vec timeoutput1 = getTime();
          outuni1 = getroot_K1_fast_Binom(0, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
	  //   arma::vec timeoutput2 = getTime();
 //printTime(timeoutput1, timeoutput2, "getroot_K1_fast_Binom");
          //double qinv = -1*q;
          outuni2 = getroot_K1_fast_Binom(0, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
	 // arma::vec timeoutput3 = getTime();
 	//printTime(timeoutput2, timeoutput3, "getroot_K1_fast_Binom");
        }else if(traitType == "survival"){
          outuni1 = getroot_K1_fast_Poi(0, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
          outuni2 = getroot_K1_fast_Poi(0, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
        }

/*	double outuni1root = outuni1["root"];
        double outuni2root = outuni2["root"];
        bool Isconverge1 = outuni1["Isconverge"];
        bool Isconverge2 = outuni2["Isconverge"];

        std::cout << "outuni1root" << outuni1root << std::endl;
        std::cout << "outuni2root" << outuni2root << std::endl;
        std::cout << "Isconverge1" << Isconverge1 << std::endl;
        std::cout << "Isconverge2" << Isconverge2 << std::endl;
*/

        Rcpp::List getSaddle;
        Rcpp::List getSaddle2;
        if(outuni1["Isconverge"]  && outuni2["Isconverge"])
        {
          if( traitType == "binary"){
		                //  arma::vec timeoutput4 = getTime();
                getSaddle  = Get_Saddle_Prob_fast_Binom(outuni1["root"], mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
		// arma::vec timeoutput5 = getTime();
                getSaddle2 = Get_Saddle_Prob_fast_Binom(outuni2["root"], mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
	//	arma::vec timeoutput6 = getTime();
	//printTime(timeoutput4, timeoutput5, "Get_Saddle_Prob_fast_Binom");
 //printTime(timeoutput5, timeoutput6, "Get_Saddle_Prob_fast_Binom");
          }else if(traitType == "survival"){
                getSaddle  = Get_Saddle_Prob_fast_Poi(outuni1["root"], mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
                getSaddle2 = Get_Saddle_Prob_fast_Poi(outuni2["root"], mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
          }
                if(getSaddle["isSaddle"]){
                        p1 = getSaddle["pval"];
                }else{
			Isconverge = false;

                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }

                if(getSaddle2["isSaddle"]){
                        p2 = getSaddle2["pval"];

                }else{
			Isconverge = false;
                        if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
                }

                if(logp){
                        pval = add_logp(p1,p2);
                }else {
                        pval = std::abs(p1)+std::abs(p2);
                }
                //Isconverge=true;
        }else {
                        //std::cout << "Error_Converge" << std::endl;
                        pval = pval_noadj;
                        Isconverge=false;
        }
	isSPAConverge = Isconverge;
        //result["pvalue"] = pval;
        //result["Isconverge"] = Isconverge;
        //return(result);
}




// [[Rcpp::export]]
double SPA_pval(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, bool & isSPAConverge){
       double pval;
        using namespace Rcpp;
        List result;
        double p1, p2;
        bool Isconverge = true;
        Rcpp::List outuni1;
        Rcpp::List outuni2;
        if( traitType == "binary"){
          outuni1 = getroot_K1_Binom(0, mu, g, q, tol);
          outuni2 = getroot_K1_Binom(0, mu, g, qinv, tol);
        }else if(traitType == "survival"){
          outuni1 = getroot_K1_Poi(0, mu, g, q, tol);
          outuni2 = getroot_K1_Poi(0, mu, g, qinv, tol);
        }


        double outuni1root = outuni1["root"];
        double outuni2root = outuni2["root"];
        bool Isconverge1 = outuni1["Isconverge"];
        bool Isconverge2 = outuni2["Isconverge"];

/*
        std::cout << "outuni1root" << outuni1root << std::endl;
        std::cout << "outuni2root" << outuni2root << std::endl;
        std::cout << "Isconverge1" << Isconverge1 << std::endl;
        std::cout << "Isconverge2" << Isconverge2 << std::endl;
*/

        Rcpp::List getSaddle;
        Rcpp::List getSaddle2;
        if(outuni1["Isconverge"]  && outuni2["Isconverge"])
        {
                if( traitType == "binary"){
                  getSaddle = Get_Saddle_Prob_Binom(outuni1["root"], mu, g, q, logp);
                  getSaddle2 = Get_Saddle_Prob_Binom(outuni2["root"], mu, g, qinv, logp);
                }else if(traitType == "survival"){
                  getSaddle = Get_Saddle_Prob_Poi(outuni1["root"], mu, g, q, logp);
                  getSaddle2 = Get_Saddle_Prob_Poi(outuni2["root"], mu, g, qinv, logp);
                }

                if(getSaddle["isSaddle"]){
                        p1 = getSaddle["pval"];
                }else{
                        Isconverge = false;
                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }
                if(getSaddle2["isSaddle"]){
                        p2 = getSaddle2["pval"];
                }else{
                        Isconverge = false;
                        if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
                }
		//std::cout << "p1 " << p1 << " p2 " << p2 << std::endl;
                if(logp)
                {
                        pval = add_logp(p1,p2);
                } else {
                        pval = std::abs(p1)+std::abs(p2);
                }
		//std::cout << "pval " << pval << std::endl;
                //Isconverge=true;
        }else {
                        //std::cout << "Error_Converge" << std::endl;
                        pval = pval_noadj;
                        Isconverge=false;
                }
        isSPAConverge = Isconverge;

	return(pval);
        //result["pvalue"] = pval;
        //result["Isconverge"] = Isconverge;
        //return(result);
}

