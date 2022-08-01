#ifndef GESECA_H
#define GESECA_H


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gesecaCpp(const NumericMatrix & E, const NumericVector & inpScores,
                    unsigned genesetSize, unsigned sampleSize, int seed, double eps,
                    double moveScale = 1.0);

#endif // GESECA_H
