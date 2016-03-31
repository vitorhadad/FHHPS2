// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Kernel regression
cx_double kreg(cx_mat Y1, mat X1, rowvec x1, double b1) {

	int n_obs = Y1.n_rows;
	int n_features = X1.n_cols;

	mat xx1 = repmat(x1, n_obs, 1);
	mat psi = (X1 - xx1)/b1;
	mat K = prod(1.0/std::sqrt(2*PI)*arma::exp(-arma::pow(psi,2)/2),1);

	return accu(Y1 % K)/accu(K);
}

// [[Rcpp::export]]
ComplexVector kreg_R(ComplexMatrix Y, NumericMatrix X, NumericMatrix x, double b)   {

    cx_mat Y1 = as<cx_mat>(Y);
    mat X1 = as<mat>(X);    
    rowvec x1 = as<rowvec>(x);
    cx_double output = kreg(Y1, X1, x1, b);
    
    return wrap(output);
}







