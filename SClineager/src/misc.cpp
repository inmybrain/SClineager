//------- Source from SClineager_Gibbs_ver3.cpp: do not edit by hand
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// (190102) source from https://github.com/RcppCore/rcpp-gallery/blob/gh-pages/src/2013-07-13-dmvnorm_arma.Rmd
arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd = false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * std::log(2.0 * M_PI);
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }
  
  if (logd == false) {
    out = exp(out);
  }
  return out;
}

arma::vec diwish_arma(arma::cube x,
                      double df,
                      arma::mat psi,
                      bool logd = false) {
  // Note that constant term 2^{df * xdim / 2} * Gamma(df / 2) is not included!
  int n = x.n_slices;
  int xdim = x.n_cols;
  arma::vec out(n);
  
  double sign;
  double constant = df + xdim + 1;
  for (int i = 0; i < n; i++) {
    double ld_x;
    arma::log_det(ld_x, sign, x.slice(i));
    out(i) = -0.5 * (constant * ld_x + accu(inv_sympd(x.slice(i)) % psi));
  }
  double ld_psi;
  log_det(ld_psi, sign, psi);
  out += 0.5 * df * ld_psi;
  
  if (logd == false) {
    out = exp(out);
  }
  return out;
}


