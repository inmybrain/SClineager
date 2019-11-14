//------- Source from SClineager_Gibbs_ver2.cpp: do not edit by hand
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include "misc.h"
double getPostLike(arma::mat k_mat,
                   arma::mat transform_mat,
                   arma::mat genotype_mat,
                   arma::uvec is_high_cover,
                   arma::vec mu,
                   arma::mat sigma,
                   arma::mat psi,
                   int dfreedom,
                   bool logd){
  int g_nr = genotype_mat.n_rows;
  int k_nc = k_mat.n_cols;
  double loglike;
  
  loglike = accu(log(normpdf(transform_mat.elem(is_high_cover), 
                             genotype_mat(is_high_cover),
                             k_mat.elem(is_high_cover))));
    for(int j = 0; j < g_nr; j++) {
      rowvec muvec(k_nc);
      muvec.fill(mu(j));
      loglike += dmvnrm_arma(genotype_mat.row(j),
                             muvec,
                             sigma,
                             logd)(0);
    }
    cube Sigma(k_nc, k_nc, 1);
  Sigma.slice(0) = sigma;
  loglike += diwish_arma(Sigma, dfreedom, psi, logd)(0);
  return loglike;
}

//' @name sclineager_gibbs
//' @title An internal function for Gibbs sampling
//' @description An internal function for Gibbs sampling wrapped by \link{sclineager_internal}
//' @seealso \link{sclineager_internal}
//' @exports
// [[Rcpp::export]]
Rcpp::List sclineager_gibbs(
    arma::mat psi,
    arma::mat k_mat,
    arma::mat transform_mat,
    arma::vec mu,
    int dfreedom,
    arma::mat sigma,
    int max_iter = 100,
    bool save = false
){
  double tick = 0.2;
  
  mat genotype_mat = transform_mat;
  mat genotype_mat_avg = zeros(genotype_mat.n_rows, genotype_mat.n_cols), sigma_avg = zeros(sigma.n_rows, sigma.n_cols);
  cube genotype_mat_all(transform_mat.n_rows, transform_mat.n_cols, max_iter);
  cube sigma_all(sigma.n_rows, sigma.n_cols, max_iter);
  uvec is_low_cover = find_nonfinite(transform_mat);
  uvec is_high_cover = find_finite(transform_mat);
  int g_nr = genotype_mat.n_rows;
  int k_nc = k_mat.n_cols;
  vec loglike(max_iter);
  
  for(int iter = 0; iter < max_iter; iter++) {
    if(any(iter == regspace(1, int(max_iter * tick), max_iter))){
      printf("%3.f%%...", 100 * double(iter) / double(max_iter));
    }
    mat inverse_sigma = inv_sympd(sigma);
    mat scale_mat = psi;
    
    for(int j = 0; j < g_nr; j++) {
      mat tmp = transform_mat;
      tmp.elem(is_low_cover).fill(0);
      mat inverse_Kj = zeros(k_nc, k_nc);
      for(int k = 0; k < k_nc; k++) {
        if(conv_to<double>::from(k_mat.row(j).col(k)) > 0.0){
          inverse_Kj(k, k) = 1.0 / conv_to<double>::from(k_mat.row(j).col(k));
        }
      }
      mat C = inv(inverse_Kj + inverse_sigma);
      genotype_mat.row(j) = (C * (sum(inverse_sigma, 1) * mu(j) + inverse_Kj * tmp.row(j).t()) +
        chol(C) * randn(k_nc, 1)).t();
      
      scale_mat += kron(genotype_mat.row(j).t() - mu(j), 
                        genotype_mat.row(j) - mu(j));
    }
    // Update covariance matrix
    sigma = iwishrnd(scale_mat, dfreedom + genotype_mat.n_rows);
    
    if(save){
      genotype_mat_all.slice(iter) = genotype_mat;
      sigma_all.slice(iter) = sigma;
    }
    if(iter > round(max_iter / 2.0)){
      genotype_mat_avg += genotype_mat;
      sigma_avg += sigma;
    }
    
    // Posterior likelihood
    loglike(iter) = getPostLike(k_mat, transform_mat, genotype_mat, 
                                is_high_cover, mu, sigma, psi, 
                                dfreedom, true);
  } // iter
  genotype_mat_avg /= double (max_iter - round(max_iter / 2.0));
  sigma_avg /= double (max_iter - round(max_iter / 2.0));
  if(!save){
    genotype_mat_all.fill(datum::nan);
    sigma_all.fill(datum::nan);
  }
  return List::create(
    Named("genotype_mat") = genotype_mat_avg,
    Named("genotype_mat_all") = genotype_mat_all,
    Named("sigma") = sigma_avg,
    Named("sigma_all") = sigma_all,
    Named("loglike") = loglike
  );
}
