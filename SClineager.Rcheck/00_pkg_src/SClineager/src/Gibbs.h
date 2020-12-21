//------- Source from SClineager_Gibbs_ver3.cpp: do not edit by hand
#ifndef Gibbs
#define Gibbs

double getPostLike(arma::mat k_mat,
                   arma::mat transform_mat,
                   arma::mat genotype_mat,
                   arma::uvec is_high_cover,
                   arma::vec mu,
                   arma::mat sigma,
                   arma::mat psi,
                   int dfreedom,
                   bool logd)
;
Rcpp::List sclineager_gibbs(
    arma::mat psi,
    arma::mat k_mat,
    arma::mat transform_mat,
    arma::vec mu,
    int dfreedom,
    arma::mat sigma,
    int max_iter = 100,
    bool save = false,
    bool loglike_save = false  
)
;
#endif
