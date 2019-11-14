//------- Source from SClineager_Gibbs_ver2.cpp: do not edit by hand
#ifndef misc
#define misc

arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd = false) 
;
arma::vec diwish_arma(arma::cube x,
                      double df,
                      arma::mat psi,
                      bool logd = false) 
;
#endif
