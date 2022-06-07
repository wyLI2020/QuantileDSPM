#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//////////////////////////////////////// Two-stage Hybrid Bootstrap ////////////////////////////////////////

// [[Rcpp::export]]
arma::vec fc_estbth_2SLS_TS(int Me, int Mpi, arma::vec tauseq, int N, int T, int q, int p, arma::vec y0, arma::mat Z, arma::mat X, arma::mat W1, arma::mat W2, arma::mat secA1, arma::mat secA2, arma::vec rhoest_hc, arma::cube Best_hc, arma::mat phiest_hc, arma::mat etaest_hc, arma::mat thetaest_hc, double rho_ini) {
  // fr_estbth_2SLS_TS in DSPM_BTS_2SLS_fr.R, estimation in the b-th step of the two-stage hybrid bootstrapping procedure
  Environment OPT = Environment::namespace_env("QuantileDSPM");
  Function f = OPT["fr_estbth_2SLS_TS"];
  NumericVector estB_mid = f(Named("Me")=Me, Named("Mpi")=Mpi, Named("tauseq")=tauseq, Named("N")=N, Named("T")=T, Named("q")=q, Named("p")=p, Named("y0")=y0, Named("Z")=Z, Named("X")=X, Named("W1")=W1, Named("W2")=W2, Named("secA1")=secA1, Named("secA2")=secA2, Named("rhoest_hc")=rhoest_hc, Named("Best_hc")=Best_hc, Named("phiest_hc")=phiest_hc, Named("etaest_hc")=etaest_hc, Named("thetaest_hc")=thetaest_hc, Named("rho_ini")=rho_ini);
  arma::vec estB(estB_mid.begin(), estB_mid.size(), false);
  return estB;
}

// [[Rcpp::export]]
arma::mat fc_estBnum_2SLS_TS(int Bnum, int dimALL, int Me, int Mpi, arma::vec tauseq, int N, int T, int q, int p, arma::vec y0, arma::mat Z, arma::mat X, arma::mat W1, arma::mat W2, arma::mat secA1, arma::mat secA2, arma::vec rhoest_hc, arma::cube Best_hc, arma::mat phiest_hc, arma::mat etaest_hc, arma::mat thetaest_hc, double rho_ini) {
  // estimates in Bnum bootstrap, the two-stage hybrid bootstrapping procedure
  arma::mat estBnum(dimALL, Bnum); estBnum.fill(0.0);
  for(int b = 0; b < Bnum; b++) {
    estBnum.col(b) = fc_estbth_2SLS_TS(Me, Mpi, tauseq, N, T, q, p, y0, Z, X, W1, W2, secA1, secA2, rhoest_hc, Best_hc, phiest_hc, etaest_hc, thetaest_hc, rho_ini);
  }
  return estBnum;
}


