#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat fc_asmat(arma::vec vec1, int nrow, int ncol) {
  // Fill matrix with elements of vector
  arma::mat vec1_mat(nrow*ncol, 1); vec1_mat.col(0) = vec1;
  vec1_mat.reshape(nrow, ncol);
  return vec1_mat;
}

// [[Rcpp::export]]
arma::vec fc_asvec(arma::mat mat1) {
  // Matrix straighten
  int nrow = mat1.n_rows;
  int ncol = mat1.n_cols;
  mat1.reshape(nrow*ncol, 1);
  return mat1.col(0);
}

// [[Rcpp::export]]
arma::mat fc_iota(int n) {
  // \iota_{n}, n*1 matrix of ones
  arma::mat iota_n(n, 1); iota_n.fill(1.0);
  return iota_n;
}

// [[Rcpp::export]]
arma::mat fc_I(int n) {
  // I_{n}, n*n identity matrix
  arma::mat I_n(n,n); I_n.eye(n,n);
  return I_n;
}

// [[Rcpp::export]]
arma::mat fc_B(double rho, int N, arma::mat W2) {
  // B(\rho) = I_{N} - \rho W_{2}
  arma::mat I_N = fc_I(N);
  arma::mat B = I_N - rho * W2;
  return B;
}

// // [[Rcpp::export]]
// arma::mat fc_S(double rho, int N, int T, arma::mat W2) {
//   // S(\rho) = I_{T} \otimes B(\rho)
//   arma::mat I_T = fc_I(T);
//   arma::mat B = fc_B(rho, N, W2);
//   arma::mat S = kron(I_T, B);
//   return S;
// }

// [[Rcpp::export]]
arma::mat fc_MatrixMultip(arma::mat M1, arma::mat M2) {
  // matrix multiplication
  arma::mat multipM = M1 * M2;
  return multipM;
}

////////////////////////////// Estimations //////////////////////////////

// [[Rcpp::export]]
arma::vec fc_phitilde_2SLS(arma::vec wrD, arma::mat wrH, arma::mat wrH_pred) {
  // In pre-estimation, calculate a consistent estimator of \phi = (\alpha, \beta, \gamma^{\prime})^{\prime}, \widetilde{\phi}, by the two step least square (2SLS) method
  arma::vec phitilde = (wrH_pred.t() * wrH).i() * wrH_pred.t() * wrD;
  return phitilde;
}

// [[Rcpp::export]]
arma::vec fc_etaest(arma::vec phiest, arma::vec Y, arma::mat Zbreve) {
  // In pre-estimation and stage 2, calculate the approximate of \eta, \widetilde{\eta} or \widehat{\eta}, using \widetilde{\phi} or \widehat{\phi}, respectively
  arma::vec etaest = Y - Zbreve * phiest;
  return etaest;
}

// [[Rcpp::export]]
arma::vec fc_Vtilde(arma::vec etatilde, int N, int T) {
  // In pre-estimation, calculate the approximate of V, \widetilde{V}, using \widetilde{\eta}
  arma::mat etatilde_mat = fc_asmat(etatilde, N, T);
  arma::mat iota_T = fc_iota(T);
  arma::vec Vtilde = etatilde - kron(iota_T, mean(etatilde_mat, 1));
  return Vtilde;
}

// [[Rcpp::export]]
arma::mat fc_secA1(arma::mat W2, int N) {
  // In Stage 1, A_{1} = I_{T} \otimes (W_{2}^{\prime} W_{2} - diagW2) \triangleq (I_{T} \otimes secA1)
  arma::mat diagW2(N,N); diagW2.fill(0.0); // diagW2 \triangleq diag\{w_{2,.1}^{\prime}w_{2,.1}, ..., w_{2,.N}^{\prime}w_{2,.N}\}
  for(int i = 0; i < N; i++) {
    arma::vec ithcol = W2.col(i);
    diagW2(i,i) = dot(ithcol, ithcol);
  }
  arma::mat secA1 = W2.t() * W2 - diagW2; // secA1 \triangleq (W_{2}^{\prime} W_{2} - diagW2)
  return secA1;
}

// [[Rcpp::export]]
arma::mat fc_Lambda1lambda1(arma::vec Vtilde, int N, int T, arma::mat W2, arma::mat secA1, arma::mat secA2) {
  // In Stage 1, c(\Lambda_{1}, \lambda_{1})
  arma::mat Vtilde_mat = fc_asmat(Vtilde, N, T);
  arma::mat M_Vtilde_mat = W2 * Vtilde_mat; // M \widetilde{V} = (I_{T} \otimes W_{2}) \widetilde{V}
  arma::mat A1_Vtilde_mat = secA1 * Vtilde_mat; // A_{1} \widetilde{V} = (I_{T} \otimes secA1) \widetilde{V}
  arma::mat A2_Vtilde_mat = secA2 * Vtilde_mat; // A_{1} \widetilde{V} = (I_{T} \otimes secA2) \widetilde{V}
  arma::mat A1_M_Vtilde_mat = secA1 * M_Vtilde_mat; // A_{1} M \widetilde{V}
  arma::mat A2_M_Vtilde_mat = secA2 * M_Vtilde_mat; // A_{2} M \widetilde{V}
  arma::mat A2T_M_Vtilde_mat = secA2.t() * M_Vtilde_mat; // A_{2}^{\prime} M \widetilde{V}
  arma::vec M_Vtilde = fc_asvec(M_Vtilde_mat); 
  arma::vec A1_Vtilde = fc_asvec(A1_Vtilde_mat); 
  arma::vec A2_Vtilde = fc_asvec(A2_Vtilde_mat); 
  arma::vec A1_M_Vtilde = fc_asvec(A1_M_Vtilde_mat); 
  arma::vec A2_M_Vtilde = fc_asvec(A2_M_Vtilde_mat); 
  arma::vec A2T_M_Vtilde = fc_asvec(A2T_M_Vtilde_mat); 
  // arma::mat Lambda1(2,2); Lambda1.fill(0.0); // matrix \Lambda_{1}
  // arma::vec lambda1(2); // vector \lambda_{1}
  // Lambda1(0,0) = dot(Vtilde, A1_M_Vtilde) * 2.0 / (N*T*1.0);
  // Lambda1(1,0) = dot(Vtilde, A2_M_Vtilde+A2T_M_Vtilde) / (N*T*1.0);
  // Lambda1(0,1) = - dot(M_Vtilde, A1_M_Vtilde) / (N*T*1.0);
  // Lambda1(1,1) = - dot(M_Vtilde, A2_M_Vtilde) / (N*T*1.0);
  // lambda1(0) = dot(Vtilde, A1_Vtilde) / (N*T*1.0);
  // lambda1(1) = dot(Vtilde, A2_Vtilde) / (N*T*1.0);
  arma::mat Lambda1lambda1(2,3); Lambda1lambda1.fill(0.0);
  Lambda1lambda1(0,0) = dot(Vtilde, A1_M_Vtilde) * 2.0 / (N*T*1.0);
  Lambda1lambda1(1,0) = dot(Vtilde, A2_M_Vtilde+A2T_M_Vtilde) / (N*T*1.0);
  Lambda1lambda1(0,1) = - dot(M_Vtilde, A1_M_Vtilde) / (N*T*1.0);
  Lambda1lambda1(1,1) = - dot(M_Vtilde, A2_M_Vtilde) / (N*T*1.0);
  Lambda1lambda1(0,2) = dot(Vtilde, A1_Vtilde) / (N*T*1.0);
  Lambda1lambda1(1,2) = dot(Vtilde, A2_Vtilde) / (N*T*1.0);
  return Lambda1lambda1;
}

// [[Rcpp::export]]
arma::mat fc_Sigmavarepsilonest(arma::mat Best, arma::vec Vtilde, int N, int T) {
  // In Stage 1, calculate the estimate of \Sigma_{\varepsilon}
  arma::mat Vtilde_mat = fc_asmat(Vtilde, N, T);
  arma::mat varepsilonest_mat = Best * Vtilde_mat; 
  arma::vec sigma2varepsilonest = var(varepsilonest_mat, 0, 1);
  arma::mat Sigmavarepsilonest(N,N); Sigmavarepsilonest.fill(0.0); Sigmavarepsilonest.diag() = sigma2varepsilonest;
  return Sigmavarepsilonest;
}

// [[Rcpp::export]]
arma::mat fc_Sigmagest(arma::mat Sigma_varepsilonest, arma::mat secA1, arma::mat secA2, int N) {
  // In Stage 1, calculate the estimate of \Sigma_{g}
  arma::mat Sigmagest(2,2); Sigmagest.fill(0.0);
  double sum00 = 0.0; double sum01 = 0.0; double sum11 = 0.0; 
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      sum00 = sum00 + (secA1(i,j) + secA1(j,i)) * (secA1(i,j) + secA1(j,i)) * Sigma_varepsilonest(i,i) * Sigma_varepsilonest(j,j);
      sum01 = sum01 + (secA1(i,j) + secA1(j,i)) * (secA2(i,j) + secA2(j,i)) * Sigma_varepsilonest(i,i) * Sigma_varepsilonest(j,j);
      sum11 = sum11 + (secA2(i,j) + secA2(j,i)) * (secA2(i,j) + secA2(j,i)) * Sigma_varepsilonest(i,i) * Sigma_varepsilonest(j,j);
    }
  }
  Sigmagest(0,0) = sum00 / (2.0*N);
  Sigmagest(0,1) = sum01 / (2.0*N);
  Sigmagest(1,0) = sum01 / (2.0*N);
  Sigmagest(1,1) = sum11 / (2.0*N);
  return Sigmagest;
}

// [[Rcpp::export]]
arma::vec fc_kappacheck(arma::mat Lambda1, arma::vec lambda1) {
  // In Stage 1, calculate the linear least square estimator of \kappa = (\rho, \rho^{2})^{\prime}, \widecheck{\kappa}, using \widetilde{V}
  arma::vec kappacheck = Lambda1.i() * lambda1;
  return kappacheck;
}

// [[Rcpp::export]]
arma::vec fc_phihat(arma::vec etatilde, arma::mat Bhat, arma::mat Sigma_varepsilonhat, int N, int T, arma::vec Y, arma::mat Zbreve) {
  // In Stage 1, calculate the GLS estimator of \phi, \widehat{\phi}, using \widetilde{\theta}
  int dimphi = Zbreve.n_cols;
  arma::mat etatilde_mat = fc_asmat(etatilde, N, T);
  arma::vec thetatilde = mean(etatilde_mat, 1); // \widetilde{\theta}
  arma::mat iota_T = fc_iota(T);
  arma::vec Ybartilde = Y - kron(iota_T, thetatilde); // \widetilde{\overline{Y}} \triangleq (Y - \iota_{T} \otimes \widetilde{\theta})
  arma::mat Ybartilde_mat = fc_asmat(Ybartilde, N, T);
  arma::mat secOmegahatInverse = Bhat.t() * Sigma_varepsilonhat.i() * Bhat; // \widehat{\Omega}^{-1} = I_{T} \otimes (B^{\prime}(\widehat{\rho} \widehat{\Sigma}_{\varepsilon}^{-1} B(\widehat{\rho}))
  arma::mat OmegahatInverse_Ybartilde_mat = secOmegahatInverse * Ybartilde_mat; // \widehat{\Omega}^{-1} (Y - \iota_{T} \otimes \widetilde{\theta})
  arma::vec OmegahatInverse_Ybartilde = fc_asvec(OmegahatInverse_Ybartilde_mat);
  arma::mat OmegahatInverse_Zbreve(N*T, dimphi); OmegahatInverse_Zbreve.fill(0.0); // \widehat{\Omega}^{-1} \breve{Z}
  for(int l = 0; l < dimphi; l++) {
    arma::mat Zbrevel_mat = fc_asmat(Zbreve.col(l), N, T); // the l-th column of \breve{Z}
    arma::mat OmegahatInverse_Zbrevel_mat = secOmegahatInverse * Zbrevel_mat;
    arma::vec OmegahatInverse_Zbrevel = fc_asvec(OmegahatInverse_Zbrevel_mat);
    OmegahatInverse_Zbreve.col(l) = OmegahatInverse_Zbrevel;
  }
  arma::vec phihat = (Zbreve.t() * OmegahatInverse_Zbreve).i() * Zbreve.t() * OmegahatInverse_Ybartilde; // \widehat{\phi}
  return phihat;
}

////////////////////////////// ASD by Theoretical Results //////////////////////////////

// [[Rcpp::export]]
double IQR_rcpp(arma::vec x){
  // IQR_rcpp - IQR() in R, interquartile range of the x values
  Function f("IQR");
  NumericVector result = f(Named("x") = x);
  return result(0);
}

// [[Rcpp::export]]
double dnorm_rcpp(double x){
  // dnorm_rcpp - dnorm() in R, Gaussian kernel function
  Function f("dnorm");
  NumericVector result = f(Named("x") = x);
  return result(0);
}

// [[Rcpp::export]]
arma::vec fc_Xihattauseq(arma::vec tauseq, int N, arma::vec thetahat, arma::mat Xwith1, arma::mat varphihat) { // Gaussian kernel
  // In ASD of \widehat{\varphi}(\tau), estimates of \Xi(\tau) = f_{\xi(\tau)}(0|X) where \xi(\tau) = \theta - X \varphi_{0}(\tau), for \tau in tauseq
  int taunum = tauseq.n_elem;
  arma::vec Xihattauseq(taunum); // the j-th element is \widehat{Xi}(\tau_{j})
  for (int j = 0; j < taunum; j++) {
    // double tau = tauseq(j); // \tau
    arma::vec varphitauhat = varphihat.col(j); // \widehat{\varphi}(\tau)
    arma::vec xitauhat = thetahat - Xwith1 * varphitauhat; // \widehat{Xi}(\tau) = \widehat{\theta} - X \widehat{\varphi}(\tau)
    // bandwidth h
    double s = stddev(xitauhat); // standard deviation of the \widehat{Xi}(\tau)
    double Rhat = IQR_rcpp(xitauhat); // interquartile of the \widehat{Xi}(\tau)
    arma::vec vecFORmin(2); vecFORmin(0) = s; vecFORmin(1) = Rhat/1.34;
    double h = 0.9 * pow((N*1.0), -0.2) * min(vecFORmin); 
    // K((0 - \widehat{xi}_{i}(\tau)) / h), i=1,...,N
    arma::vec K(N); // the i-th element is K((0 - \widehat{\xi}_{i}(\tau)) / h)
    for (int i = 0; i < N; i++) {
      K(i) = dnorm_rcpp((0.0 - xitauhat(i)) / h);
    }
    // \widehat{\Xi}(\tau)
    double Xitauhat = sum(K) / (N*h*1.0); // \widehat{\Xi}(\tau) = \frac{1}{Nh} \sum_{i=1}^{N} K((0 - \widehat{\xi}_{i}(\tau)) / h)
    Xihattauseq(j) = Xitauhat;
  }
  return Xihattauseq;
}

// [[Rcpp::export]]
arma::mat fc_ASDvarphihat(arma::vec tauseq, int N, int p, arma::vec thetahat, arma::mat X, arma::mat varphihat){
  // estimate of ASD of \widehat{\varphi}(\tau), for \tau in tauseq
  int taunum = tauseq.n_elem;
  int dim = 1+p; 
  arma::mat iota_N = fc_iota(N);
  arma::mat Xwith1 = join_rows(iota_N, X);
  arma::vec Xihattauseq = fc_Xihattauseq(tauseq, N, thetahat, Xwith1, varphihat); // the j-th element is \widehat{\Xi}(\tau_j)
  arma::mat DN(dim,dim); DN.fill(0.0); // D_{N} = \frac{1}{N} \sum_{i=1}^{N} x_{i} x_{i}^{\prime}
  for (int i = 0; i < N; i++) {
    arma::mat xiprime = Xwith1.row(i);
    arma::mat xi = xiprime.t();
    DN = DN + (xi * xiprime / (N*1.0));
  }
  arma::mat asd_varphihat(dim, taunum); asd_varphihat.fill(0.0); // the j-th column is the estimate of ASD of \widehat{\varphi}(\tau_j)
  for (int j = 0; j < taunum; j++) {
    double tau = tauseq(j);
    double Xitauhat = Xihattauseq(j);
    arma::mat CVM_varphitauhat = DN.i() * tau * (1.0-tau) / (Xitauhat * Xitauhat * N * 1.0); // estimate of the covariance matrix of \widehat{\varphi}(\tau_j), that is \frac{1}{N} \tau_j(1-\tau_j) \widehat{Xi}^{-2}(\tau_j) D_{N}^{-1}
    arma::vec asd2_varphitauhat = CVM_varphitauhat.diag(); // square of estimated ASD
    arma::vec asd_varphitauhat = sqrt(asd2_varphitauhat); // estimated ASD of \widehat{\varphi}(\tau)
    asd_varphihat.col(j) = asd_varphitauhat;
  }
  return asd_varphihat;
}

////////////////////////////// ASD by Two-stage Hybrid Bootstrap //////////////////////////////

// [[Rcpp::export]]
arma::mat f_Ystarmat_Generation_TS(int N, int T, int q, arma::vec y0, double rhohat, arma::vec phihat, arma::vec thetahat, arma::vec varepsilonstar, arma::cube Zarray, arma::mat W1, arma::mat W2) {
  // the matrix of (y^{\star}_{0}, y^{\star}_{1}, ..., y^{\star}_{T})
  double alphahat = phihat(0);
  double betahat = phihat(1);
  arma::vec gammahat = phihat.subvec(2, q+1);
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::mat Bhat = I_N - rhohat * W2;
  arma::mat varepsilonstar_mat(N*T, 1); varepsilonstar_mat.col(0) = varepsilonstar;
  varepsilonstar_mat.reshape(N, T);
  // Y^{\star}
  arma::mat Ystar_mat(N, T+1); Ystar_mat.fill(0.0);
  Ystar_mat.col(0) = y0;
  for(int t = 1; t < T+1; t++) {
    Ystar_mat.col(t) = thetahat + alphahat * Ystar_mat.col(t-1) + betahat * W1 * Ystar_mat.col(t-1) + Zarray.slice(t-1) * gammahat + Bhat.i() * varepsilonstar_mat.col(t-1);
  }
  return Ystar_mat;
}

// [[Rcpp::export]]
arma::vec fc_varepsilonstar_TS(arma::vec estar, int N, int T, arma::mat Bhat, arma::vec etahat) {
  // \varepsilon^{\star} with \varepsilon^{\star}_{it} = \widehat{\varepsilon}_{it} e^{\star}_{it}
  arma::vec Vhat = fc_Vtilde(etahat, N, T); // \widehat{V}, using \widehat{\eta}
  arma::mat Vhat_mat = fc_asmat(Vhat, N, T);
  arma::mat varepsilonhat_mat = Bhat * Vhat_mat;
  arma::vec varepsilonhat = fc_asvec(varepsilonhat_mat); // \widehat{\varepsilon}
  arma::vec varepsilonstar = varepsilonhat % estar; // \varepsilon^{\star}, \varepsilon^{\star}_{it} = \widehat{\varepsilon}_{it} e^{\star}_{it}
  return varepsilonstar;
}


