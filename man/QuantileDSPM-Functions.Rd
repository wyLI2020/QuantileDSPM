\name{QuantileDSPM-Functions}
\alias{IQR_rcpp}
\alias{LossFUN_rho}
\alias{dnorm_rcpp}
\alias{f_Ystarmat_Generation_TS}
\alias{fc_ASDvarphihat}
\alias{fc_B}
\alias{fc_I}
\alias{fc_Lambda1lambda1}
\alias{fc_MatrixMultip}
\alias{fc_Sigmagest}
\alias{fc_Sigmavarepsilonest}
\alias{fc_Vtilde}
\alias{fc_Xihattauseq}
\alias{fc_asmat}
\alias{fc_asvec}
\alias{fc_estBnum_2SLS_TS}
\alias{fc_estbth_2SLS_TS}
\alias{fc_etaest}
\alias{fc_iota}
\alias{fc_kappacheck}
\alias{fc_phihat}
\alias{fc_phitilde_2SLS}
\alias{fc_secA1}
\alias{fc_varepsilonstar_TS}
\alias{fr_asd_rhocheck}
\alias{fr_asd_rhohatfo}
\alias{fr_estDSPM_ASDBTS_2SLS}
\alias{fr_estbth_2SLS_TS}
\alias{fr_phitilde_2SLS}
\alias{fr_rhocheck}
\alias{fr_rhohatfo}
\alias{fr_varphihat}
\alias{fr_varphihatB_TS}
\alias{noZ_IQR_rcpp}
\alias{noZ_dnorm_rcpp}
\alias{noZ_f_Ystarmat_Generation_TS}
\alias{noZ_fc_ASDvarphihat}
\alias{noZ_fc_B}
\alias{noZ_fc_I}
\alias{noZ_fc_Lambda1lambda1}
\alias{noZ_fc_MatrixMultip}
\alias{noZ_fc_Sigmagest}
\alias{noZ_fc_Sigmavarepsilonest}
\alias{noZ_fc_Vtilde}
\alias{noZ_fc_Xihattauseq}
\alias{noZ_fc_asmat}
\alias{noZ_fc_asvec}
\alias{noZ_fc_estBnum_2SLS_TS}
\alias{noZ_fc_estbth_2SLS_TS}
\alias{noZ_fc_etaest}
\alias{noZ_fc_iota}
\alias{noZ_fc_kappacheck}
\alias{noZ_fc_phihat}
\alias{noZ_fc_phitilde_2SLS}
\alias{noZ_fc_secA1}
\alias{noZ_fc_varepsilonstar_TS}
\alias{noZ_fr_asd_rhocheck}
\alias{noZ_fr_asd_rhohatfo}
\alias{noZ_fr_estDSPM_ASDBTS_2SLS}
\alias{noZ_fr_estbth_2SLS_TS}
\alias{noZ_fr_phitilde_2SLS}
\alias{noZ_fr_rhocheck}
\alias{noZ_fr_rhohatfo}
\alias{noZ_fr_varphihat}
\alias{noZ_fr_varphihatB_TS}
\title{Estimation and bootstrapping procedures}
\description{
  Estimation and bootstrapping procedures
}
\usage{
  fr_estDSPM_ASDBTS_2SLS(tauseq, Y, Ylag1, Z, X, W1, W2, rho_ini, Clevel, Bnum, Me, Mpi);
  noZ_fr_estDSPM_ASDBTS_2SLS(tauseq, Y, Ylag1, X, W1, W2, rho_ini, Clevel, Bnum, Me, Mpi)
}
\arguments{
  \item{tauseq}{vector, the quantile(s) to be estimated}
  \item{Y}{\emph{NT}-dimensional vector, response}
  \item{Ylag1}{\emph{NT}-dimensional vector, the first lag of response}
  \item{Z}{(\emph{NT}, q) matrix, time-varying regressors}
  \item{X}{(\emph{N}, p) matrix, time-invariant regressors}
  \item{W1}{(\emph{N}, \emph{N}) matrix, the spatial weights matrix}
  \item{W2}{(\emph{N}, \emph{N}) matrix, the spatial weights matrix}
  \item{rho_ini}{scalar, the initial value for parameter}
  \item{Clevel}{scalar, the confidence level}
  \item{Bnum}{integer, the number of bootstrap samples}
  \item{Me}{1 or 2, the choice for random weights in the wild bootstrap}
  \item{Mpi}{1 or 2, the choice for random weights in the random weights bootstrap}
}
\value{
  \code{fr_estDSPM_ASDBTS_2SLS} and \code{noZ_fr_estDSPM_ASDBTS_2SLS} return matrices.
}
\details{
  With the preproccess in following "Examples", we can get the estimators of the DSP model.
}
\examples{
  library(QuantileDSPM)
  data("Dataset_CityAirQuality")
  W <- Dataset_CityAirQuality$W
  Y <- Dataset_CityAirQuality$Y
  Ylag1 <- Dataset_CityAirQuality$Ylag1
  X <- Dataset_CityAirQuality$X
  Z <- Dataset_CityAirQuality$Z
  set.seed(6)
  result_mat <- fr_estDSPM_ASDBTS_2SLS(tauseq=c(0.25,0.75), Y, Ylag1, Z, X, W, W, rho_ini=0.5, Clevel=0.1, Bnum=1000, Me=1, Mpi=1)
  result <- matrix(NA, nrow = nrow(result_mat), ncol = 4)
  result[,c(1,2)] <- result_mat[,c(1,3)]
  result[,3] <- result[,1] / result[,2]
  result[,4] <- 2*(1-pnorm(abs(result[,1] / result[,2])))
  rownames(result) <- rownames(result_mat); colnames(result) <- c("est", "ASD", "z", "p-value")
  round(result, 3)
  
  library(QuantileDSPM)
  data("Dataset_HousingPrices")
  W <- Dataset_HousingPrices$W
  Y <- Dataset_HousingPrices$Y
  Ylag1 <- Dataset_HousingPrices$Ylag1
  X <- Dataset_HousingPrices$X
  Z <- Dataset_HousingPrices$Z
  set.seed(6)
  result_mat <- noZ_fr_estDSPM_ASDBTS_2SLS(tauseq=c(0.1, 0.9), Y, Ylag1, X, W, W, rho_ini=0.5, Clevel=0.1, Bnum=1000, Me=1, Mpi=1)
  result <- matrix(NA, nrow = nrow(result_mat), ncol = 4)
  result[,c(1,2)] <- result_mat[,c(1,3)]
  result[,3] <- result[,1] / result[,2]
  result[,4] <- 2*(1-pnorm(abs(result[,1] / result[,2])))
  rownames(result) <- rownames(result_mat); colnames(result) <- c("est", "ASD", "z", "p-value")
  round(result, 3)
}
% \references{
  % See the documentation for Armadillo, and RcppArmadillo, for more details.
% }
\author{Wenyu Li}

