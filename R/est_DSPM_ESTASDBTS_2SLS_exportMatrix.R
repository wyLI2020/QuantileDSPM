# library(Rcpp)
# library(RcppArmadillo)
# library(quantreg)
# library(parallel)
# source('DSPM_EST_fr.R')
# Rcpp::sourceCpp('DSPM_ESTASDBTS_2SLS_fc.cpp')
# source('DSPM_ESTASDBTS_2SLS_fr.R')
# source('DSPM_BTS_2SLS_fr.R')
# Rcpp::sourceCpp('DSPM_BTS_2SLS_fc.cpp')

########################################### Estimation ###########################################
fr_estDSPM_ASDBTS_2SLS <- function(tauseq, Y, Ylag1, Z, X, W1, W2, rho_ini, Clevel, Bnum, Me, Mpi) {
  ## Estimation of \rho, \phi and \varphi(\tau), that is \widehat{\rho}, \check{\rho}, \widetilde{\phi}, \widehat{\phi}, and \widehat{\varphi}(\tau)
  N <- nrow(X)
  T <- nrow(Z) / N
  q <- ncol(Z)
  p <- ncol(X)
  Zbreve <- cbind(Ylag1, as.vector(fc_MatrixMultip(W1, matrix(Ylag1, N, T))), Z) # \breve{Z} = (Y_{-1}, (\iota_{T} \otimes W_{1})Y_{-1}, Z)
  ############################## Estimation Procedure ##############################
  ### Pre-estimation
  phitilde <- fr_phitilde_2SLS(N, T, Y, Ylag1, Z, W1) # A consistent estimator of \phi = (\alpha, \beta, \gamma^{\prime})^{\prime}, \widetilde{\phi}, by the two step least square (2SLS) method
  etatilde <- fc_etaest(phitilde, Y, Zbreve) # A approximate of \eta, \widetilde{\eta}, using \widetilde{\phi}
  Vtilde <- fc_Vtilde(etatilde, N, T) # A approximate of V, \widetilde{V}, using \widetilde{\eta}
  ### Stage 1
  #### \widehat{\rho}_{fo} and \check{\rho}
  secA1 <- fc_secA1(W2, N) # A_{1} \triangleq (I_{T} \otimes secA1)
  secA2 <- W2 # A_{2} = M = I_{T} \otime W_{2}
  Lambda1lambda1 <- fc_Lambda1lambda1(Vtilde, N, T, W2, secA1, secA2)
  Lambda1 <- Lambda1lambda1[1:2,1:2] # \Lambda_{1}
  lambda1 <- Lambda1lambda1[,3] # \lambda_{1}
  rhohatfo <- fr_rhohatfo(N, T, W2, secA1, secA2, Lambda1, lambda1, Vtilde, rho_ini) # The feasible optimal GMM (FOGMM) estimator of \rho, \widehat{\rho}_{fo}, using \widetilde{V} and weighting matrix \widetilde{\Sigma}_{g}^{-1}
  rhocheck <- fr_rhocheck(Lambda1, lambda1) # The linear least square estimator of \rho, \check{\rho}
  ##### ASD of \widehat{\rho}_{fo} and \check{\rho} by theoretical results
  Bhatfo <- fc_B(rhohatfo, N, W2) # B(\widehat{\rho}_{fo})
  Sigma_varepsilonhat <- fc_Sigmavarepsilonest(Bhatfo, Vtilde, N, T) # \widehat{\Sigma}_{\varepsilon} with \widehat{\varepsilon} = S(\widehat{\rho}_{fo}) \widetilde{V}
  asd_rhohatfo <- fr_asd_rhohatfo(rhohatfo, Sigma_varepsilonhat, N, T, secA1, secA2, Lambda1) # estimate of ASD of \widehat{\rho}_{fo}
  Bcheck <- fc_B(rhocheck, N, W2) # B(\check{\rho})
  Sigma_varepsiloncheck <- fc_Sigmavarepsilonest(Bcheck, Vtilde, N, T) # \check{\Sigma}_{\varepsilon} with \check{\varepsilon} = S(\check{\rho}) \widetilde{V}
  asd_rhocheck <- fr_asd_rhocheck(Sigma_varepsiloncheck, N, T, secA1, secA2, Lambda1) # estimate of ASD of \check{\rho}
  #### \widehat{\phi} and \check{\phi}, based on \widehat{\rho}_{fo} and \check{\rho}, respectively
  phihat <- fc_phihat(etatilde, Bhatfo, Sigma_varepsilonhat, N, T, Y, Zbreve) # The GLS estimator of \phi, \widehat{\phi}, using \widetilde{\theta} and \widehat{\rho}_{fo}
  phicheck <- fc_phihat(etatilde, Bcheck, Sigma_varepsiloncheck, N, T, Y, Zbreve) # The GLS estimator of \phi, \check{\phi}, using \widetilde{\theta} and \check{\rho}
  ### Stage 2
  #### \widehat{\varphi}(\tau) based on \widehat{\phi}
  etahat <- fc_etaest(phihat, Y, Zbreve) # A approximate of \eta, \widehat{\eta}, using \widehat{\phi}
  etahat_mat <- matrix(etahat, nrow = N, ncol = T)
  thetahat <- rowMeans(etahat_mat) # \widehat{\theta}
  varphihat <- fr_varphihat(tauseq, p, X, thetahat) # \widehat{\varphi}(\tau) for \tau in tauseq
  #### \check{\varphi}{\tau} based on \check{\phi}
  etacheck <- fc_etaest(phicheck, Y, Zbreve) # A approximate of \eta, \check{\eta}, using \check{\phi}
  etacheck_mat <- matrix(etacheck, nrow = N, ncol = T)
  thetacheck <- rowMeans(etacheck_mat) # \check{\theta}
  varphicheck <- fr_varphihat(tauseq, p, X, thetacheck) # \check{\varphi}(\tau) for \tau in tauseq
  #### ASD of \widehat{\varphi}(\tau) and \check{\varphi}(\tau) by theoretical results
  asd_varphihat <- fc_ASDvarphihat(tauseq, N, p, thetahat, X, varphihat) # the j-th column is the estimated ASD of \widehat{\varphi}(\tau_j)
  asd_varphicheck <- fc_ASDvarphihat(tauseq, N, p, thetacheck, X, varphicheck) # the j-th column is the estimated ASD of \check{\varphi}(\tau_j)
  ### estimates and ASD by theoretical results
  estALL <- c(phitilde, rhohatfo, phihat, as.vector(varphihat), phitilde, rhocheck, phicheck, as.vector(varphicheck))
  asdALL <- c(rep(NA, 2+q), asd_rhohatfo, rep(NA, 2+q), as.vector(asd_varphihat), rep(NA, 2+q), asd_rhocheck, rep(NA, 2+q), as.vector(asd_varphicheck))
  dimALL <- length(estALL)
  ############################## Two-stage Hybrid Bootstrap ##############################
  ### ASD by the two-stage hybrid bootstrapping procedure
  y0 <- Ylag1[1:N]
  rhoest_hc <- c(rhohatfo, rhocheck)
  Best_hc <- array(NA, dim = c(N,N,2)); Best_hc[,,1] <- Bhatfo; Best_hc[,,2] <- Bcheck
  phiest_hc <- matrix(NA, nrow = 2+q, ncol = 2); phiest_hc[,1] <- phihat; phiest_hc[,2] <- phicheck
  etaest_hc <- matrix(NA, nrow = N*T, ncol = 2); etaest_hc[,1] <- etahat; etaest_hc[,2] <- etacheck
  thetaest_hc <- matrix(NA, nrow = N, ncol = 2); thetaest_hc[,1] <- thetahat; thetaest_hc[,2] <- thetacheck
  estALL_BTS_Bnum <- fc_estBnum_2SLS_TS(Bnum, dimALL, Me, Mpi, tauseq, N, T, q, p, y0, Z, X, W1, W2, secA1, secA2, rhoest_hc, Best_hc, phiest_hc, etaest_hc, thetaest_hc, rho_ini)
  asdALL_BTS <- apply(estALL_BTS_Bnum, MARGIN = 1, FUN = sd)
  #### confidence interval
  ##### standard interval
  SIlower_BTS <- estALL - qnorm(1-Clevel/2) * asdALL_BTS # the lower bounds of standard intervals
  SIupper_BTS <- estALL + qnorm(1-Clevel/2) * asdALL_BTS # the upper bounds of standard intervals
  ##### percentile interval
  PIlower_BTS <- apply(estALL_BTS_Bnum, MARGIN = 1, FUN = quantile, probs = Clevel/2, na.rm = TRUE) # the lower bounds of percentile intervals
  PIupper_BTS <- apply(estALL_BTS_Bnum, MARGIN = 1, FUN = quantile, probs = 1-Clevel/2, na.rm = TRUE) # the upper bounds of percentile intervals
  ############################## Results ##############################
  ### results
  results <- c(estALL, asdALL, asdALL_BTS, SIlower_BTS, SIupper_BTS, PIlower_BTS, PIupper_BTS)
  results_mat <- matrix(results, nrow = dimALL, ncol = 7)
  namephitilde <- paste0(c("alpha", "beta", paste0("gamma", 1:q)), "tilde")
  namephihat <- paste0(c("alpha", "beta", paste0("gamma", 1:q)), "hat")
  namephicheck <- paste0(c("alpha", "beta", paste0("gamma", 1:q)), "check")
  namevarphihat <- paste0("varphihat", rep(0:p, length(tauseq)), "_tau", rep(tauseq, each = p+1))
  namevarphicheck <- paste0("varphicheck", rep(0:p, length(tauseq)), "_tau", rep(tauseq, each = p+1))
  rownames(results_mat) <- c(namephitilde, "rhohatfo", namephihat, namevarphihat, namephitilde, "rhocheck", namephicheck, namevarphicheck)
  colnames(results_mat) <- c("est", "ASD", "ASD_BTS", "SIL_BTS", "SIR_BTS", "PIL_BTS", "PIR_BTS")
  return(results_mat)
}


