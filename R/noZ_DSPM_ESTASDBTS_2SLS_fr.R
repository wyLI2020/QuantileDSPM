############################## Estimations ##############################

noZ_fr_phitilde_2SLS <- function(N, T, Y, Ylag1, W1) {
  # \widetilde{\phi} by 2SLS
  wrY <- Y[(2*N+1):(N*T)] # Y^{\wr}
  wrYlag1 <- Y[(N+1):(N*(T-1))] # Y^{\wr}_{-1}
  wrYlag2 <- Y[1:(N*(T-2))] # Y^{\wr}_{-2}
  wrYlag3 <- Ylag1[1:(N*(T-2))] # Y^{\wr}_{-3}
  wrD <- wrY - wrYlag1 # D^{\wr} = Y^{\wr} - Y^{\wr}_{-1}
  wrH1 <- wrYlag1 - wrYlag2 # the 1st section of H^{\wr}, D^{\wr}_{-1} = Y^{\wr}_{-1} - Y^{\wr}_{-2}
  wrH2 <- as.vector(noZ_fc_MatrixMultip(W1, matrix(wrH1, N, T-2))) # the 2nd section of H^{\wr}, (\iota_{T-2} \otimes W_{1}) D^{\wr}_{-1}
  wrH <- cbind(wrH1, wrH2) # H^{\wr} = (D^{\wr}_{-1}, (\iota_{T-2} \otimes W_{1}) D^{\wr}_{-1}, Z^{\wr} - Z^{\wr}_{-1})
  IV_1 <- wrYlag2 - wrYlag3 # the 1st section of instrumental variable matrix, D^{\wr}_{-2} = Y^{\wr}_{-2} - Y^{\wr}_{-3}
  IV_2 <- as.vector(noZ_fc_MatrixMultip(W1, matrix(IV_1, N, T-2))) # the 2nd section of instrumental variable matrix, (\iota_{T-2} \otimes W_{1}) D^{\wr}_{-2}
  IV_df <- as.data.frame(cbind(IV_1, IV_2)) # instrumental variable matrix to predict H^{\wr} (specifically, D^{\wr}_{-1}), R^{\wr} = H^{\wr}_{-1} = (D^{\wr}_{-2}, (\iota_{T-2} \otimes W_{1}) D^{\wr}_{-2}, Z^{\wr}_{-1} - Z^{\wr}_{-2})
  lmwrH1 <- lm(wrH1 ~ ., data = IV_df)
  wrH1_pred <- predict(lmwrH1, IV_df) # \widetilde{D}^{\wr}_{-1}, predictor of D^{\wr}_{-1} by instrumental variables
  wrH2_pred <- as.vector(noZ_fc_MatrixMultip(W1, matrix(wrH1_pred, N, T-2))) # (\iota_{T-2} \otimes W_{1}) \widetilde{D}^{\wr}_{-1}
  wrH_pred <- cbind(wrH1_pred, wrH2_pred) # predictor of H^{\wr}, \widetilde{H}^{\wr} = (\widetilde{D}^{\wr}_{-1}, (\iota_{T-2} \otimes W_{1}) \widetilde{D}^{\wr}_{-1}, Z^{\wr} - Z^{\wr}_{-1})
  phitilde <- noZ_fc_phitilde_2SLS(wrD, wrH, wrH_pred) # A consistent estimator of \phi = (\alpha, \beta, \gamma^{\prime})^{\prime}, \widetilde{\phi}, by the two step least square (2SLS) method
  return(phitilde)
}

noZ_fr_rhohatfo <- function(N, T, W2, secA1, secA2, Lambda1, lambda1, Vtilde, rho_ini) {
  # \widehat{\rho}_{fo}
  rhohatI <- optim(par = rho_ini, fn = LossFUN_rho, method = "L-BFGS-B", Lambda1=Lambda1, lambda1=lambda1, aTa=diag(2))$par # The GMM estimator of \rho, \widehat{\rho}_{I} = \widehat{\rho}(I_{2}), using \widetilde{V} and weighting matrix I_{2}
  BhatI <- noZ_fc_B(rhohatI, N, W2) # B(\widehat{\rho}_{I})
  Sigma_varepsilontilde <- noZ_fc_Sigmavarepsilonest(BhatI, Vtilde, N, T) # \widetilde{\Sigma}_{\varepsilon} with \widetilde{\varepsilon} = S(\widehat{\rho}_{I}) \widetilde{V}
  Sigma_gtilde <- noZ_fc_Sigmagest(Sigma_varepsilontilde, secA1, secA2, N) # \widetilde{\Sigma}_{g}
  rhohatfo <- optim(par = rho_ini, fn = LossFUN_rho, method = "L-BFGS-B", Lambda1=Lambda1, lambda1=lambda1, aTa=solve(Sigma_gtilde))$par # The feasible optimal GMM (FOGMM) estimator of \rho, \widehat{\rho}_{fo}, using \widetilde{V} and weighting matrix \widetilde{\Sigma}_{g}^{-1}
  return(rhohatfo)
}

noZ_fr_rhocheck <- function(Lambda1, lambda1) {
  # \check{\rho}
  kappacheck <- noZ_fc_kappacheck(Lambda1, lambda1) # The linear least square estimator of \kappa = (\rho, \rho^{2})^{\prime}, \check{\kappa}, using \widetilde{V}
  rhocheck <- kappacheck[1] # \check{\rho}
  return(rhocheck)
}

noZ_fr_varphihat <- function(tauseq, p, X, thetahat) {
  # \widehat{\varphi}(\tau) for \tau in tauseq
  varphihat <- matrix(nrow = 1+p, ncol = length(tauseq)) # \widehat{\varphi}(\tau) for \tau in tauseq
  if (length(tauseq) == 1) {
    varphihat[,1] <- rq(thetahat ~ X, tau = tauseq)$coefficients
  } else {
    varphihat <- rq(thetahat ~ X, tau = tauseq)$coefficients
  }
  return(varphihat)
}

############################## ASD by Theoretical Results ##############################

noZ_fr_asd_rhohatfo <- function(rhohatfo, Sigma_varepsilonhat, N, T, secA1, secA2, Lambda1) {
  # ASD of \widehat{\rho}_{fo} by theoretical results
  Sigma_ghat <- noZ_fc_Sigmagest(Sigma_varepsilonhat, secA1, secA2, N) # \widehat{\Sigma}_{g}
  dghat <- Lambda1 %*% c(1, 2*rhohatfo) # \widehat{d}_{g}, where d_{g} = \Lambda (1, 2 \rho_{0})^{\prime}
  asd_rhohatfo <- sqrt(1/crossprod(dghat, solve(Sigma_ghat)%*%dghat) / (N*T)) # estimate of ASD of \widehat{\rho}_{fo}
  return(asd_rhohatfo)
}

noZ_fr_asd_rhocheck <- function(Sigma_varepsiloncheck, N, T, secA1, secA2, Lambda1) {
  # ASD of \check{\rho} by theoretical results
  Sigma_gcheck <- noZ_fc_Sigmagest(Sigma_varepsiloncheck, secA1, secA2, N) # \check{\Sigma}_{g}
  cov_kappacheck <- solve(crossprod(Lambda1, solve(Sigma_gcheck)%*%Lambda1)) / (N*T) # estimate of covariance matrix of \check{\kappa}
  asd_rhocheck <- sqrt(cov_kappacheck[1,1]) # estimate of ASD of \check{\rho}
  return(asd_rhocheck)
}

############################## ASD by Two-stage Hybrid Bootstrap ##############################

noZ_fr_varphihatB_TS <- function(pistar, tauseq, p, X, thetahatB) {
  # \widehat{\varphi}_{B}(\tau) for \tau in tauseq, by random weights bootstrap
  varphihatB <- matrix(nrow = 1+p, ncol = length(tauseq)) # \widehat{\varphi}(\tau) for \tau in tauseq
  if (length(tauseq) == 1) {
    varphihatB[,1] <- rq(thetahatB ~ X, tau = tauseq, weights = pistar)$coefficients
  } else {
    varphihatB <- rq(thetahatB ~ X, tau = tauseq, weights = pistar)$coefficients
  }
  return(varphihatB)
}


