############################## Two-stage Hybrid Bootstrap ##############################

fr_estbth_2SLS_TS <- function(Me, Mpi, tauseq, N, T, q, p, y0, Z, X, W1, W2, secA1, secA2, rhoest_hc, Best_hc, phiest_hc, etaest_hc, thetaest_hc, rho_ini) {
  ## estimation in the b-th step of the two-stage hybrid bootstrapping procedure
  ### preparation
  #### random weights {e^{\star}_{it}, i=1,...,N, t=1,...,T} in Stage 1, i.i.d. with zero mean, one variance and finite fourth moment
  if (Me == 1) { # If Me=1, we use the standard normal distribution weights
    estar <- rnorm(N*T)
  } else if (Me == 2) { # If Me=2, we use the two-point distribution, takes the value of -1 or 1 with probability 0.5
    estar <- rbinom(N*T, 1, 0.5)*2 - 1
  }
  #### random weights {\pi_{i}, i=1,...,N} in Stage 2, i.i.d. non-negative with mean and variance both equal to 1
  if (Mpi == 1) { # If Mpi=1, we use the standard exponential distribution weights
    pistar <- rexp(N)
  } else if (Mpi == 2) { # If Mpi=2, we use the Rademacher distribution, which takes the value of 0 or 2 with probability 0.5
    pistar <- rbinom(N, 1, 0.5)*2
  }
  #### transform Z to Zarray
  Zarraymid <- array(as.vector(t(Z)), dim = c(q, N, T))
  Zarray <- aperm(Zarraymid,c(2,1,3)) # dim = c(N, q, T)
  ############################## based on rhohat ##############################
  ### Stage 1, , the recursive-design wild bootstrap
  #### bootstrap samples
  rhohatfo <- rhoest_hc[1]; Bhatfo <- Best_hc[,,1]; phihat <- phiest_hc[,1]; etahat <- etaest_hc[,1]; thetahat <- thetaest_hc[,1]
  varepsilonstar_h <- fc_varepsilonstar_TS(estar, N, T, Bhatfo, etahat) # \varepsilon^{\star}, \varepsilon^{\star}_{it} = \widehat{\varepsilon}_{it} e^{\star}_{it}
  Ystar_h_mat <- f_Ystarmat_Generation_TS(N, T, q, y0, rhohatfo, phihat, thetahat, varepsilonstar_h, Zarray, W1, W2) 
  Ystar_h <- as.vector(Ystar_h_mat[,2:(T+1)]) # Y^{\star}
  Ystarlag1_h <- as.vector(Ystar_h_mat[,1:T]) # Y^{\star}_{-1}
  Zbrevestar_h <- cbind(Ystarlag1_h, as.vector(fc_MatrixMultip(W1, matrix(Ystarlag1_h, N, T))), Z) # \breve{Z}^{\star} = (Y^{\star}_{-1}, (\iota_{T} \otimes W_{1})Y^{\star}_{-1}, Z)
  #### bootstrap parameter estimates
  ##### Pre-estimation
  phitildeB_h <- fr_phitilde_2SLS(N, T, Ystar_h, Ystarlag1_h, Z, W1) # \widetilde{\phi}_{B}
  etatildeB_h <- fc_etaest(phitildeB_h, Ystar_h, Zbrevestar_h) # \widetilde{\eta}_{B}
  VtildeB_h <- fc_Vtilde(etatildeB_h, N, T) # \widetilde{V}_{B}
  ##### \widehat{\rho}_{foB}
  LambdaBlambdaB_h <- fc_Lambda1lambda1(VtildeB_h, N, T, W2, secA1, secA2)
  LambdaB_h <- LambdaBlambdaB_h[1:2,1:2] # \Lambda_{B}
  lambdaB_h <- LambdaBlambdaB_h[,3] # \lambda_{B}
  rhohatfoB <- fr_rhohatfo(N, T, W2, secA1, secA2, LambdaB_h, lambdaB_h, VtildeB_h, rho_ini) # \widehat{\rho}_{foB}
  #### \widehat{\phi}_{B}
  BhatfoB <- fc_B(rhohatfoB, N, W2) # B(\widehat{\rho}_{foB})
  Sigma_varepsilonhatB <- fc_Sigmavarepsilonest(BhatfoB, VtildeB_h, N, T) # \widehat{\Sigma}_{\varepsilon B} with \widehat{\varepsilon}_{B} = S(\widehat{\rho}_{foB}) \widetilde{V}_{B}
  phihatB <- fc_phihat(etatildeB_h, BhatfoB, Sigma_varepsilonhatB, N, T, Ystar_h, Zbrevestar_h) # \widehat{\phi}_{B}
  ### Stage 2, the random weight bootstrap
  etahatB <- fc_etaest(phihatB, Ystar_h, Zbrevestar_h) # \widehat{\eta}_{B}
  etahatB_mat <- matrix(etahatB, nrow = N, ncol = T)
  thetahatB <- rowMeans(etahatB_mat) # \widehat{\theta}_{B}
  varphihatB <- fr_varphihatB_TS(pistar, tauseq, p, X, thetahatB) # \widehat{\varphi}_{B}(\tau) for \tau in tauseq
  ############################## based on rhocheck ##############################
  ### Stage 1, , the recursive-design wild bootstrap
  #### bootstrap samples
  rhocheck <- rhoest_hc[2]; Bcheck <- Best_hc[,,2]; phicheck <- phiest_hc[,2]; etacheck <- etaest_hc[,2]; thetacheck <- thetaest_hc[,2]
  varepsilonstar_c <- fc_varepsilonstar_TS(estar, N, T, Bcheck, etacheck) # \varepsilon^{\star}, \varepsilon^{\star}_{it} = \check{\varepsilon}_{it} e^{\star}_{it}
  Ystar_c_mat <- f_Ystarmat_Generation_TS(N, T, q, y0, rhocheck, phicheck, thetacheck, varepsilonstar_c, Zarray, W1, W2)
  Ystar_c <- as.vector(Ystar_c_mat[,2:(T+1)]) # Y^{\star}
  Ystarlag1_c <- as.vector(Ystar_c_mat[,1:T]) # Y^{\star}_{-1}
  Zbrevestar_c <- cbind(Ystarlag1_c, as.vector(fc_MatrixMultip(W1, matrix(Ystarlag1_c, N, T))), Z) # \breve{Z}^{\star} = (Y^{\star}_{-1}, (\iota_{T} \otimes W_{1})Y^{\star}_{-1}, Z)
  #### bootstrap parameter estimates
  ##### Pre-estimation
  phitildeB_c <- fr_phitilde_2SLS(N, T, Ystar_c, Ystarlag1_c, Z, W1) # \widetilde{\phi}_{B}
  etatildeB_c <- fc_etaest(phitildeB_c, Ystar_c, Zbrevestar_c) # \widetilde{\eta}_{B}
  VtildeB_c <- fc_Vtilde(etatildeB_c, N, T) # \widetilde{V}_{B}
  ##### \check{\rho}_{B}
  LambdaBlambdaB_c <- fc_Lambda1lambda1(VtildeB_c, N, T, W2, secA1, secA2)
  LambdaB_c <- LambdaBlambdaB_c[1:2,1:2] # \Lambda_{B}
  lambdaB_c <- LambdaBlambdaB_c[,3] # \lambda_{B}
  rhocheckB <- fr_rhocheck(LambdaB_c, lambdaB_c) # \check{\rho}_{B}
  #### \check{\phi}_{B}
  BcheckB <- fc_B(rhocheckB, N, W2) # B(\check{\rho}_{B})
  Sigma_varepsiloncheckB <- fc_Sigmavarepsilonest(BcheckB, VtildeB_c, N, T) # \check{\Sigma}_{\varepsilon B} with \check{\varepsilon}_{B} = S(\check{\rho}_{B}) \widetilde{V}_{B}
  phicheckB <- fc_phihat(etatildeB_c, BcheckB, Sigma_varepsiloncheckB, N, T, Ystar_c, Zbrevestar_c) # \check{\phi}_{B}
  ### Stage 2, the random weight bootstrap
  etacheckB <- fc_etaest(phicheckB, Ystar_c, Zbrevestar_c) # \check{\eta}_{B}
  etacheckB_mat <- matrix(etacheckB, nrow = N, ncol = T)
  thetacheckB <- rowMeans(etacheckB_mat) # \check{\theta}_{B}
  varphicheckB <- fr_varphihatB_TS(pistar, tauseq, p, X, thetacheckB) # \check{\varphi}_{B}(\tau) for \tau in tauseq
  ### results
  estB <- c(phitildeB_h, rhohatfoB, phihatB, as.vector(varphihatB), phitildeB_c, rhocheckB, phicheckB, as.vector(varphicheckB))
  return(estB)
}


