# QuantileDSPM

## Installation

```R
#install.packages("devtools")
library(devtools)
install_github("wyLI2020/QuantileDSPM")
```

## Usage

```R
fr_estDSPM_ASDBTS_2SLS(tauseq, Y, Ylag1, Z, X, W1, W2, rho_ini, Clevel, Bnum, Me, Mpi)
```

- **tauseq**:  vector, the quantile(s) to be estimated
- **Y**:  $NT$​-dimensional vector, response
- **Ylag1**:  $NT$-dimensional vector, the first lag of response
- **Z**:  $(NT, q)$ matrix, time-varying regressors
- **X**:  $(N, p)$ matrix, time-invariant regressors
- **W1**:  $(N, N)$ matrix, the spatial weights matrix
- **W2**:  $(N, N)$ matrix, the spatial weights matrix
- **rho_ini**:  scalar, the initial value for parameter $\rho$ 
- **Clevel**:  scalar, the confidence level 
- **Bnum**:  integer, the number of bootstrap samples
- **Me**:  1 or 2, the choice for random weights $\{e^{*}_{it}\}$  in the wild bootstrap
- **Mpi**:  1 or 2, the choice for random weights $\{\pi^{*}_{i}\}$ in the random weights bootstrap

```R
noZ_fr_estDSPM_ASDBTS_2SLS(tauseq, Y, Ylag1, X, W1, W2, rho_ini, Clevel, Bnum, Me, Mpi)
```

- **tauseq**:  vector, the quantile(s) to be estimated
- **Y**:  $NT$​-dimensional vector, response
- **Ylag1**:  $NT$-dimensional vector, the first lag of response
- **X**:  $(N, p)$ matrix, time-invariant regressors
- **W1**:  $(N, N)$ matrix, the spatial weights matrix
- **W2**:  $(N, N)$ matrix, the spatial weights matrix
- **rho_ini**:  scalar, the initial value for parameter $\rho$ 
- **Clevel**:  scalar, the confidence level 
- **Bnum**:  integer, the number of bootstrap samples
- **Me**:  1 or 2, the choice for random weights $\{e^{*}_{it}\}$
- **Mpi**:  1 or 2, the choice for random weights $\{\pi^{*}_{i}\}$

## Example

```R
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
```

```R
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

```

