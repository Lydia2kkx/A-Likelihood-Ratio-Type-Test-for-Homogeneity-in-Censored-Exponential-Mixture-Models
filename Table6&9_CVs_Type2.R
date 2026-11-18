#Simulated critical values for two real data applications 
M <- 500000
N <- 50
R <- round(N*0.9)
R2 <- round(N*0.8)
R3 <- round(N*0.7)
R4 <- round(N*0.6)
R5 <- round(N*0.5)

TSs_fun <- function(j){
  TSs_1 <- NULL
  TSs_2 <- NULL
  TSs_3 <- NULL
  TSs_4 <- NULL
  TSs_5 <- NULL
  TSs_6 <- NULL
  set.seed(j)
  obs <- sort(rexp(N))
  Hs <- NULL
  for (k in 1:(N-1)) {
    H <- -k*log(k) - (N-k)*log(N-k) + k*log(sum(obs[1:k])) + (N-k)*(log(sum(obs[(k+1):N])))
    Hs <- c(Hs, H)
  }
  Hmin <- min(Hs)
  TSs_1 <- - N*log(N) + N*log(sum(obs)) - Hmin
  obs_type2 <- c(obs[1:(R)], rep(obs[(R)], N - R))
  d <- c(rep(1, R), rep(0, N - R))
  Hs <- NULL
  for (k in 1:(R-1)) {
    H <- -k*log(k) - (R-k)*log(R-k) + k*log(sum(obs_type2[1:k])) + (R-k)*(log(sum(obs_type2[(k+1):N])))
    Hs <- c(Hs, H)
  }
  Hmax <- - (R-1)*log(R-1) + log(sum(obs_type2[(R):N])) + (R-1)*log(sum(obs_type2[1:(R-1)]))
  Hmin <- min(Hs, Hmax)
  TSs_2 <- - R*log(R) + R*log(sum(obs_type2)) - Hmin
  
  obs_type2_2 <- c(obs[1:(R2)], rep(obs[(R2)], N - R2))
  d <- c(rep(1, R2), rep(0, N - R2))
  Hs <- NULL
  for (k in 1:(R2-1)) {
    H <- -k*log(k) - (R2-k)*log(R2-k) + k*log(sum(obs_type2_2[1:k])) + (R2-k)*(log(sum(obs_type2_2[(k+1):N])))
    Hs <- c(Hs, H)
  }
  Hmax <- - (R2-1)*log(R2-1) + log(sum(obs_type2_2[(R2):N])) + (R2-1)*log(sum(obs_type2_2[1:(R2-1)]))
  Hmin <- min(Hs, Hmax)
  TSs_3 <- - R2*log(R2) + R2*log(sum(obs_type2_2)) - Hmin
  obs_type2_3 <- c(obs[1:(R3)], rep(obs[(R3)], N - R3))
  d <- c(rep(1, R3), rep(0, N - R3))
  Hs <- NULL
  for (k in 1:(R3-1)) {
    H <- -k*log(k) - (R3-k)*log(R3-k) + k*log(sum(obs_type2_3[1:k])) + (R3-k)*(log(sum(obs_type2_3[(k+1):N])))
    Hs <- c(Hs, H)
  }
  Hmax <- - (R3-1)*log(R3-1) + log(sum(obs_type2_3[(R3):N])) + (R3-1)*log(sum(obs_type2_3[1:(R3-1)]))
  Hmin <- min(Hs, Hmax)
  TSs_4 <- - R3*log(R3) + R3*log(sum(obs_type2_3)) - Hmin
  obs_type2_4 <- c(obs[1:(R4)], rep(obs[(R4)], N - R4))
  d <- c(rep(1, R4), rep(0, N - R4))
  Hs <- NULL
  for (k in 1:(R4-1)) {
    H <- -k*log(k) - (R4-k)*log(R4-k) + k*log(sum(obs_type2_4[1:k])) + (R4-k)*(log(sum(obs_type2_4[(k+1):N])))
    Hs <- c(Hs, H)
  }
  Hmax <- - (R4-1)*log(R4-1) + log(sum(obs_type2_4[(R4):N])) + (R4-1)*log(sum(obs_type2_4[1:(R4-1)]))
  Hmin <- min(Hs, Hmax)
  TSs_5 <- - R4*log(R4) + R4*log(sum(obs_type2_4)) - Hmin
  obs_type2_5 <- c(obs[1:(R5)], rep(obs[(R5)], N - R5))
  d <- c(rep(1, R5), rep(0, N - R5))
  Hs <- NULL
  for (k in 1:(R5-1)) {
    H <- -k*log(k) - (R5-k)*log(R5-k) + k*log(sum(obs_type2_5[1:k])) + (R5-k)*(log(sum(obs_type2_5[(k+1):N])))
    Hs <- c(Hs, H)
  }
  Hmax <- - (R5-1)*log(R5-1) + log(sum(obs_type2_5[(R5):N])) + (R5-1)*log(sum(obs_type2_5[1:(R5-1)]))
  Hmin <- min(Hs, Hmax)
  TSs_6 <- - R5*log(R5) + R5*log(sum(obs_type2_5)) - Hmin
  return(data.frame(TSs_1 = TSs_1, TSs_2 = TSs_2, TSs_3 = TSs_3,
                    TSs_4 = TSs_4, TSs_5 = TSs_5, TSs_6 = TSs_6))
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cl <- 24
registerDoParallel(cl)
#stop cluster
#stopCluster(cl)

tictoc::tic()
result_list <- foreach(i=1:M) %dopar% TSs_fun(i)
exectime <- tictoc::toc(quiet=TRUE)
print("foreach")
exectime$toc - exectime$tic

result_df <- do.call("rbind.data.frame", result_list)
any(result_df < 0)
save.image("Table6&9_CVs_Type2.RData")

alpha <- 0.05
apply(result_df, 2, quantile, 1-alpha)
