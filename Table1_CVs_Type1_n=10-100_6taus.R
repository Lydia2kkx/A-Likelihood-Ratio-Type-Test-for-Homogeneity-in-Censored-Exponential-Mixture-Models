M <- 500000
N <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100)
C1 <- -log(0.1)
C2 <- -log(0.2)
C3 <- -log(0.3)
C4 <- -log(0.4)
C5 <- -log(0.5)

TSs_fun <- function(x){
  TSs_1 <- NULL
  TS1s <- NULL
  TS2s <- NULL
  TS3s <- NULL
  TS4s <- NULL
  TS5s <- NULL
  seeds <- NULL
  j <- 1
  iter_counter <- 1
  while (j <= M) {
    set.seed(iter_counter)
    obs <- sort(rexp(x))
    if(!all(sum(obs < C1), sum(obs < C2), sum(obs < C3), sum(obs < C4), sum(obs < C5))){
      iter_counter <- iter_counter + 1
      next
    }else{
      seeds <- c(seeds, iter_counter)
      Hs <- NULL
      for (k in 1:(x-1)) {
        H <- -k*log(k) - (x-k)*log(x-k) + k*log(sum(obs[1:k])) + (x-k)*(log(sum(obs[(k+1):x])))
        Hs <- c(Hs, H)
      }
      Hmin <- min(Hs)
      TS <- - x*log(x) + x*log(sum(obs)) - Hmin
      TSs_1 <- c(TSs_1, TS)
      
      obs_censored1 <- apply(cbind(obs, C1), 1, min)
      d1 <- c(rep(1, sum(obs < C1)), rep(0, (x - sum(obs < C1))))
      m1 <- sum(d1)
      H1s <- NULL
      for (k in 1:(m1-1)) {
        H <- -k*log(k) - (m1-k)*log(m1-k) + k*log(sum(obs_censored1[1:k])) + (m1-k)*(log(sum(obs_censored1[(k+1):x])))
        H1s <- c(H1s, H)
      }
      H1min <- min(H1s)
      TS1 <- -m1*log(m1) + m1*log(sum(obs_censored1)) - H1min
      TS1s <- c(TS1s, TS1)
      
      obs_censored2 <- apply(cbind(obs, C2), 1, min)
      d2 <- c(rep(1, sum(obs < C2)), rep(0, (x - sum(obs < C2))))
      m2 <- sum(d2)
      H2s <- NULL
      for (k in 1:(m2-1)) {
        H <- -k*log(k) - (m2-k)*log(m2-k) + k*log(sum(obs_censored2[1:k])) + (m2-k)*(log(sum(obs_censored2[(k+1):x])))
        H2s <- c(H2s, H)
      }
      H2min <- min(H2s)
      TS2 <- -m2*log(m2) + m2*log(sum(obs_censored2)) - H2min
      TS2s <- c(TS2s, TS2)
      
      obs_censored3 <- apply(cbind(obs, C3), 1, min)
      d3 <- c(rep(1, sum(obs < C3)), rep(0, (x - sum(obs < C3))))
      m3 <- sum(d3)
      H3s <- NULL
      for (k in 1:(m3-1)) {
        H <- -k*log(k) - (m3-k)*log(m3-k) + k*log(sum(obs_censored3[1:k])) + (m3-k)*(log(sum(obs_censored3[(k+1):x])))
        H3s <- c(H3s, H)
      }
      H3min <- min(H3s)
      TS3 <- -m3*log(m3) + m3*log(sum(obs_censored3)) - H3min
      TS3s <- c(TS3s, TS3)
      
      obs_censored4 <- apply(cbind(obs, C4), 1, min)
      d4 <- c(rep(1, sum(obs < C4)), rep(0, (x - sum(obs < C4))))
      m4 <- sum(d4)
      H4s <- NULL
      for (k in 1:(m4-1)) {
        H <- -k*log(k) - (m4-k)*log(m4-k) + k*log(sum(obs_censored4[1:k])) + (m4-k)*(log(sum(obs_censored4[(k+1):x])))
        H4s <- c(H4s, H)
      }
      H4min <- min(H4s)
      TS4 <- -m4*log(m4) + m4*log(sum(obs_censored4)) - H4min
      TS4s <- c(TS4s, TS4)
      
      obs_censored5 <- apply(cbind(obs, C5), 1, min)
      d5 <- c(rep(1, sum(obs < C5)), rep(0, (x - sum(obs < C5))))
      m5 <- sum(d5)
      H5s <- NULL
      for (k in 1:(m5-1)) {
        H <- -k*log(k) - (m5-k)*log(m5-k) + k*log(sum(obs_censored5[1:k])) + (m5-k)*(log(sum(obs_censored5[(k+1):x])))
        H5s <- c(H5s, H)
      }
      H5min <- min(H5s)
      TS5 <- -m5*log(m5) + m5*log(sum(obs_censored5)) - H5min
      TS5s <- c(TS5s, TS5)
    }
    iter_counter <- iter_counter + 1
    j <- j + 1
  }
  
  TS1s <- ifelse(TS1s < 0, 0, TS1s)
  TS2s <- ifelse(TS2s < 0, 0, TS2s)
  TS3s <- ifelse(TS3s < 0, 0, TS3s)
  TS4s <- ifelse(TS4s < 0, 0, TS4s)
  TS5s <- ifelse(TS5s < 0, 0, TS5s)
  return(list(TSs_1 = TSs_1, TS1s = TS1s, TS2s = TS2s,
              TS3s = TS3s, TS4s = TS4s, TS5s = TS5s, seeds = seeds))
}

parallel.function <- function(ind){
  result <- TSs_fun(x = N[ind])
}

library(foreach)
library(doParallel)
cl <- 14
registerDoParallel(cl)

tictoc::tic()
result_list <- foreach(i=1:length(N)) %dopar% parallel.function(i)
exectime <- tictoc::toc(quiet=TRUE)

print("foreach")
exectime$toc - exectime$tic

save.image("Table1_CVs_Type1_n=10_100_6taus.RData")

alpha <- 0.05
TS0 <- sapply(result_list, '[', 'TSs_1')
CVs0 <- sapply(TS0, quantile, 1-alpha)

TS1 <- sapply(result_list, '[', 'TS1s')
CVs1 <- sapply(TS1, quantile, 1-alpha)

TS2 <- sapply(result_list, '[', 'TS2s') 
CVs2 <- sapply(TS2, quantile, 1-alpha)

TS3 <- sapply(result_list, '[', 'TS3s') 
CVs3 <- sapply(TS3, quantile, 1-alpha)

TS4 <- sapply(result_list, '[', 'TS4s') 
CVs4 <- sapply(TS4, quantile, 1-alpha)

TS5 <- sapply(result_list, '[', 'TS5s') 
CVs5 <- sapply(TS5, quantile, 1-alpha)

AllCVs <- data.frame(n = N, CVs0 = CVs0, CVs1 = CVs1, CVs2 = CVs2, CVs3 = CVs3, CVs4 = CVs4, CVs5 = CVs5)
AllCVs



