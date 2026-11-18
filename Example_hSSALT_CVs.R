#####Simulated critical values for Type-I censoring with N = 35, tau = 8 for the first stress level
M <- 500000
N <- 35
tau <- 8
C1 <- qexp(pexp(tau, rate = 1/exp(3.5)))

TSs_fun <- function(x, C1){
  TSs_1 <- NULL
  TS1s <- NULL
  seeds <- NULL
  seed <- 1
  j <- 1
  while (j <= M) {
    set.seed(seed)
    obs <- sort(rexp(x))
    obs_censored1 <- apply(cbind(obs, C1), 1, min)
    d1 <- c(rep(1, sum(obs <= C1)), rep(0, (x - sum(obs <= C1))))
    m1 <- sum(d1)
    if(m1 <= 1){
      seed <- seed + 1
      next
    }else{
      Hs <- NULL
      for (k in 1:(x-1)) {
        H <- -k*log(k) - (x-k)*log(x-k) + k*log(sum(obs[1:k])) + (x-k)*(log(sum(obs[(k+1):x])))
        Hs <- c(Hs, H)
      }
      Hmin <- min(Hs)
      TS <- - x*log(x) + x*log(sum(obs)) - Hmin
      TSs_1 <- c(TSs_1, TS)
      
      H1s <- NULL
      for (k in 1:(m1-1)) {
        H <- -k*log(k) - (m1-k)*log(m1-k) + k*log(sum(obs_censored1[1:k])) + (m1-k)*(log(sum(obs_censored1[(k+1):x])))
        H1s <- c(H1s, H)
      }
      H1min <- min(H1s)
      TS1 <- -m1*log(m1) + m1*log(sum(obs_censored1)) - H1min
      TS1s <- c(TS1s, TS1)
      seeds <- c(seeds, seed)
      seed <- seed + 1
      j <- j + 1
    }
  }
  TS1s <- ifelse(TS1s < 0, 0, TS1s)
  return(list(TSs_1 = TSs_1, TS1s = TS1s, Seeds = seeds))
}

tictoc::tic()
TSs_ls <- TSs_fun(x = N, C1 = C1)
exectime <- tictoc::toc(quiet = TRUE)
exectime$toc - exectime$tic

print("CV for hSSALT Type1")
save.image("Example_hSSALT_CVs_Type1.RData")

alpha <- 0.05
##No censoring
sapply(TSs_ls["TSs_1"], quantile, 1-alpha)
##With censoring
sapply(TSs_ls["TS1s"], quantile, 1-alpha)

#####Simulated critical values for Type-II censoring with N = 27, R = 22 for the second stress level
M <- 500000
N <- 27
R <- 22

TSs_fun <- function(ind){
  TSs_1 <- NULL
  TSs_2 <- NULL
  for (j in 1:M) {
    set.seed(j)
    obs <- sort(rexp(N[ind]))
    Hs <- NULL
    for (k in 1:(N[ind]-1)) {
      H <- -k*log(k) - (N[ind]-k)*log(N[ind]-k) + k*log(sum(obs[1:k])) + (N[ind]-k)*(log(sum(obs[(k+1):N[ind]])))
      Hs <- c(Hs, H)
    }
    Hmin <- min(Hs)
    TS <- - N[ind]*log(N[ind]) + N[ind]*log(sum(obs)) - Hmin
    TSs_1 <- c(TSs_1, TS)
    obs_type2 <- c(obs[1:(R[ind])], rep(obs[(R[ind])], N[ind] - R[ind]))
    d <- c(rep(1, R[ind]), rep(0, N[ind] - R[ind]))
    Hs <- NULL
    for (k in 1:(R[ind]-1)) {
      H <- -k*log(k) - (R[ind]-k)*log(R[ind]-k) + k*log(sum(obs_type2[1:k])) + (R[ind]-k)*(log(sum(obs_type2[(k+1):N[ind]])))
      Hs <- c(Hs, H)
    }
    Hmax <- - (R[ind]-1)*log(R[ind]-1) + log(sum(obs_type2[(R[ind]):N[ind]])) + (R[ind]-1)*log(sum(obs_type2[1:(R[ind]-1)]))
    Hmin <- min(Hs, Hmax)
    TS <- - R[ind]*log(R[ind]) + R[ind]*log(sum(obs_type2)) - Hmin
    TSs_2 <- c(TSs_2, TS)
  }
  TSs_1 <- ifelse(TSs_1 < 0, 0, TSs_1)
  TSs_2 <- ifelse(TSs_2 < 0, 0, TSs_2)
  return(list(TSs_1 = TSs_1, TSs_2 = TSs_2))
}

tictoc::tic()
TSs_ls <- TSs_fun(1)
exectime <- tictoc::toc(quiet = TRUE)
exectime$toc - exectime$tic

print("CV for hSSALT Type2")
save.image("Example_hSSALT_CVs_Type2.RData")

alpha <- 0.05
##No censoring
sapply(TSs_ls["TSs_1"], quantile, 1-alpha)
##With censoring
sapply(TSs_ls["TSs_2"], quantile, 1-alpha)
