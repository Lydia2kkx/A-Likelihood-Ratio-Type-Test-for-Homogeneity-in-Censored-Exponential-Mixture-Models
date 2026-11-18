M <- 500000
N <- c(10, 100, 1000)
R <- N*0.8

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

TSs_ls <- lapply(1:length(N), TSs_fun)

save.image("Figure1&7_CVs.RData")

alpha <- 0.05
#no censoring
TS1 <- sapply(TSs_ls, '[', 'TSs_1')
sapply(TS1, quantile, 1-alpha)

#r = 0.8n Type-II censoring
TS2 <- sapply(TSs_ls, '[', 'TSs_2') 
sapply(TS2, quantile, 1-alpha)


