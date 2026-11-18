M <- 500000
N <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100)
R1 <- round(N*0.9)
R2 <- N*0.8
R3 <- round(N*0.7)
R4 <- N*0.6
R5 <- round(N*0.5)

TSs_fun <- function(ind){
  TSs <- NULL
  TScs <- NULL
  TSs_1 <- NULL
  TSs_2 <- NULL
  TSs_3 <- NULL
  TSs_4 <- NULL
  TSs_5 <- NULL
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
    TSs <- c(TSs, TS)
    
    obs_type2 <- c(obs[1:(R1[ind])], rep(obs[(R1[ind])], N[ind] - R1[ind]))
    d <- c(rep(1, R1[ind]), rep(0, N[ind] - R1[ind]))
    TSc <- - (N[ind] - R1[ind]) + R1[ind]*log(sum(obs_type2)/R1[ind]) - sum(log(obs[1:R1[ind]]))
    TScs <- c(TScs, TSc)
    Hs <- NULL
    for (k in 1:(R1[ind]-1)) {
      H <- -k*log(k) - (R1[ind]-k)*log(R1[ind]-k) + k*log(sum(obs_type2[1:k])) + (R1[ind]-k)*(log(sum(obs_type2[(k+1):N[ind]])))
      Hs <- c(Hs, H)
    }
    Hmax <- - (R1[ind]-1)*log(R1[ind]-1) + log(sum(obs_type2[(R1[ind]):N[ind]])) + (R1[ind]-1)*log(sum(obs_type2[1:(R1[ind]-1)]))
    Hmin <- min(Hs, Hmax)
    TS <- - R1[ind]*log(R1[ind]) + R1[ind]*log(sum(obs_type2)) - Hmin
    TSs_1 <- c(TSs_1, TS)
    obs_type2_2 <- c(obs[1:(R2[ind])], rep(obs[(R2[ind])], N[ind] - R2[ind]))
    d <- c(rep(1, R2[ind]), rep(0, N[ind] - R2[ind]))
    Hs <- NULL
    for (k in 1:(R2[ind]-1)) {
      H <- -k*log(k) - (R2[ind]-k)*log(R2[ind]-k) + k*log(sum(obs_type2_2[1:k])) + (R2[ind]-k)*(log(sum(obs_type2_2[(k+1):N[ind]])))
      Hs <- c(Hs, H)
    }
    Hmax <- - (R2[ind]-1)*log(R2[ind]-1) + log(sum(obs_type2_2[(R2[ind]):N[ind]])) + (R2[ind]-1)*log(sum(obs_type2_2[1:(R2[ind]-1)]))
    Hmin <- min(Hs, Hmax)
    TS <- - R2[ind]*log(R2[ind]) + R2[ind]*log(sum(obs_type2_2)) - Hmin
    TSs_2 <- c(TSs_2, TS)
    obs_type2_3 <- c(obs[1:(R3[ind])], rep(obs[(R3[ind])], N[ind] - R3[ind]))
    d <- c(rep(1, R3[ind]), rep(0, N[ind] - R3[ind]))
    Hs <- NULL
    for (k in 1:(R3[ind]-1)) {
      H <- -k*log(k) - (R3[ind]-k)*log(R3[ind]-k) + k*log(sum(obs_type2_3[1:k])) + (R3[ind]-k)*(log(sum(obs_type2_3[(k+1):N[ind]])))
      Hs <- c(Hs, H)
    }
    Hmax <- - (R3[ind]-1)*log(R3[ind]-1) + log(sum(obs_type2_3[(R3[ind]):N[ind]])) + (R3[ind]-1)*log(sum(obs_type2_3[1:(R3[ind]-1)]))
    Hmin <- min(Hs, Hmax)
    TS <- - R3[ind]*log(R3[ind]) + R3[ind]*log(sum(obs_type2_3)) - Hmin
    TSs_3 <- c(TSs_3, TS)
    obs_type2_4 <- c(obs[1:(R4[ind])], rep(obs[(R4[ind])], N[ind] - R4[ind]))
    d <- c(rep(1, R4[ind]), rep(0, N[ind] - R4[ind]))
    Hs <- NULL
    for (k in 1:(R4[ind]-1)) {
      H <- -k*log(k) - (R4[ind]-k)*log(R4[ind]-k) + k*log(sum(obs_type2_4[1:k])) + (R4[ind]-k)*(log(sum(obs_type2_4[(k+1):N[ind]])))
      Hs <- c(Hs, H)
    }
    Hmax <- - (R4[ind]-1)*log(R4[ind]-1) + log(sum(obs_type2_4[(R4[ind]):N[ind]])) + (R4[ind]-1)*log(sum(obs_type2_4[1:(R4[ind]-1)]))
    Hmin <- min(Hs, Hmax)
    TS <- - R4[ind]*log(R4[ind]) + R4[ind]*log(sum(obs_type2_4)) - Hmin
    TSs_4 <- c(TSs_4, TS)
    obs_type2_5 <- c(obs[1:(R5[ind])], rep(obs[(R5[ind])], N[ind] - R5[ind]))
    d <- c(rep(1, R5[ind]), rep(0, N[ind] - R5[ind]))
    Hs <- NULL
    for (k in 1:(R5[ind]-1)) {
      H <- -k*log(k) - (R5[ind]-k)*log(R5[ind]-k) + k*log(sum(obs_type2_5[1:k])) + (R5[ind]-k)*(log(sum(obs_type2_5[(k+1):N[ind]])))
      Hs <- c(Hs, H)
    }
    Hmax <- - (R5[ind]-1)*log(R5[ind]-1) + log(sum(obs_type2_5[(R5[ind]):N[ind]])) + (R5[ind]-1)*log(sum(obs_type2_5[1:(R5[ind]-1)]))
    Hmin <- min(Hs, Hmax)
    TS <- - R5[ind]*log(R5[ind]) + R5[ind]*log(sum(obs_type2_5)) - Hmin
    TSs_5 <- c(TSs_5, TS)
  }
  TSs_1 <- ifelse(TSs_1 < 0, 0, TSs_1)
  TSs_2 <- ifelse(TSs_2 < 0, 0, TSs_2)
  TSs_3 <- ifelse(TSs_3 < 0, 0, TSs_3)
  TSs_4 <- ifelse(TSs_4 < 0, 0, TSs_4)
  TSs_5 <- ifelse(TSs_5 < 0, 0, TSs_5)
  return(list(TSs = TSs, TScs = TScs, TSs_1 = TSs_1, TSs_2 = TSs_2, 
              TSs_3 = TSs_3, TSs_4 = TSs_4, TSs_5 = TSs_5))
}

parallel.function <- function(ind){
  result <- TSs_fun(ind = ind)
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cl <- 14
registerDoParallel(cl)
#stop cluster
#stopCluster(cl)

tictoc::tic()
result_list <- foreach(i=1:length(N)) %dopar% parallel.function(i)
exectime <- tictoc::toc(quiet=TRUE)

print("foreach")
exectime$toc - exectime$tic

save.image("Table2_CVs_Type2_n=10_100_6rs.RData")

alpha <- 0.05
TS0 <- sapply(result_list, '[', 'TSs')
CVs0 <- sapply(TS0, quantile, 1-alpha)

TS1 <- sapply(result_list, '[', 'TSs_1')
CVs1 <- sapply(TS1, quantile, 1-alpha)

TS2 <- sapply(result_list, '[', 'TSs_2') 
CVs2 <- sapply(TS2, quantile, 1-alpha)

TS3 <- sapply(result_list, '[', 'TSs_3') 
CVs3 <- sapply(TS3, quantile, 1-alpha)

TS4 <- sapply(result_list, '[', 'TSs_4') 
CVs4 <- sapply(TS4, quantile, 1-alpha)

TS5 <- sapply(result_list, '[', 'TSs_5') 
CVs5 <- sapply(TS5, quantile, 1-alpha)

AllCVs <- data.frame(n = N, CVs0 = CVs0, CVs1 = CVs1, CVs2 = CVs2, CVs3 = CVs3, CVs4 = CVs4, CVs5 = CVs5)
AllCVs
# write.csv(AllCVs, "AllCVs.csv", row.names = FALSE)





