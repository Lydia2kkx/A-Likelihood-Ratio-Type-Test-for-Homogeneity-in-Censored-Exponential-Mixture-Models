##########Power study
N <- c(10, 100, 1000)
p <- c(0.1, 0.5, 0.9)
M <- 10000
theta <- c(1:20, seq(25, 100, 5))

power_study <- function(N, p, theta){
  TSs <- NULL
  N1s <- NULL
  for (j in 1:M) {
    set.seed(j)
    N1 <- rbinom(1, N, p)
    N2 <- N - N1    
    set.seed(j)
    group1 <- rexp(N1, rate = 1)
    set.seed(j)
    group2 <- rexp(N2, rate = 1/(theta))
    obs <- sort(c(group1, group2))
    Hs <- NULL
    for (m in 1:(N-1)) {
      H <- -m*log(m) - (N-m)*log(N-m) + m*log(sum(obs[1:m])) + (N-m)*(log(sum(obs[(m+1):N])))
      Hs <- c(Hs, H)
    }
    Hmin <- min(Hs)
    TS <- - N*log(N) + N*log(sum(obs)) - Hmin
    TSs <- c(TSs, TS)
    N1s <- c(N1s, N1)
  }
  y <- data.frame(Ind = 1:M, N = N, p = p, theta = theta, N1 = N1s, N2 = N-N1s, TestStats = TSs)
  return(y)
}

parameter_starts <- expand.grid(N, p, theta)
colnames(parameter_starts) <- c("N", "p", "theta")

parallel.function <- function(i){
  result <- power_study(N = parameter_starts[i, 1], p = parameter_starts[i, 2], 
                        theta = parameter_starts[i, 3])
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cl <- 48
registerDoParallel(cl)
#stop cluster
#stopCluster(cl)

tictoc::tic()
result_list <- foreach(i=1:nrow(parameter_starts)) %dopar% parallel.function(i)
exectime <- tictoc::toc(quiet=TRUE)
print("foreach")
exectime$toc - exectime$tic

save.image("Figure1&7_Power_study_nocensoring.RData")

#Values 6.854143  42.281236 360.221106 
#come from Figure1&7_CVs.R, output of sapply(TS1, quantile, 1-alpha)
Rejection_rate <- function(x){
  if(unique(x$N) == 10){
    y <- sum(x$TestStats > 6.854143)/M
  }else if(unique(x$N) == 100){
    y <- sum(x$TestStats > 42.281236)/M
  }else{
    y <- sum(x$TestStats > 360.221106)/M
  }
  dat <- data.frame(N = unique(x$N), p = unique(x$p), theta = unique(x$theta), Rejection_rate = y)
  return(dat)
}

reject_rate <- lapply(result_list, Rejection_rate)
reject_rate <- do.call("rbind.data.frame", reject_rate)
reject_rate$Censoring <- FALSE
write.csv(reject_rate, "Figure1&7_reject_rate_nocensoring.csv", row.names = FALSE)
