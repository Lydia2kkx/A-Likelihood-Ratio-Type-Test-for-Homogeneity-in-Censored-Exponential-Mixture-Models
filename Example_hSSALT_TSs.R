n <- 35
theta1 <- exp(3.5)
theta22 <- exp(2.0)
theta21 <- exp(-0.2)
tau <- 8

######Calculate the test statistic under the first stress level
T1 <- c(1.46, 2.22, 3.92, 4.24, 5.47, 5.60, 6.12, 6.56)
n1 <- length(T1)
T1_new <- c(T1, rep(tau, n-n1))
Hs <- NULL
for (k in 1:(n1-1)) {
  H <- -k*log(k) - (n1-k)*log(n1-k) + k*log(sum(T1_new[1:k])) + (n1-k)*(log(sum(T1_new[(k+1):(n-n1)])))
  Hs <- c(Hs, H)
}
Hmin <- min(Hs)
TS <- - n1*log(n1) + n1*log(sum(T1_new)) - Hmin
TS

######Version 1 (slightly different in TS)
######Calculate the test statistic under the second stress level
#####Below used the data that are rounded in table 3 
#####The calculated TS is the slightly different from the paper
T2 <- c(8.01, 8.10, 8.22, 8.59, 8.77, 8.80, 8.80, 8.84,
        8.90, 8.97, 9.62, 8.30, 8.98, 9.43, 9.87, 11.14, 11.85, 12.14, 13.49, 14.19, 14.33, 15.28)
R <- length(T2)
N <- n - n1
T2_new <- sort(T2 - tau)

obs_type2 <- c(T2_new, rep(T2_new[length(T2_new)], N - R))
d <- c(rep(1, R), rep(0, N - R))
Hs <- NULL
for (k in 1:(R-1)) {
  H <- -k*log(k) - (R-k)*log(R-k) + k*log(sum(obs_type2[1:k])) + (R-k)*(log(sum(obs_type2[(k+1):N])))
  Hs <- c(Hs, H)
}
TS_2 <- - R*log(R) + R*log(sum(obs_type2)) - min(Hs)  
TS_2


######Version 2
######Calculate the test statistic under the second stress level
#####Below used the data that are generated from the paper Lu & Kateri (2025), not the rounded one in table 3 
#####The calculated TS is the same as provided in the paper
T2 <- c(8.19, 8.30, 8.74, 8.98, 9.43, 9.87, 11.14, 11.76, 11.85, 12.14,
        14.04, 14.19, 14.24, 14.33, 15.28, 16.58, 16.85, 16.92, 17.80, 20.45,
        13.05, 13.49, 20.98, 21.09, 22.01, 26.34, 28.66)
n2 <- length(T2)
r <- 30
p <- 0.4
n21 <- round(p*(n-n1))
n22 <- n-n1-n21
set.seed(1234)
T22 <- sample(T2, n22)
set.seed(1234)
T21 <- sort(-theta21*log(1-runif(n21)) + tau)
T2_combine <- sort(c(T21, T22))
d <- c(rep(1, r-n1), rep(0, n-r))
cs_combine <- sort(T2_combine)[r - n1] 
t22_combine <- apply(cbind(T2_combine, cs_combine), 1, min) # observed or censored data
t22_combine_new <- t22_combine - tau
Hs <- NULL
for (k in 1:(r-n1-1)) {
  H <- -k*log(k) - (r-n1-k)*log(r-n1-k) + k*log(sum(t22_combine_new[1:k])) + (r-n1-k)*(log(sum(t22_combine_new[(k+1):(n-n1)])))
  Hs <- c(Hs, H)
}
Hmin <- min(Hs)
TS_2 <- - (r-n1)*log(r-n1) + (r-n1)*log(sum(t22_combine_new)) - Hmin
TS_2