#Weibull Models
#Data Set 8.1 Complete Data: All 50 Items Put into Use at t = 0 and Failure Times Are in Weeks
RealData <- c(0.013, 0.065, 0.111, 0.111, 0.163, 0.309, 0.426, 0.535, 0.684, 0.747,
              0.997, 1.284, 1.304, 1.647, 1.829, 2.336, 2.838, 3.269, 3.977, 3.981,
              4.520, 4.789, 4.849, 5.202, 5.291, 5.349, 5.911, 6.018, 6.427, 6.456,
              6.572, 7.023, 7.087, 7.291, 7.787, 8.596, 9.388, 10.261, 10.713, 11.658,
              13.006, 13.388, 13.842, 17.152, 17.283, 19.418, 23.471, 24.777, 32.795, 48.105)
N <- length(RealData)
plot(RealData)

plot(density(RealData))

R <- round(N*0.9)
R2 <- round(N*0.8)
R3 <- round(N*0.7)
R4 <- round(N*0.6)
R5 <- round(N*0.5)
#####Full sample
obs <- sort(RealData)
Hs <- NULL
for (k in 1:(N-1)) {
  H <- -k*log(k) - (N-k)*log(N-k) + k*log(sum(obs[1:k])) + (N-k)*(log(sum(obs[(k+1):N])))
  Hs <- c(Hs, H)
}
Hmin <- min(Hs)
TS <- - N*log(N) + N*log(sum(obs)) - Hmin
TS

obs_type2 <- c(obs[1:(R)], rep(obs[(R)], N - R))
d <- c(rep(1, R), rep(0, N - R))
Hs <- NULL
for (k in 1:(R-1)) {
  H <- -k*log(k) - (R-k)*log(R-k) + k*log(sum(obs_type2[1:k])) + (R-k)*(log(sum(obs_type2[(k+1):N])))
  Hs <- c(Hs, H)
}
Hmax <- - (R-1)*log(R-1) + log(sum(obs_type2[(R):N])) + (R-1)*log(sum(obs_type2[1:(R-1)]))
Hmin <- min(Hs, Hmax)
TS_2 <- - R*log(R) + R*log(sum(obs_type2)) - Hmin  
TS_2

obs_type2_2 <- c(obs[1:(R2)], rep(obs[(R2)], N - R2))
d <- c(rep(1, R2), rep(0, N - R2))
Hs <- NULL
for (k in 1:(R2-1)) {
  H <- -k*log(k) - (R2-k)*log(R2-k) + k*log(sum(obs_type2_2[1:k])) + (R2-k)*(log(sum(obs_type2_2[(k+1):N])))
  Hs <- c(Hs, H)
}
Hmax <- - (R2-1)*log(R2-1) + log(sum(obs_type2_2[(R2):N])) + (R2-1)*log(sum(obs_type2_2[1:(R2-1)]))
Hmin <- min(Hs, Hmax)
TS_3 <- - R2*log(R2) + R2*log(sum(obs_type2_2)) - Hmin
TS_3

obs_type2_3 <- c(obs[1:(R3)], rep(obs[(R3)], N - R3))
d <- c(rep(1, R3), rep(0, N - R3))
Hs <- NULL
for (k in 1:(R3-1)) {
  H <- -k*log(k) - (R3-k)*log(R3-k) + k*log(sum(obs_type2_3[1:k])) + (R3-k)*(log(sum(obs_type2_3[(k+1):N])))
  Hs <- c(Hs, H)
}
Hmax <- - (R3-1)*log(R3-1) + log(sum(obs_type2_3[(R3):N])) + (R3-1)*log(sum(obs_type2_3[1:(R3-1)]))
Hmin <- min(Hs, Hmax)
TS_4 <- - R3*log(R3) + R3*log(sum(obs_type2_3)) - Hmin
TS_4

obs_type2_4 <- c(obs[1:(R4)], rep(obs[(R4)], N - R4))
d <- c(rep(1, R4), rep(0, N - R4))
Hs <- NULL
for (k in 1:(R4-1)) {
  H <- -k*log(k) - (R4-k)*log(R4-k) + k*log(sum(obs_type2_4[1:k])) + (R4-k)*(log(sum(obs_type2_4[(k+1):N])))
  Hs <- c(Hs, H)
}
Hmax <- - (R4-1)*log(R4-1) + log(sum(obs_type2_4[(R4):N])) + (R4-1)*log(sum(obs_type2_4[1:(R4-1)]))
Hmin <- min(Hs, Hmax)
TS_5 <- - R4*log(R4) + R4*log(sum(obs_type2_4)) - Hmin
TS_5

obs_type2_5 <- c(obs[1:(R5)], rep(obs[(R5)], N - R5))
d <- c(rep(1, R5), rep(0, N - R5))
Hs <- NULL
for (k in 1:(R5-1)) {
  H <- -k*log(k) - (R5-k)*log(R5-k) + k*log(sum(obs_type2_5[1:k])) + (R5-k)*(log(sum(obs_type2_5[(k+1):N])))
  Hs <- c(Hs, H)
}
Hmax <- - (R5-1)*log(R5-1) + log(sum(obs_type2_5[(R5):N])) + (R5-1)*log(sum(obs_type2_5[1:(R5-1)]))
Hmin <- min(Hs, Hmax)
TS_6 <- - R5*log(R5) + R5*log(sum(obs_type2_5)) - Hmin
TS_6

##############simulated p values
###The simulated CVs can be taken from Table 2
load("Table2_CVs_Type2_n=10_100_6rs.RData")
CVs_n50 <- result_list[[9]]
p0 <- (1 + sum(CVs_n50$TSs >= TS)) / (M + 1)

p1 <- (1 + sum(CVs_n50$TSs_1 >= TS_2)) / (M + 1)

p2 <- (1 + sum(CVs_n50$TSs_2 >= TS_3)) / (M + 1)

p3 <- (1 + sum(CVs_n50$TSs_3 >= TS_4)) / (M + 1)

p4 <- (1 + sum(CVs_n50$TSs_4 >= TS_5)) / (M + 1)

p5 <- (1 + sum(CVs_n50$TSs_5 >= TS_6)) / (M + 1)

p_values <- c(p0, p1, p2, p3, p4, p5)
p_values
