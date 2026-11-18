###Simulation for Figure 2
##########Power study
N <- seq(50, 500, 50)
p <- c(0.1, 0.5, 0.9)
M <- 10000
theta <- seq(2, 8, 2)

power_study <- function(N, p, R, theta){
  TSs_type2 <- NULL
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
    obs_type2 <- c(obs[1:R], rep(obs[R], N - R))
    d <- c(rep(1, R), rep(0, N - R))
    Hs <- NULL
    for (m in 1:(R-1)) {
      H <- -m*log(m) - (R-m)*log(R-m) + m*log(sum(obs_type2[1:m])) + (R-m)*(log(sum(obs_type2[(m+1):N])))
      Hs <- c(Hs, H)
    }
    Hmin <- min(Hs)
    TS <- - R*log(R) + R*log(sum(obs_type2)) - Hmin
    TSs_type2 <- c(TSs_type2, TS)
    N1s <- c(N1s, N1)
  }
  y <- data.frame(Ind = 1:M, N = N, p = p, R = R, theta = theta, N1 = N1s, N2 = N-N1s, TestStats = TSs_type2)
  return(y)
}

####No censoring
parameter_starts <- expand.grid(N, p, theta)
colnames(parameter_starts) <- c("N", "p", "theta")
parameter_starts$R <- parameter_starts$N
parameter_starts <- parameter_starts[, c("N", "p", "R", "theta")]

parallel.function <- function(i){
  result <- power_study(N = parameter_starts[i, 1], p = parameter_starts[i, 2], 
                        R = parameter_starts[i, 3], theta = parameter_starts[i, 4])
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cl <- 7
registerDoParallel(cl)
#stop cluster
#stopCluster(cl)

tictoc::tic()
result_list <- foreach(i=1:nrow(parameter_starts)) %dopar% parallel.function(i)
exectime <- tictoc::toc(quiet=TRUE)
print("foreach")
exectime$toc - exectime$tic

#Values 23.24018  42.28124  60.75567  79.01288  97.05645 114.95872 132.78205 150.50667 168.17881 185.78753 
#come from Figure2_CVs_Type2_n=50_500.R, output of sapply(TS0, quantile, 1-alpha)
Rejection_rate <- function(x){
  if(unique(x$N) == 50){
    y <- sum(x$TestStats > 23.24018)/M   ####Values are CVs
  }
  if(unique(x$N) == 100){
    y <- sum(x$TestStats > 42.28124)/M
  }
  if(unique(x$N) == 150){
    y <- sum(x$TestStats > 60.75567)/M
  }
  if(unique(x$N) == 200){
    y <- sum(x$TestStats > 79.01288)/M
  }
  if(unique(x$N) == 250){
    y <- sum(x$TestStats > 97.05645)/M
  }
  if(unique(x$N) == 300){
    y <- sum(x$TestStats > 114.95872)/M
  }
  if(unique(x$N) == 350){
    y <- sum(x$TestStats > 132.78205)/M
  }
  if(unique(x$N) == 400){
    y <- sum(x$TestStats > 150.50667)/M
  }
  if(unique(x$N) == 450){
    y <- sum(x$TestStats > 168.17881)/M
  }
  if(unique(x$N) == 500){
    y <- sum(x$TestStats > 185.78753)/M
  }
  dat <- data.frame(N = unique(x$N), p = unique(x$p), theta = unique(x$theta), Rejection_rate = y)
  return(dat)
}

reject_rate <- lapply(result_list, Rejection_rate)
reject_rate <- do.call("rbind.data.frame", reject_rate)
write.csv(reject_rate, "Figure2_Data_reject_rate_Ns_nocensoring.csv")

####10% censoring (90% Observed, R = 0.9N)
parameter_starts <- expand.grid(N, p, theta)
colnames(parameter_starts) <- c("N", "p", "theta")
parameter_starts$R <- parameter_starts$N*0.9
parameter_starts <- parameter_starts[, c("N", "p", "R", "theta")]

parallel.function <- function(i){
  result <- power_study(N = parameter_starts[i, 1], p = parameter_starts[i, 2], 
                        R = parameter_starts[i, 3], theta = parameter_starts[i, 4])
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cl <- 7
registerDoParallel(cl)
#stop cluster
#stopCluster(cl)

tictoc::tic()
result_list <- foreach(i=1:nrow(parameter_starts)) %dopar% parallel.function(i)
exectime <- tictoc::toc(quiet=TRUE)
print("foreach")
exectime$toc - exectime$tic

#Values 23.63557  43.13427  62.07424  80.79190  99.30405 117.66439 135.95359 154.12759 172.27693 190.35293 
#come from Figure2_CVs_Type2_n=50_500.R, output of sapply(TS1, quantile, 1-alpha)
Rejection_rate <- function(x){
  if(unique(x$N) == 50){
    y <- sum(x$TestStats > 23.63557)/M
  }
  if(unique(x$N) == 100){
    y <- sum(x$TestStats > 43.13427)/M
  }
  if(unique(x$N) == 150){
    y <- sum(x$TestStats > 62.07424)/M
  }
  if(unique(x$N) == 200){
    y <- sum(x$TestStats > 80.7919)/M
  }
  if(unique(x$N) == 250){
    y <- sum(x$TestStats > 99.30405)/M
  }
  if(unique(x$N) == 300){
    y <- sum(x$TestStats > 117.66439)/M
  }
  if(unique(x$N) == 350){
    y <- sum(x$TestStats > 135.95359)/M
  }
  if(unique(x$N) == 400){
    y <- sum(x$TestStats > 154.12759)/M
  }
  if(unique(x$N) == 450){
    y <- sum(x$TestStats > 172.27693)/M
  }
  if(unique(x$N) == 500){
    y <- sum(x$TestStats > 190.35293)/M
  }
  dat <- data.frame(N = unique(x$N), p = unique(x$p), theta = unique(x$theta), Rejection_rate = y)
  return(dat)
}

reject_rate <- lapply(result_list, Rejection_rate)
reject_rate <- do.call("rbind.data.frame", reject_rate)
write.csv(reject_rate, "Figure2_Data_reject_rate_Ns_p9.csv")

####20% censoring (80% Observed, R = 0.8N)
parameter_starts <- expand.grid(N, p, theta)
colnames(parameter_starts) <- c("N", "p", "theta")
parameter_starts$R <- parameter_starts$N*0.8
parameter_starts <- parameter_starts[, c("N", "p", "R", "theta")]

parallel.function <- function(i){
  result <- power_study(N = parameter_starts[i, 1], p = parameter_starts[i, 2], 
                        R = parameter_starts[i, 3], theta = parameter_starts[i, 4])
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cl <- 7
registerDoParallel(cl)
#stop cluster
#stopCluster(cl)

tictoc::tic()
result_list <- foreach(i=1:nrow(parameter_starts)) %dopar% parallel.function(i)
exectime <- tictoc::toc(quiet=TRUE)
print("foreach")
exectime$toc - exectime$tic

#Values 24.37308  45.01415  65.61357  86.43106 107.34887 128.34672 149.37031 170.40555 191.43157 212.40529 
#come from Figure2_CVs_Type2_n=50_500.R, output of sapply(TS2, quantile, 1-alpha)
Rejection_rate <- function(x){
  if(unique(x$N) == 50){
    y <- sum(x$TestStats > 24.37308)/M
  }
  if(unique(x$N) == 100){
    y <- sum(x$TestStats > 45.01415)/M
  }
  if(unique(x$N) == 150){
    y <- sum(x$TestStats > 65.61357)/M
  }
  if(unique(x$N) == 200){
    y <- sum(x$TestStats > 86.43106)/M
  }
  if(unique(x$N) == 250){
    y <- sum(x$TestStats > 107.34887)/M
  }
  if(unique(x$N) == 300){
    y <- sum(x$TestStats > 128.34672)/M
  }
  if(unique(x$N) == 350){
    y <- sum(x$TestStats > 149.37031)/M
  }
  if(unique(x$N) == 400){
    y <- sum(x$TestStats > 170.40555)/M
  }
  if(unique(x$N) == 450){
    y <- sum(x$TestStats > 191.43157)/M
  }
  if(unique(x$N) == 500){
    y <- sum(x$TestStats > 212.40529)/M
  }
  dat <- data.frame(N = unique(x$N), p = unique(x$p), theta = unique(x$theta), Rejection_rate = y)
  return(dat)
}

reject_rate <- lapply(result_list, Rejection_rate)
reject_rate <- do.call("rbind.data.frame", reject_rate)
write.csv(reject_rate, "Figure2_Data_reject_rate_Ns_p8.csv")

####30% censoring (70% Observed, R = 0.7N)
parameter_starts <- expand.grid(N, p, theta)
colnames(parameter_starts) <- c("N", "p", "theta")
parameter_starts$R <- parameter_starts$N*0.7
parameter_starts <- parameter_starts[, c("N", "p", "R", "theta")]

parallel.function <- function(i){
  result <- power_study(N = parameter_starts[i, 1], p = parameter_starts[i, 2], 
                        R = parameter_starts[i, 3], theta = parameter_starts[i, 4])
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cl <- 7
registerDoParallel(cl)
#stop cluster
#stopCluster(cl)

tictoc::tic()
result_list <- foreach(i=1:nrow(parameter_starts)) %dopar% parallel.function(i)
exectime <- tictoc::toc(quiet=TRUE)
print("foreach")
exectime$toc - exectime$tic

#Values 26.77391  52.86021  79.01113 105.12447 131.16301 157.17610 182.02250 209.10127 235.00510 260.90508 
#come from Figure2_CVs_Type2_n=50_500.R, output of sapply(TS3, quantile, 1-alpha)
Rejection_rate <- function(x){
  if(unique(x$N) == 50){
    y <- sum(x$TestStats > 26.77391)/M
  }
  if(unique(x$N) == 100){
    y <- sum(x$TestStats > 52.86021)/M
  }
  if(unique(x$N) == 150){
    y <- sum(x$TestStats > 79.01113)/M
  }
  if(unique(x$N) == 200){
    y <- sum(x$TestStats > 105.12447)/M
  }
  if(unique(x$N) == 250){
    y <- sum(x$TestStats > 131.16301)/M
  }
  if(unique(x$N) == 300){
    y <- sum(x$TestStats > 157.1761)/M
  }
  if(unique(x$N) == 350){
    y <- sum(x$TestStats > 182.0225)/M
  }
  if(unique(x$N) == 400){
    y <- sum(x$TestStats > 209.10127)/M
  }
  if(unique(x$N) == 450){
    y <- sum(x$TestStats > 235.0051)/M
  }
  if(unique(x$N) == 500){
    y <- sum(x$TestStats > 260.90508)/M
  }
  dat <- data.frame(N = unique(x$N), p = unique(x$p), theta = unique(x$theta), Rejection_rate = y)
  return(dat)
}

reject_rate <- lapply(result_list, Rejection_rate)
reject_rate <- do.call("rbind.data.frame", reject_rate)
write.csv(reject_rate, "Figure2_Data_reject_rate_Ns_p7.csv")

####40% censoring (60% Observed, R = 0.6N)
parameter_starts <- expand.grid(N, p, theta)
colnames(parameter_starts) <- c("N", "p", "theta")
parameter_starts$R <- parameter_starts$N*0.6
parameter_starts <- parameter_starts[, c("N", "p", "R", "theta")]

parallel.function <- function(i){
  result <- power_study(N = parameter_starts[i, 1], p = parameter_starts[i, 2], 
                        R = parameter_starts[i, 3], theta = parameter_starts[i, 4])
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cl <- 7
registerDoParallel(cl)
#stop cluster
#stopCluster(cl)

tictoc::tic()
result_list <- foreach(i=1:nrow(parameter_starts)) %dopar% parallel.function(i)
exectime <- tictoc::toc(quiet=TRUE)
print("foreach")
exectime$toc - exectime$tic

#Values 29.20504  58.42193  87.53809 116.57302 145.54423 174.53299 203.41812 232.26286 261.14656 289.97696 
#come from Figure2_CVs_Type2_n=50_500.R, output of sapply(TS4, quantile, 1-alpha)
Rejection_rate <- function(x){
  if(unique(x$N) == 50){
    y <- sum(x$TestStats > 29.20504)/M
  }
  if(unique(x$N) == 100){
    y <- sum(x$TestStats > 58.42193)/M
  }
  if(unique(x$N) == 150){
    y <- sum(x$TestStats > 87.53809)/M
  }
  if(unique(x$N) == 200){
    y <- sum(x$TestStats > 116.57302)/M
  }
  if(unique(x$N) == 250){
    y <- sum(x$TestStats > 145.54423)/M
  }
  if(unique(x$N) == 300){
    y <- sum(x$TestStats > 174.53299)/M
  }
  if(unique(x$N) == 350){
    y <- sum(x$TestStats > 203.41812)/M
  }
  if(unique(x$N) == 400){
    y <- sum(x$TestStats > 232.26286)/M
  }
  if(unique(x$N) == 450){
    y <- sum(x$TestStats > 261.14656)/M
  }
  if(unique(x$N) == 500){
    y <- sum(x$TestStats > 289.97696)/M
  }
  dat <- data.frame(N = unique(x$N), p = unique(x$p), theta = unique(x$theta), Rejection_rate = y)
  return(dat)
}

reject_rate <- lapply(result_list, Rejection_rate)
reject_rate <- do.call("rbind.data.frame", reject_rate)
write.csv(reject_rate, "Figure2_Data_reject_rate_Ns_p6.csv")

##############Plot Figure 2
reject_nocensoring <- read.csv("Figure2_Data_reject_rate_Ns_nocensoring.csv")
reject_nocensoring$CRate <- 0
reject_0.9 <- read.csv("Figure2_Data_reject_rate_Ns_p9.csv")
reject_0.9$CRate <- "10%"
reject_0.8 <- read.csv("Figure2_Data_reject_rate_Ns_p8.csv")
reject_0.8$CRate <- "20%"
reject_0.7 <- read.csv("Figure2_Data_reject_rate_Ns_p7.csv")
reject_0.7$CRate <- "30%"
reject_0.6 <- read.csv("Figure2_Data_reject_rate_Ns_p6.csv")
reject_0.6$CRate <- "40%"
reject_all <- rbind(reject_nocensoring, reject_0.9, reject_0.8, reject_0.7, reject_0.6)
reject_all <- reject_all[, -1]
colnames(reject_all)[1] <- "n"
colnames(reject_all)

custom_labeller1 <- setNames(paste0("p = ", unique(reject_all$p)), unique(reject_all$p))
custom_labeller2 <- setNames(paste0("theta == ", unique(reject_all$theta)), unique(reject_all$theta))

library(ggplot2)
ggplot(reject_all, aes(x = n, y = Rejection_rate, col = CRate)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() +
  facet_grid(theta ~ p, labeller = labeller(p = custom_labeller1, theta = as_labeller(custom_labeller2, label_parsed))) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",  
        plot.title = element_text(size=12, face="bold", hjust = 0.5),
        panel.grid.major = element_line(color = "gray", size = 0.2),  
        panel.grid.minor = element_blank(), 
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.width = grid::unit(1.5, "cm"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)
  ) +
  labs(y = "Power", color = "Censoring rate") + 
  scale_color_viridis_d()

ggsave(filename = "Figure2_Power_Ns_CRates.pdf", width = 8, height = 8, units = "in", device = cairo_pdf)
