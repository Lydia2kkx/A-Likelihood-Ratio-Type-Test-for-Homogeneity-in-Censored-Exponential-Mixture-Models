##########Power study for Weibull distribution & Exponential
N <- seq(50, 500, 50)
p <- c(0.1, 0.5, 0.9)
M <- 10000
shape <- c(seq(0.1, 0.9, 0.2), c(1, 2, 5, 10))
scale <- c(0.5, 2, 5, 10)

power_study <- function(N, p, R, shape, scale){
  TSs_type2 <- NULL
  N1s <- NULL
  for (j in 1:M) {
    set.seed(j)
    N1 <- rbinom(1, N, p)
    N2 <- N - N1
    set.seed(j)
    group1 <- rweibull(N1, shape = 1, scale = scale)
    set.seed(j)
    group2 <- rweibull(N2, shape = shape, scale = scale)
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
  y <- data.frame(Ind = 1:M, N = N, p = p, R = R, shape = shape, scale = scale, N1 = N1s, N2 = N-N1s, TestStats = TSs_type2)
  return(y)
}

parameter_starts <- expand.grid(N, p, shape, scale)
colnames(parameter_starts) <- c("N", "p", "shape", "scale")
parameter_starts$R <- parameter_starts$N*0.8
parameter_starts <- parameter_starts[, c("N", "p", "R", "shape", "scale")]

parallel.function <- function(i){
  result <- power_study(N = parameter_starts[i, 1], p = parameter_starts[i, 2], 
                        R = parameter_starts[i, 3], shape = parameter_starts[i, 4], 
                        scale = parameter_starts[i, 5])
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cl <- 24
registerDoParallel(cl)
#stop cluster
#stopCluster(cl)

tictoc::tic()
result_list <- foreach(i=1:nrow(parameter_starts)) %dopar% parallel.function(i)
exectime <- tictoc::toc(quiet=TRUE)
save.image("Figure8&9_Power_study_Weibull&Exp.RData")

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
  dat <- data.frame(N = unique(x$N), p = unique(x$p), shape = unique(x$shape), scale = unique(x$scale), Rejection_rate = y)
  return(dat)
}

reject_rate <- lapply(result_list, Rejection_rate)
reject_rate <- do.call("rbind.data.frame", reject_rate)
#write.csv(reject_rate, "reject_rate_Ns_nocensoring.csv")
reject_rate$group <- paste0(reject_rate$shape, reject_rate$scale1)
library(ggplot2)
reject_rate1 <- reject_rate[reject_rate$shape < 1, ]

colnames(reject_rate1)[1] <- "n"
colnames(reject_rate1)

custom_labeller1 <- setNames(paste0("p = ", unique(reject_rate1$p)), unique(reject_rate1$p))
custom_labeller2 <- setNames(paste0("theta == ", unique(reject_rate1$scale)), unique(reject_rate1$scale))

ggplot(reject_rate1, aes(x = n, y = Rejection_rate, color = as.factor(shape), group = group)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() +
  facet_grid(scale ~ p, labeller = labeller(p = custom_labeller1, scale = as_labeller(custom_labeller2, label_parsed))) +
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
  labs(y = "Power", color = "\u03B2 ") +
  ylim(0, 1) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))+ 
  scale_color_viridis_d()

ggsave(filename = "Figure8_Power_study_Weibull.pdf", width = 8, height = 8, units = "in", device = cairo_pdf)


reject_rate2 <- reject_rate[reject_rate$shape > 1,]

colnames(reject_rate2)[1] <- "n"
colnames(reject_rate2)

custom_labeller1 <- setNames(paste0("p = ", unique(reject_rate2$p)), unique(reject_rate2$p))
custom_labeller2 <- setNames(paste0("theta == ", unique(reject_rate2$scale)), unique(reject_rate2$scale))

ggplot(reject_rate2, aes(x = n, y = Rejection_rate, color = as.factor(shape))) + 
  geom_line() + 
  geom_point() + 
  theme_bw() +
  facet_grid(scale ~ p, labeller = labeller(p = custom_labeller1, scale = as_labeller(custom_labeller2, label_parsed))) +
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
  labs(y = "Power", color = "\u03B2 ") +
  ylim(0, 0.03)+ 
  scale_color_viridis_d()

ggsave(filename = "Figure9_Power_study_Weibull.pdf", width = 8, height = 8, units = "in", device = cairo_pdf)

