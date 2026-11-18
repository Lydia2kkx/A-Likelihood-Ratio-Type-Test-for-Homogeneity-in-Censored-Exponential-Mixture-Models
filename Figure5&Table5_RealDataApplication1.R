# Load required packages
library(ggplot2)
library(dplyr)
library(mixtools)
# RealData and model parameters
RealData <- c(0.013, 0.065, 0.111, 0.111, 0.163, 0.309, 0.426, 0.535, 0.684, 0.747,
              0.997, 1.284, 1.304, 1.647, 1.829, 2.336, 2.838, 3.269, 3.977, 3.981,
              4.520, 4.789, 4.849, 5.202, 5.291, 5.349, 5.911, 6.018, 6.427, 6.456,
              6.572, 7.023, 7.087, 7.291, 7.787, 8.596, 9.388, 10.261, 10.713, 11.658,
              13.006, 13.388, 13.842, 17.152, 17.283, 19.418, 23.471, 24.777, 32.795, 48.105)
# Model parameters
d <- 1
####Multiple initial values for the EM fit
lambda <- seq(0.1, 0.9, 0.1)
rate1 <- 1 / seq(0.1, 1, 0.1)
rate2 <- 1 / seq(2, 20, 2) 
em_grid <- expand.grid(
  lambda = lambda, 
  rate1  = rate1, 
  rate2  = rate2
)

fits <- vector("list", nrow(em_grid))
logliks <- rep(NA, nrow(em_grid))

for (i in seq_len(nrow(em_grid))) {
  lambda <- em_grid$lambda[i]
  rate1  <- em_grid$rate1[i]
  rate2  <- em_grid$rate2[i]
  
  fit <- try(expRMM_EM(RealData, d, lambda = c(lambda, 1 - lambda), rate   = c(rate1, rate2), k = 2), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    fits[[i]] <- fit
    logliks[i] <- fit$loglik
  }
}

####Select the one with the largest log-likelihood
best_idx <- which.max(logliks)
best_model_em <- fits[[best_idx]]
best_model_em$lambda
1/best_model_em$rate

expmixture <- function(x, rate, p) {
  p[1]*pexp(x, rate[1]) + p[2]*pexp(x, rate[2])
}
ks.test(RealData, expmixture, rate = best_model_em$rate, p = best_model_em$lambda)

library(fitdistrplus)
fit.weibull <- fitdist(RealData, "weibull")
fit.weibull
fit.weibull$loglik
ks.test(RealData, "pweibull", shape = fit.weibull$estimate[1], scale = fit.weibull$estimate[2])
plot(fit.weibull)

fit.exp <- fitdist(RealData, "exp")
fit.exp
fit.exp$loglik
ks.test(RealData, "pexp", rate = fit.exp$estimate)
plot(fit.exp)

fit.weibull.aic <- 2*2-2*fit.weibull$loglik
fit.weibull.aic
fit.weibull.bic <- log(length(RealData))*2-2*fit.weibull$loglik
fit.weibull.bic

fit.exp.aic <- 2*1-2*fit.exp$loglik
fit.exp.aic
fit.exp.bic <- log(length(RealData))*1-2*fit.exp$loglik
fit.exp.bic

fit.exp.mixture.aic <- 2*3-2*best_model_em$loglik
fit.exp.mixture.aic
fit.exp.mixture.bic <- log(length(RealData))*3-2*best_model_em$loglik
fit.exp.mixture.bic

# Create a data frame for the RealData
df <- data.frame(x = RealData)

# Create a sequence of x values for the density curves
x_seq <- seq(0, max(RealData), by = 0.001)

# Calculate the densities based on the model parameters
mix_curve <- best_model_em$lambda[1] * dexp(x_seq, best_model_em$rate[1]) + 
  best_model_em$lambda[2] * dexp(x_seq, best_model_em$rate[2])
weibull_density <- dweibull(x_seq, shape = fit.weibull$estimate[1], scale = fit.weibull$estimate[2])
exp_density <- dexp(x_seq, rate = fit.exp$estimate)

# Create a data frame for the density curves
df_curves <- data.frame(
  x = x_seq,
  mix_curve = mix_curve,
  weibull_density = weibull_density,
  exp_density = exp_density
)

# Figure 5 (left)
ggplot(df, aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 15, fill = NA, color = "black") +
  geom_line(data = df_curves, aes(x = x, y = mix_curve, color = "Exponential Mixture"), size = 1.2) +
  geom_line(data = df_curves, aes(x = x, y = weibull_density, color = "Weibull"), size = 1.2) +
  geom_line(data = df_curves, aes(x = x, y = exp_density, color = "Exponential"), size = 1.2) +
  scale_color_manual(values = c(
    "Exponential" = "#440154FF", "Weibull" = "#21908CFF", "Exponential Mixture" = "#FDE725FF")) +
  # ggtitle("Histogram and fitted densities for three models") +
  labs(x = "Data", y = "Density", color = "Models") +
  theme(
    plot.background = element_rect(fill = NA, color = NA),  
    panel.background = element_rect(fill = "white", color = NA),    
    panel.grid.major = element_blank(),                             
    panel.grid.minor = element_blank(),                             
    axis.line = element_line(color = "black"),                      
    plot.title = element_text(hjust = 0.5),          
    legend.position = c(0.95, 0.95),                                   
    legend.justification = c(1, 1),                                 
    legend.box.background = element_blank(), 
    legend.box.spacing = unit(0.2, "cm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.width = grid::unit(1.5, "cm"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12)
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.1))

ggsave(filename = "Figure5_Subplot1.pdf", width = 6, height = 5, units = "in", device = cairo_pdf)

# Figure 5 (right) Data preparation
n <- length(RealData)
censor_levels <- c(0, 0.1, 0.2)

grid <- seq(min(RealData), max(RealData), length.out = 300)

km_all <- data.frame()
censor_points <- data.frame()
plot_data <- data.frame()

for (c in censor_levels) {
  r <- floor(n * (1 - c))
  Fmax <- r / n
  ordered_data <- sort(RealData)
  observed <- ordered_data[1:r]
  censored_true <- ordered_data[(r + 1):n]
  censor_time <- max(ordered_data[r], ordered_data[r + 1])
  time <- c(observed, rep(censor_time, n - r))
  status <- c(rep(1, r), rep(0, n - r))
  
  # ---------------- KM ----------------
  km_fit <- survfit(Surv(time, status) ~ 1)
  km_df <- data.frame(Failure = km_fit$time, km_cdf = pmin(1 - km_fit$surv, Fmax),
                      censor = paste0(round(c * 100), "%"), Type = "Observed")
  km_all <- rbind(km_all, km_df)
  
  # ---------------- censored TRUE values ----------------
  if (n - r > 0) {
    cens_df <- data.frame(Failure = censored_true, km_cdf = rep(Fmax, n - r),
                          censor = paste0(round(c * 100), "%"), Type = "Censored")
    censor_points <- rbind(censor_points, cens_df)
  }
  
  # ---------------- exponential ----------------
  lambda_exp <- 1 / mean(time)
  exp_cdf <- pexp(grid, rate = lambda_exp)
  
  # ---------------- Weibull ----------------
  fit.weib <- fitdist(time, "weibull")
  weib_cdf <- pweibull(grid, shape = fit.weib$estimate["shape"], scale = fit.weib$estimate["scale"])
  
  # ---------------- EM mixture ----------------
  em_results <- lapply(seq_len(nrow(em_grid)), function(i) {
    lambda <- em_grid$lambda[i]
    rate1  <- em_grid$rate1[i]
    rate2  <- em_grid$rate2[i]
    fit <- tryCatch(
      expRMM_EM(c(observed, rep(censor_time, n - r)), d = c(rep(1, r), rep(0, n - r)),
                lambda = c(lambda, 1 - lambda), rate = c(rate1, rate2), k = 2), 
      error = function(e) NULL)
    
    if (!is.null(fit) && length(fit$loglik) > 0) {
      list(fit = fit, loglik = tail(fit$loglik, 1))
    } else {
      NULL
    }
  })
  
  em_results <- Filter(Negate(is.null), em_results)
  
  if (length(em_results) > 0) {
    best_fit <- em_results[[which.max(sapply(em_results, `[[`, "loglik"))]]$fit
    mix_cdf <- best_fit$lambda[1] * pexp(grid, best_fit$rate[1]) +
      best_fit$lambda[2] * pexp(grid, best_fit$rate[2])
  } else {
    mix_cdf <- rep(NA, length(grid))
  }
  tmp <- data.frame(Failure = grid, exp_curve = exp_cdf, weib_curve = weib_cdf,
                    mix_curve = mix_cdf, censor = paste0(round(c * 100), "%"))
  plot_data <- rbind(plot_data, tmp)
}

km_fit1 <- survfit(Surv(RealData, rep(1, n)) ~ 1)
plot_points <- data.frame(
  Failure = unique(RealData),
  km_cdf = 1 - km_fit1$surv,
  censor = paste0(0, "%"),
  Type = "Observed"
)
plot_points <- rbind(plot_points, data.frame(Failure = unique(RealData),
                                             km_cdf = 1 - km_fit1$surv,
                                             censor = paste0(10, "%"),
                                             Type = c(rep("Observed", 44), rep("Censored", 5))))
plot_points <- rbind(plot_points, data.frame(Failure = unique(RealData),
                                             km_cdf = 1 - km_fit1$surv,
                                             censor = paste0(20, "%"),
                                             Type = c(rep("Observed", 39), rep("Censored", 10))))
plot_points$Type <- factor(plot_points$Type, levels = c("Observed", "Censored"))

# Figure 5 (right)
ggplot() +
  # observed failures
  geom_point(data = plot_points, aes(x = Failure, y = km_cdf, color = Type, shape = Type),
             size = 2) +
  # censoring threshold line
  geom_vline(data = data.frame(censor = paste0(round(censor_levels * 100), "%"),
                               xint = sapply(censor_levels, function(c) {
                                 r <- floor(n * (1 - c))
                                 max(sort(RealData)[r], sort(RealData)[r + 1])
                                 })), aes(xintercept = xint), linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") + 
  # exponential
  geom_line(data = plot_data, aes(x = Failure, y = exp_curve, color = "Exponential"),
            linewidth = 1) +
  # weibull
  geom_line(data = plot_data, aes(x = Failure, y = weib_curve, color = "Weibull"),
            linewidth = 1) +
  # mixture
  geom_line(data = plot_data, aes(x = Failure, y = mix_curve, color = "Exponential Mixture"),
            linewidth = 1) +
  facet_wrap(~ censor, scales = "free_x") +
  labs(x = "Failure Time", y = "CDF") +
  scale_color_manual(values = c("Observed" = "black", "Censored" = "gray60", "Exponential" = "#440154FF",
                                "Exponential Mixture" = "#FDE725FF", "Weibull" = "#21908CFF"),
                     breaks = c("Observed", "Censored", "Exponential", "Weibull", "Exponential Mixture")) +
  scale_shape_manual(values = c("Observed" = 16, "Censored" = 17), guide = "none") +
  guides(color = guide_legend(override.aes = list(shape = c(16, 17, NA, NA, NA),
                                                  linetype = c("blank", "blank", "solid", "solid", "solid")))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.width = grid::unit(1.5, "cm"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))

ggsave(filename = "Figure5_Subplot2.pdf", width = 11, height = 5, units = "in", device = cairo_pdf)
