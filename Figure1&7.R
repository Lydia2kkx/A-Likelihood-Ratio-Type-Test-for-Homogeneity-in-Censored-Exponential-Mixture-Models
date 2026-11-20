library(ggplot2)

dat_censoring <- read.csv("Figure1&7_reject_rate_censoring.csv")
dat_nocensoring <- read.csv("Figure1&7_reject_rate_nocensoring.csv")

###MCSE of the power estimates
M <- 10000
dat_censoring$MCSE <- sqrt(dat_censoring$Rejection_rate * (1 - dat_censoring$Rejection_rate) / M)
max(dat_censoring$MCSE)
dat_nocensoring$MCSE <- sqrt(dat_nocensoring$Rejection_rate * (1 - dat_nocensoring$Rejection_rate) / M)
max(dat_nocensoring$MCSE)

dat <- rbind(dat_censoring, dat_nocensoring)
colnames(dat)[1] <- "n"
colnames(dat)

custom_labeller1 <- setNames(paste0("n = ", unique(dat$n)), unique(dat$n))
custom_labeller2 <- setNames(paste0("p = ", unique(dat$p)), unique(dat$p))

###Figure for all theta2
ggplot(dat, aes(x = theta, y = Rejection_rate, col = Censoring)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(n ~ p, labeller = labeller(n = custom_labeller1, p = custom_labeller2)) +
  theme(
    strip.background = element_blank(),
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
  )  +
  # ggtitle(bquote("Proportion of rejections of" ~ H[0])) +
  labs(x = expression(theta[2]), y = "Power") +
  scale_color_viridis_d()

ggsave(filename = "Figure1_Power_study.pdf", width = 8, height = 8, units = "in", device = cairo_pdf)


###Figure for small theta2
dat_smalltheta <- dat[dat$theta <= 20, ]
ggplot(dat_smalltheta, aes(x = theta, y = Rejection_rate, col = Censoring)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(n ~ p, labeller = labeller(n = custom_labeller1, p = custom_labeller2)) +
  theme(
    strip.background = element_blank(),
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
  # ggtitle(bquote("Proportion of rejections of" ~ H[0] ~ "for small" ~ theta)) +
  labs(x = expression(theta[2]), y = "Power")+
  scale_color_viridis_d()

ggsave(filename = "Figure7_Power_study.pdf", width = 8, height = 8, units = "in", device = cairo_pdf)
