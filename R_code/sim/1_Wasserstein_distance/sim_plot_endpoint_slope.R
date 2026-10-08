###############################################################################
# Plot simulated Wasserstein distances against the rates from the paper
#
# Input:
#   output/W1.Rdata
#
# Output:
#   plots/ma1_gaussian.pdf
#   plots/ma1_t.pdf
#   plots/ar1_gaussian.pdf
#   plots/ar1_t.pdf
#
# and corresponding PNG files.
###############################################################################

rm(list = ls(all = TRUE))

library(ggplot2)
library(patchwork)


###############################################################################
# Load merged simulation results
###############################################################################

load("output/W1.Rdata")

if (!exists("W1")) {
  stop("The file output/W1.Rdata does not contain an object called W1.")
}

dir.create("plots", showWarnings = FALSE)


###############################################################################
# For robustness: make sure beta is available also for the MA(1) cases
#
# MA schedules:
#   schedule 1: a_n = n^(-1/2), beta = -1/2
#   schedule 2: a_n = 0.5,      beta = 0
#   schedule 3: a_n = n^(1/4),  beta = 1/4
###############################################################################

ind <- W1$mdl == "MA1" & (is.na(W1$beta) | !is.finite(W1$beta))

W1$beta[ind & W1$schedule_id == 1] <- -1/2
W1$beta[ind & W1$schedule_id == 2] <- 0
W1$beta[ind & W1$schedule_id == 3] <- 1/4


###############################################################################
# Average the independently estimated Wasserstein distances over b
#
# Each row of W1_plot corresponds to one:
#
#   model x innovation x schedule x n
#
# combination.
###############################################################################

W1_plot <- aggregate(
  W1 ~ mdl + schedule_id + schedule_name + beta + innov + n + k,
  data = W1,
  FUN = mean
)


###############################################################################
# Also compute Monte Carlo SD and SE across b.
#
# These are not used in the four figures below, but are useful diagnostics.
###############################################################################

W1_sd <- aggregate(
  W1 ~ mdl + schedule_id + schedule_name + beta + innov + n + k,
  data = W1,
  FUN = sd
)

names(W1_sd)[names(W1_sd) == "W1"] <- "W1_sd"


W1_B <- aggregate(
  W1 ~ mdl + schedule_id + schedule_name + beta + innov + n + k,
  data = W1,
  FUN = length
)

names(W1_B)[names(W1_B) == "W1"] <- "B"


W1_plot <- merge(
  W1_plot,
  W1_sd,
  by = c("mdl", "schedule_id", "schedule_name",
         "beta", "innov", "n", "k")
)

W1_plot <- merge(
  W1_plot,
  W1_B,
  by = c("mdl", "schedule_id", "schedule_name",
         "beta", "innov", "n", "k")
)

W1_plot$W1_se <- W1_plot$W1_sd / sqrt(W1_plot$B)


###############################################################################
# Rates from the revised paper
###############################################################################

# MA(1): the corrected standardized upper bound is O(n^(-1/2)),
# uniformly in theta. Use this reference rate for all three schedules.

rate_MA1 <- function(n, beta, schedule_id) {
  
  rate <- numeric(length(n))
  
  # theta_n = n^(-1/2): beta = -1/2
  ind_1 <- schedule_id == 1
  rate[ind_1] <- n[ind_1]^(-1/2)
  
  # fixed theta_n = 0.5: represented by beta = 0
  ind_2 <- schedule_id == 2
  rate[ind_2] <- n[ind_2]^(-1/2)
  
  # theta_n = n^(1/4): beta = 1/4
  ind_3 <- schedule_id == 3
  rate[ind_3] <- n[ind_3]^(-1/2)
  
  return(rate)
}


# ---------------------------------------------------------------------------
# AR(1)
#
# The paper gives
#
#   sup_{phi in A_n} W_n(P_phi)
#     = O( sqrt(log(n * nu_n)^3 / (n * nu_n)) ).
#
# In the simulations
#
#   nu_n = n^(-beta),
#
# hence
#
#   rate = sqrt(((1 - beta) * log(n))^3 / n^(1 - beta)).
# ---------------------------------------------------------------------------

rate_AR1 <- function(n, beta) {
  
  nv <- n^(1 - beta)
  sqrt(log(nv)^3 / nv)
  
}


###############################################################################
# Add the relevant rate to every row
###############################################################################

W1_plot$rate <- NA_real_


# MA(1)
ind <- W1_plot$mdl == "MA1"

W1_plot$rate[ind] <- rate_MA1(
  n           = W1_plot$n[ind],
  beta        = W1_plot$beta[ind],
  schedule_id = W1_plot$schedule_id[ind]
)


# AR(1)
ind <- W1_plot$mdl == "AR1"

W1_plot$rate[ind] <- rate_AR1(
  n    = W1_plot$n[ind],
  beta = W1_plot$beta[ind]
)


###############################################################################
# Checks
###############################################################################

if (any(!is.finite(W1_plot$rate))) {
  stop("At least one theoretical rate is not finite.")
}

if (any(W1_plot$rate <= 0)) {
  stop("At least one theoretical rate is non-positive.")
}

if (any(W1_plot$W1 <= 0)) {
  stop("At least one mean Wasserstein distance is non-positive.")
}


###############################################################################
# Estimate log-log slopes
#
# For each schedule separately, fit
#
#   log(mean W1) = intercept + slope * log(rate).
#
# A slope close to one means that the empirical W1 changes proportionally
# to the theoretical rate over the simulated range.
###############################################################################

get_slopes <- function(dat) {
  
  schedule_ids <- sort(unique(dat$schedule_id))
  
  slopes <- sapply(
    schedule_ids,
    function(j) {
      
      aux <- dat[dat$schedule_id == j, ]
      
      fit <- lm(log(W1) ~ log(rate), data = aux)
      
      unname(coef(fit)[2])
    }
  )
  
  names(slopes) <- schedule_ids
  
  return(slopes)
}


###############################################################################
# Plot MA(1)
###############################################################################

make_MA_plot <- function(dat) {
  
  ggplot(
    dat,
    aes(
      x = rate,
      y = W1,
      colour = factor(schedule_id),
      linetype = factor(schedule_id),
      shape = factor(schedule_id),
      group = factor(schedule_id)
    )
  ) +
    
    # Thin lines connecting the simulated points
    geom_line(
      linewidth = 0.45,
      alpha = 0.45
    ) +
    
    # Simulated Wasserstein distances
    geom_point(
      size = 1.8
    ) +
    
    scale_x_log10() +
    scale_y_log10() +
    
    scale_colour_manual(
      values = c(
        "1" = "orange",
        "2" = "green4",
        "3" = "blue"
      ),
      breaks = c("1", "2", "3"),
      labels = expression(
        theta[n] == n^{-1/2},
        theta[n] == 0.5,
        theta[n] == n^{1/4}
      ),
      name = NULL
    ) +
    
    scale_linetype_manual(
      values = c(
        "1" = "dashed",
        "2" = "solid",
        "3" = "dotted"
      ),
      breaks = c("1", "2", "3"),
      labels = expression(
        theta[n] == n^{-1/2},
        theta[n] == 0.5,
        theta[n] == n^{1/4}
      ),
      name = NULL
    ) +
    
    scale_shape_manual(
      values = c(
        "1" = 16,
        "2" = 17,
        "3" = 15
      ),
      breaks = c("1", "2", "3"),
      labels = expression(
        theta[n] == n^{-1/2},
        theta[n] == 0.5,
        theta[n] == n^{1/4}
      ),
      name = NULL
    ) +
    
    labs(
      x = expression(r[n]),
      y = expression(hat(W)[1])
    ) +
    
    theme_bw(base_size = 15) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}


###############################################################################
# Plot AR(1)
###############################################################################

make_AR_plot <- function(dat) {
  
  dat <- dat[order(dat$schedule_id, dat$n), ]
  
  ggplot(
    dat,
    aes(
      x = rate,
      y = W1,
      colour = factor(schedule_id),
      linetype = factor(schedule_id),
      shape = factor(schedule_id),
      group = factor(schedule_id)
    )
  ) +
    
    # Thin lines connecting the simulated points
    geom_path(
      linewidth = 0.45,
      alpha = 0.45
    ) +
    
    # Simulated Wasserstein distances
    geom_point(
      size = 1.8
    ) +
    
    scale_x_log10() +
    scale_y_log10() +
    
    scale_colour_manual(
      values = c(
        "1" = "orange",
        "2" = "green4",
        "3" = "blue"
      ),
      breaks = c("1", "2", "3"),
      labels = expression(
        nu[n] == n^{-0.1},
        nu[n] == n^{-0.2},
        nu[n] == n^{-0.3}
      ),
      name = NULL
    ) +
    
    scale_linetype_manual(
      values = c(
        "1" = "dashed",
        "2" = "solid",
        "3" = "dotted"
      ),
      breaks = c("1", "2", "3"),
      labels = expression(
        nu[n] == n^{-0.1},
        nu[n] == n^{-0.2},
        nu[n] == n^{-0.3}
      ),
      name = NULL
    ) +
    
    scale_shape_manual(
      values = c(
        "1" = 16,
        "2" = 17,
        "3" = 15
      ),
      breaks = c("1", "2", "3"),
      labels = expression(
        nu[n] == n^{-0.1},
        nu[n] == n^{-0.2},
        nu[n] == n^{-0.3}
      ),
      name = NULL
    ) +
    
    labs(
      x = expression(r[n]),
      y = expression(hat(W)[1])
    ) +
    
    theme_bw(base_size = 15) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}


###############################################################################
# Diagnostic: local-quadratic endpoint derivative estimates
# for different relative bandwidths
###############################################################################

endpoint_slope <- function(dat, q) {
  
  if (all(dat$mdl == "AR1")) {
    dat <- dat[dat$n^(1 - dat$beta) > exp(3), ]
  }
  
  dat <- dat[order(dat$n), ]
  
  x <- log(dat$rate)
  y <- log(dat$W1)
  
  i0 <- which.max(dat$n)
  x0 <- x[i0]
  
  z <- x - x0
  
  # Relative bandwidth in x-space
  h <- q * diff(range(x))
  
  w <- exp(
    -0.5 * (z / h)^2
  )
  
  fit <- lm(
    y ~ z + I(z^2),
    weights = w
  )
  
  unname(coef(fit)["z"])
}


endpoint_slopes_all <- function(dat, q_s) {
  
  res <- do.call(
    rbind,
    lapply(sort(unique(dat$schedule_id)), function(i_schedule) {
      
      dat_schedule <- dat[
        dat$schedule_id == i_schedule,
      ]
      
      slopes <- sapply(
        q_s,
        function(q) endpoint_slope(dat_schedule, q)
      )
      
      data.frame(
        schedule_id = i_schedule,
        bandwidth = q_s,
        inv_rel_bandwidth = 1 / q_s,
        slope = slopes
      )
    })
  )
  
  rownames(res) <- NULL
  res
}


###############################################################################
# Plot local-quadratic endpoint slope against inverse relative bandwidth
###############################################################################

make_endpoint_slope_plot <- function(
    dat,
    schedule_ids,
    y_breaks,
    y_limits,
    show_x = TRUE
) {
  
  aux <- dat[
    dat$schedule_id %in% schedule_ids,
  ]
  
  p <- ggplot(
    aux,
    aes(
      x = inv_rel_bandwidth,
      y = slope,
      colour = factor(schedule_id),
      linetype = factor(schedule_id),
      group = factor(schedule_id)
    )
  ) +
    
    geom_hline(
      yintercept = 1,
      colour = "black",
      linewidth = 0.6,
      alpha = 0.7
    ) +
    
    geom_line(
      linewidth = 0.8
    ) +
    
    scale_x_log10() +
    scale_y_continuous(
      breaks = y_breaks
    ) +
    
    coord_cartesian(
      ylim = y_limits
    ) +
    
    scale_colour_manual(
      values = c(
        "1" = "orange",
        "2" = "green4",
        "3" = "blue"
      ),
      name = NULL
    ) +
    
    scale_linetype_manual(
      values = c(
        "1" = "dashed",
        "2" = "solid",
        "3" = "dotted"
      ),
      name = NULL
    ) +
    
    labs(
      x = expression(q^{-1}),
      y = "estimated endpoint slope"
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  if (!show_x) {
    p <- p +
      labs(x = NULL) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  return(p)
}


###############################################################################
# Extract the four experiments
###############################################################################

dat_ma_gaussian <- W1_plot[
  W1_plot$mdl == "MA1" &
    W1_plot$innov == "N(0,1)",
]

dat_ma_t <- W1_plot[
  W1_plot$mdl == "MA1" &
    W1_plot$innov == "t_7",
]

dat_ar_gaussian <- W1_plot[
  W1_plot$mdl == "AR1" &
    W1_plot$innov == "N(0,1)",
]

dat_ar_t <- W1_plot[
  W1_plot$mdl == "AR1" &
    W1_plot$innov == "t_9",
]

q_s <- seq(0.25, 6, by = 0.1)

slopes_ma_gaussian <- endpoint_slopes_all(
  dat_ma_gaussian,
  q_s
)

slopes_ma_t <- endpoint_slopes_all(
  dat_ma_t,
  q_s
)

slopes_ar_gaussian <- endpoint_slopes_all(
  dat_ar_gaussian,
  q_s
)

slopes_ar_t <- endpoint_slopes_all(
  dat_ar_t,
  q_s
)

###############################################################################
# Check that all four experiments are present
###############################################################################

if (nrow(dat_ma_gaussian) == 0) {
  stop("No MA(1) Gaussian results found.")
}

if (nrow(dat_ma_t) == 0) {
  stop("No MA(1) t_7 results found.")
}

if (nrow(dat_ar_gaussian) == 0) {
  stop("No AR(1) Gaussian results found.")
}

if (nrow(dat_ar_t) == 0) {
  stop("No AR(1) t_9 results found.")
}


###############################################################################
# Create plots
###############################################################################

p_ma_gaussian <- make_MA_plot(dat_ma_gaussian)
p_ma_t        <- make_MA_plot(dat_ma_t)

p_ma_slope_gaussian <- make_endpoint_slope_plot(
  slopes_ma_gaussian,
  schedule_ids = c(1, 2, 3),
  y_breaks = c(1, 1.25, 1.5, 1.75, 2),
  y_limits = c(0.9, 2.2),
  show_x = TRUE
)

p_ma_slope_t <- make_endpoint_slope_plot(
  slopes_ma_t,
  schedule_ids = c(1, 2, 3),
  y_breaks = c(1, 1.25, 1.5, 1.75, 2),
  y_limits = c(0.9, 2.2),
  show_x = TRUE
)

p_ar_gaussian <- make_AR_plot(dat_ar_gaussian)
p_ar_t        <- make_AR_plot(dat_ar_t)

p_ar_slope_gaussian <- make_endpoint_slope_plot(
  slopes_ar_gaussian,
  schedule_ids = c(1, 2, 3),
  y_breaks = c(1, 1.25, 1.5, 1.75, 2),
  y_limits = c(0.95, 2.05),
  show_x = TRUE
)

p_ar_slope_t <- make_endpoint_slope_plot(
  slopes_ar_t,
  schedule_ids = c(1, 2, 3),
  y_breaks = c(1, 1.25, 1.5, 1.75, 2),
  y_limits = c(0.95, 2.05),
  show_x = TRUE
)

p_ma_gaussian_combined <-
  p_ma_gaussian /
  p_ma_slope_gaussian +
  plot_layout(heights = c(2, 1))

p_ma_t_combined <-
  p_ma_t /
  p_ma_slope_t +
  plot_layout(heights = c(2, 1))

p_ar_gaussian_combined <-
  p_ar_gaussian /
  p_ar_slope_gaussian +
  plot_layout(heights = c(2, 1))

p_ar_t_combined <-
  p_ar_t /
  p_ar_slope_t +
  plot_layout(heights = c(2, 1))

###############################################################################
# Save PDF versions
###############################################################################

ggsave(
  "plots/ma1_gaussian.pdf",
  p_ma_gaussian_combined,
  width = 7,
  height = 5 * 1.5
)

ggsave(
  "plots/ma1_t.pdf",
  p_ma_t_combined,
  width = 7,
  height = 5 * 1.5
)

ggsave(
  "plots/ar1_gaussian.pdf",
  plot = p_ar_gaussian_combined,
  width = 7,
  height = 5 * 1.5
)

ggsave(
  "plots/ar1_t.pdf",
  plot = p_ar_t_combined,
  width = 7,
  height = 5 * 1.5
)

