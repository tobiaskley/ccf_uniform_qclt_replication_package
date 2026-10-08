# This script simulates Wasserstein distances for one configuration from cfg.
#
# Author: KLEY
###############################################################################

## (1) to run on cluster:
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(task_id)) {
  stop("SLURM_ARRAY_TASK_ID is not set.")
}


#  ### (2) to run locally
#  rm(list = ls(all=TRUE))
#  task_id <- 1

###############################################################################
# Start timing the complete task
###############################################################################

task_start_time <- proc.time()[["elapsed"]]


## (1)&(2) always run

library(parallel)

source("sim_params.R")


###############################################################################
# Detach parameters from configuration
###############################################################################

mdl           <- cfg[[task_id]]$mdl
schedule_id   <- cfg[[task_id]]$schedule_id
schedule_name <- cfg[[task_id]]$schedule_name
innov         <- cfg[[task_id]]$innov
innov_name    <- cfg[[task_id]]$innov_name
kappa2        <- cfg[[task_id]]$kappa2
kappa4        <- cfg[[task_id]]$kappa4
b_start       <- cfg[[task_id]]$b_start
b_end         <- cfg[[task_id]]$b_end

if (mdl == "AR1") {
  beta <- cfg[[task_id]]$beta
} else {
  beta <- NA_real_
}


###############################################################################
# Some checks
###############################################################################

if (length(k_s) != 1 || k_s[1] != 1) {
  stop("This simulation code is written for k = 1.")
}

k <- 1


###############################################################################
# Coefficient a_n
###############################################################################

get_a_n <- function(mdl, n, schedule_id, beta = NA_real_) {
  
  if (mdl == "MA1") {
    
    if (schedule_id == 1) {
      return(n^(-1/2))
    }
    
    if (schedule_id == 2) {
      return(0.5)
    }
    
    if (schedule_id == 3) {
      return(n^(1/4))
    }
    
    stop("Unknown MA(1) schedule.")
    
  }
  
  
  if (mdl == "AR1") {
    
    return(1 - n^(-beta))
    
  }
  
  
  stop("Unknown model.")
}


###############################################################################
# True autocovariance and asymptotic variance
###############################################################################

gamma_MA1 <- function(a, kappa2) {
  
  # lag k = 1
  return(kappa2 * a)
  
}


Sigma_MA1 <- function(a, kappa2, kappa4) {
  
  # asymptotic variance for lag k = 1
  return(kappa2^2 * (1 + 5*a^2 + a^4) + kappa4*a^2)
  
}


gamma_AR1 <- function(a, k, kappa2) {
  
  return(kappa2 * a^abs(k) / (1 - a^2))
  
}


Sigma_AR1 <- function(a, k, kappa2, kappa4) {
  
  k <- abs(k)
  
  term1 <- kappa2^2 *
    (1 + a^2 +
       a^(2*k) * (1 + a^2 + 2*k*(1 - a^2))) /
    (1 - a^2)^3
  
  term2 <- kappa4 * a^(2*k) / (1 - a^2)^2
  
  return(term1 + term2)
  
}


###############################################################################
# Simulation of R empirical lag-1 autocovariances
#
# The empirical autocovariance is
#
#   hat gamma(1) = (1/n) sum_{t=2}^n X_t X_{t-1},
#
# i.e. no demeaning and denominator n, exactly as in the Python simulation.
###############################################################################

###############################################################################
# Simulation of R empirical lag-1 autocovariances
#
# Common random numbers are used across the different n within a batch.
# For each fixed n, the R simulated processes still have the correct law.
###############################################################################

simulate_gammahat_MA1_batch <- function(n_vec, a_vec, R, innov) {
  
  K <- length(n_vec)
  
  # One length-R vector of accumulated autocovariances for each n
  gamma_sum <- lapply(seq_len(K), function(j) numeric(R))
  
  # Common innovations across the K different n-values
  eps_tm2 <- innov(R)
  eps_tm1 <- innov(R)
  
  max_n <- max(n_vec)
  
  for (t in 2:max_n) {
    
    eps_t <- innov(R)
    
    # Only processes whose sample size has not yet been reached are active
    active <- which(n_vec >= t)
    
    for (j in active) {
      
      a <- a_vec[j]
      
      X_tm1 <- eps_tm1 + a * eps_tm2
      X_t   <- eps_t   + a * eps_tm1
      
      gamma_sum[[j]] <- gamma_sum[[j]] + X_t * X_tm1
    }
    
    eps_tm2 <- eps_tm1
    eps_tm1 <- eps_t
  }
  
  for (j in seq_len(K)) {
    gamma_sum[[j]] <- gamma_sum[[j]] / n_vec[j]
  }
  
  return(gamma_sum)
}


simulate_gammahat_AR1_batch <- function(n_vec, a_vec, R, innov, warm_up) {
  
  K <- length(n_vec)
  
  # Current AR state for each n
  X <- lapply(seq_len(K), function(j) numeric(R))
  
  # Accumulated empirical autocovariance for each n
  gamma_sum <- lapply(seq_len(K), function(j) numeric(R))
  
  
  ###########################################################################
  # Common burn-in innovations
  ###########################################################################
  
  for (t in 1:warm_up) {
    
    eps_t <- innov(R)
    
    for (j in seq_len(K)) {
      X[[j]] <- a_vec[j] * X[[j]] + eps_t
    }
  }
  
  
  ###########################################################################
  # X_1
  ###########################################################################
  
  eps_t <- innov(R)
  
  for (j in seq_len(K)) {
    X[[j]] <- a_vec[j] * X[[j]] + eps_t
  }
  
  
  ###########################################################################
  # X_2, ..., X_max(n)
  ###########################################################################
  
  max_n <- max(n_vec)
  
  for (t in 2:max_n) {
    
    eps_t <- innov(R)
    
    active <- which(n_vec >= t)
    
    for (j in active) {
      
      X_new <- a_vec[j] * X[[j]] + eps_t
      
      gamma_sum[[j]] <-
        gamma_sum[[j]] + X_new * X[[j]]
      
      X[[j]] <- X_new
    }
  }
  
  
  for (j in seq_len(K)) {
    gamma_sum[[j]] <- gamma_sum[[j]] / n_vec[j]
  }
  
  return(gamma_sum)
}


###############################################################################
# Initialise L'Ecuyer's RNG streams
###############################################################################

RNGkind("L'Ecuyer-CMRG")
set.seed(12345)

# fetch seed for this task_id's substream
s <- .Random.seed

for (i in 1:task_id) {
  s <- nextRNGStream(s)
}

.Random.seed <- s


###############################################################################
# Standard-normal quantiles used to approximate W_1
###############################################################################

u <- (2*(1:R) - 1)/(2*R)
z <- qnorm(u)


###############################################################################
# Storage
###############################################################################

W1_gamma <- array(
  dim = c(b_end - b_start + 1,
          length(n_s),
          length(k_s))
)

KS_gamma <- array(
  dim = c(b_end - b_start + 1,
          length(n_s),
          length(k_s))
)

# These are useful for checking the simulations later
a_n_s     <- numeric(length(n_s))
gamma_n_s <- numeric(length(n_s))
Sigma_n_s <- numeric(length(n_s))


###############################################################################
# Precompute coefficients, true autocovariances and asymptotic variances
###############################################################################

for (i_n in seq_along(n_s)) {
  
  n <- n_s[i_n]
  
  a_n <- get_a_n(
    mdl         = mdl,
    n           = n,
    schedule_id = schedule_id,
    beta        = beta
  )
  
  a_n_s[i_n] <- a_n
  
  
  if (mdl == "MA1") {
    
    gamma_n_s[i_n] <- gamma_MA1(
      a      = a_n,
      kappa2 = kappa2
    )
    
    Sigma_n_s[i_n] <- Sigma_MA1(
      a      = a_n,
      kappa2 = kappa2,
      kappa4 = kappa4
    )
  }
  
  
  if (mdl == "AR1") {
    
    gamma_n_s[i_n] <- gamma_AR1(
      a      = a_n,
      k      = k,
      kappa2 = kappa2
    )
    
    Sigma_n_s[i_n] <- Sigma_AR1(
      a      = a_n,
      k      = k,
      kappa2 = kappa2,
      kappa4 = kappa4
    )
  }
}


###############################################################################
# Simulation
###############################################################################

# We process the 20 sample sizes in batches.
# Innovations are shared across n within a batch, but each n has its own a_n.
n_batch_size <- 10

n_batches <- split(
  seq_along(n_s),
  ceiling(seq_along(n_s) / n_batch_size)
)


###############################################################################
# Timing storage
###############################################################################

# Total runtime of each batch, including:
#   - simulation of gamma_hat
#   - standardisation
#   - sorting
#   - W1
#   - KS
#
# Rows correspond to b and columns to batches.
time_batch_sec <- matrix(
  NA_real_,
  nrow = b_end - b_start + 1,
  ncol = length(n_batches)
)


###############################################################################
# Main simulation loop
###############################################################################

for (b in b_start:b_end) {
  
  b_index <- b - b_start + 1
  
  for (i_batch in seq_along(n_batches)) {
    
    ind <- n_batches[[i_batch]]
    
    n_vec <- n_s[ind]
    a_vec <- a_n_s[ind]
    
    
    ###########################################################################
    # Start timing this batch
    ###########################################################################
    
    batch_start_time <- proc.time()[["elapsed"]]
    
    
    ###########################################################################
    # Simulate empirical autocovariances for all n in this batch
    ###########################################################################
    
    if (mdl == "MA1") {
      
      gamma_hat_batch <- simulate_gammahat_MA1_batch(
        n_vec = n_vec,
        a_vec = a_vec,
        R     = R,
        innov = innov
      )
      
    } else if (mdl == "AR1") {
      
      gamma_hat_batch <- simulate_gammahat_AR1_batch(
        n_vec   = n_vec,
        a_vec   = a_vec,
        R       = R,
        innov   = innov,
        warm_up = warm_up
      )
      
    } else {
      
      stop("Unknown model.")
      
    }
    
    
    ###########################################################################
    # Compute W1 and KS separately for each n in this batch
    ###########################################################################
    
    for (j in seq_along(ind)) {
      
      i_n <- ind[j]
      n   <- n_s[i_n]
      
      gamma_hat <- gamma_hat_batch[[j]]
      
      
      #########################################################################
      # Standardised statistic
      #########################################################################
      
      X <- sqrt(n) *
        (gamma_n_s[i_n] - gamma_hat) /
        sqrt(Sigma_n_s[i_n])
      
      X <- sort(X)
      
      
      #########################################################################
      # Wasserstein distance
      #########################################################################
      
      W1_gamma[b_index, i_n, 1] <-
        mean(abs(X - z))
      
      
      #########################################################################
      # Kolmogorov distance
      #########################################################################
      
      Phi_X <- pnorm(X)
      
      KS_gamma[b_index, i_n, 1] <- max(
        abs(Phi_X - (1:R) / R),
        abs(Phi_X - (0:(R - 1)) / R)
      )
      
      
      #########################################################################
      # Display result for this n
      #########################################################################
      
      cat(
        "task = ", task_id,
        ", model = ", mdl,
        ", innovation = ", innov_name,
        ", schedule = ", schedule_name,
        ", b = ", b,
        ", n = ", n,
        ", W1 = ", W1_gamma[b_index, i_n, 1],
        ", KS = ", KS_gamma[b_index, i_n, 1],
        "\n",
        sep = ""
      )
    }
    
    
    ###########################################################################
    # Runtime of this batch
    ###########################################################################
    
    time_batch_sec[b_index, i_batch] <-
      proc.time()[["elapsed"]] - batch_start_time
    
    cat(
      "Batch ", i_batch,
      " completed: n = ",
      paste(n_vec, collapse = ", "),
      "; time = ",
      round(time_batch_sec[b_index, i_batch], 2),
      " sec = ",
      round(time_batch_sec[b_index, i_batch] / 60, 3),
      " min\n\n",
      sep = ""
    )
    
    
    ###########################################################################
    # Release large simulation objects before next batch
    ###########################################################################
    
    rm(gamma_hat_batch, gamma_hat, X, Phi_X)
    gc(verbose = FALSE)
  }
}


###############################################################################
# Total runtime of this task
###############################################################################

time_task_sec <-
  proc.time()[["elapsed"]] - task_start_time

cat(
  "============================================================\n",
  "TASK ", task_id, " COMPLETED\n",
  "model = ", mdl,
  ", innovation = ", innov_name,
  ", schedule = ", schedule_name, "\n",
  "TOTAL TASK RUNTIME: ",
  round(time_task_sec, 2),
  " sec = ",
  round(time_task_sec / 60, 3),
  " min = ",
  round(time_task_sec / 3600, 3),
  " h\n",
  "============================================================\n",
  sep = ""
)


###############################################################################
# Save results
###############################################################################

save(
  W1_gamma,
  KS_gamma,
  mdl,
  schedule_id,
  schedule_name,
  beta,
  innov_name,
  kappa2,
  kappa4,
  b_start,
  b_end,
  a_n_s,
  gamma_n_s,
  Sigma_n_s,
  n_s,
  k_s,
  R,
  warm_up,
  n_batch_size,
  time_task_sec,
  time_batch_sec,
  file = paste("output/W1_", task_id, ".Rdata", sep = "")
)
