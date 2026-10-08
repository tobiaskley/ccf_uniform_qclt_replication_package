# This contains the parameters with which sim.R is run.
# In particular the list cfg 
# 
# Author: KLEY
###############################################################################

# how many independent replications per estimation of W1
R <- 10^7

# how many independent estimations of W1
B <- 50

warm_up <- 500

n_s <- c(10, 13, 19, 26, 37, 51, 71, 98, 136, 189,
         263, 365, 506, 702, 974, 1351, 1874, 2599, 3605, 5000)

k_s <- 1

innov_s <- list(function(n) {rnorm(n, 0, 1)},
                function(n) {sqrt(5/7) * rt(n, df = 7)},
                function(n) {sqrt(7/9) * rt(n, df = 9)}
)

innov_name_s <- c("N(0,1)", "t_7", "t_9")
kappa2_s <- c(1, 1, 1)
kappa4_s <- c(0, 2, 6/5)

# MA(1) coefficient schedules:
#   1: a_n = n^(-1/2)
#   2: a_n = 0.5
#   3: a_n = n^(1/4)
ma_schedule_name_s <- c("n^-1/2", "0.5", "n^1/4")


# AR(1) coefficient schedules:
#   a_n = 1 - v_n,  v_n = n^(-beta)
ar_beta_s <- c(0.1, 0.2, 0.3)
ar_schedule_name_s <- c("n^-0.1", "n^-0.2", "n^-0.3")

cfg <- list()

for (b in 1:B) {
  
  ############################################################
  # MA(1)
  #
  # innovations:
  #   i_innov = 1: Gaussian
  #   i_innov = 2: standardized t_7
  ############################################################
  
  for (i_innov in c(1, 2)) {
    
    for (i_schedule in 1:length(ma_schedule_name_s)) {
      
      cfg <- c(cfg, list(list(
        mdl           = "MA1",
        schedule_id   = i_schedule,
        schedule_name = ma_schedule_name_s[i_schedule],
        innov         = innov_s[[i_innov]],
        innov_name    = innov_name_s[i_innov],
        kappa2        = kappa2_s[i_innov],
        kappa4        = kappa4_s[i_innov],
        b_start       = b,
        b_end         = b
      )))
      
    }
  }
  
  
  ############################################################
  # AR(1)
  #
  # innovations:
  #   i_innov = 1: Gaussian
  #   i_innov = 3: standardized t_9
  ############################################################
  
  for (i_innov in c(1, 3)) {
    
    for (i_schedule in 1:length(ar_beta_s)) {
      
      cfg <- c(cfg, list(list(
        mdl           = "AR1",
        schedule_id   = i_schedule,
        schedule_name = ar_schedule_name_s[i_schedule],
        beta          = ar_beta_s[i_schedule],
        innov         = innov_s[[i_innov]],
        innov_name    = innov_name_s[i_innov],
        kappa2        = kappa2_s[i_innov],
        kappa4        = kappa4_s[i_innov],
        b_start       = b,
        b_end         = b
      )))
      
    }
  }
  
}
