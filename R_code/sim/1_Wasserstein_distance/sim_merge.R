# This script merges the results computed on the cluster
# into one data frame and saves it as W1.Rdata
#
# Author: KLEY
###############################################################################

rm(list = ls(all = TRUE))

source("sim_params.R")


###############################################################################
# Data frame for merged results
###############################################################################

W1 <- data.frame(
  mdl            = character(),
  cfg_id         = integer(),
  schedule_id    = integer(),
  schedule_name  = character(),
  beta           = numeric(),
  innov          = character(),
  b              = integer(),
  n              = integer(),
  k              = integer(),
  a_n            = numeric(),
  gamma_n        = numeric(),
  Sigma_n        = numeric(),
  W1             = numeric(),
  KS             = numeric(),
  batch_id       = integer(),
  time_batch_sec = numeric(),
  time_task_sec  = numeric(),
  stringsAsFactors = FALSE
)


###############################################################################
# Load and merge individual tasks
###############################################################################

files_completed <- 0

for (task_id in 1:length(cfg)) {
  
  file_name <- paste0("output/W1_", task_id, ".Rdata")
  
  if (file.exists(file_name)) {
    
    load(file_name)
    
    # exponent defining the coefficient schedule
    if (mdl == "MA1") {
      beta_out <- c(-1/2, 0, 1/4)[schedule_id]
    } else {
      beta_out <- beta
    }
    
    for (i_n in 1:length(n_s)) {
      
      batch_id <- ceiling(i_n / 10)
      
      W1 <- rbind(
        W1,
        data.frame(
          mdl            = mdl,
          cfg_id         = task_id,
          schedule_id    = schedule_id,
          schedule_name  = schedule_name,
          beta           = beta_out,
          innov          = innov_name,
          b              = b_start,
          n              = n_s[i_n],
          k              = k_s[1],
          a_n            = a_n_s[i_n],
          gamma_n        = gamma_n_s[i_n],
          Sigma_n        = Sigma_n_s[i_n],
          W1             = W1_gamma[1, i_n, 1],
          KS             = KS_gamma[1, i_n, 1],
          batch_id       = batch_id,
          time_batch_sec = time_batch_sec[1, batch_id],
          time_task_sec  = time_task_sec,
          stringsAsFactors = FALSE
        )
      )
      
    }
    
    cat("task id:", task_id, "- loaded", file_name, "\n")
    
    files_completed <- files_completed + 1
    
  } else {
    
    cat("task id:", task_id, "- FILE MISSING:", file_name, "\n")
    
  }
}


###############################################################################
# Summary and save
###############################################################################

cat("\n")
cat("Files completed:", files_completed, "of", length(cfg), "\n")
cat("Number of rows in merged data frame:", nrow(W1), "\n")

save(W1, file = "output/W1.Rdata")


###############################################################################
# Summary and save
###############################################################################


task_times <- unique(
  W1[, c(
    "cfg_id",
    "mdl",
    "innov",
    "schedule_id",
    "schedule_name",
    "b",
    "time_task_sec"
  )]
)

# Maximum task runtime
max(task_times$time_task_sec)

# In minutes
max(task_times$time_task_sec) / 60

# In hours
max(task_times$time_task_sec) / 3600
