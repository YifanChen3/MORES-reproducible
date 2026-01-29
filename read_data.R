# Set the directory containing your RData files
  folder_path <- "/Users/chenyifan/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/our_cv_final_low/laplace_marginal/correlated"

  # Get the list of RData files
  rdata_files <- list.files(folder_path, pattern = "\\.RData$", full.names = TRUE)

  # Initialize lists to store extracted variables
  fnormB_values <- list()
  f1B_values <- list()
  fnormbeta_values <- list()
  f1beta_values <- list()
  time_values <- list()

  # Load each RData file and extract Fnorm and Onorm
  for (file in rdata_files) {
    env <- new.env()
    load(file, envir = env)

    # Convert the environment to a list
    data_list <- as.list(env)

    # Extract variables if they exist
    fnormB_values <- c(fnormB_values, data_list$result.sim$fnorm.B)
    f1B_values <- c(f1B_values, data_list$result.sim$F1.B)
    fnormbeta_values <- c(fnormbeta_values, data_list$result.sim$fnorm.beta)
    f1beta_values <- c(f1beta_values, data_list$result.sim$F1.beta)
    time_values <- c(time_values,data_list$result.sim$time)

  }

  # Convert lists to numeric vectors
  fnormB_values <- unlist(fnormB_values)
  f1B_values <- unlist(f1B_values)
  fnormbeta_values <- unlist(fnormbeta_values)
  f1beta_values <- unlist(f1beta_values)
  time_values <- unlist(time_values)

  # Compute mean and variance
  fnormB_mean <- mean(fnormB_values, na.rm = TRUE)
  fnormB_var <- sd(fnormB_values, na.rm = TRUE)

  f1B_mean <- mean(f1B_values, na.rm = TRUE)
  f1B_var <- sd(f1B_values, na.rm = TRUE)

  fnormbeta_mean <- mean(fnormbeta_values, na.rm = TRUE)
  fnormbeta_var <- sd(fnormbeta_values, na.rm = TRUE)

  f1beta_mean <- mean(f1beta_values, na.rm = TRUE)
  f1beta_var <- sd(f1beta_values, na.rm = TRUE)

  time_mean <- mean(time_values,na.rm = T)
  time_var <- sd(time_values,na.rm = T)

  fnormB_mean
  fnormB_var
  f1B_mean
  f1B_var
  fnormbeta_mean
  fnormbeta_var
  f1beta_mean
  f1beta_var
  time_mean
  time_var





