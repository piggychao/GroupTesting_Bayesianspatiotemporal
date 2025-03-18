############################################################
#
# Bayesian Spatio-temporal Regression Models for
#   Group Testing Data Under DT Protocol
#     With Unknown Se&Sp 
#
#   VAR Model
# 
############################################################

library(parallel)
library(MASS)
library(BayesLogit)
library(coda)
library(invgamma)
library(truncnorm)
library(spam)

rm(list = ls())
### Call Rcpp function for sampling y tilde
Rcpp::sourceCpp('../sampletildey_DT.cpp')

nm_cores <- 28

simu_func <- function(seed){
  set.seed(seed)
  ## Part I. Data Generation
  ## 1. Generate spatial-temporal random effects
  ### Create an adjacency matrix - W matrix
  ### Use the South Carolina State as an example 
  # de <- counties(state = "SC") #download a shapefile of counties directly into R
  # de.nb <- poly2nb(as(de, "Spatial")) #convert shapefile into neighbor/spatial object 
  # de.Wmatrix <-  nb2mat(de.nb, style = "B") #create an adjacency matrix - W matrix
  
  ### Alternative: import W matrix directly
  de.Wmatrix <-  readRDS("Wmatrix.rds") #create an adjacency matrix - W matrix
  
  S_r <- dim(de.Wmatrix)[1] #number of regions
  W_mtx <- matrix(as.numeric(de.Wmatrix), nrow = S_r)
  W_mtx_spar <- as.spam(W_mtx)
  D_mtx <- diag(apply(W_mtx, 1, sum)) #calculate D matrix
  D_mtx_spar <- as.spam(D_mtx)
  rho <- 0.9 #true rho in random effects
  tau_sqrt_tr <- 0.5 #true tau_sqrt in random effects
  sigma_r <- tau_sqrt_tr * solve(D_mtx_spar - rho*W_mtx_spar)
  
  Tps <- 10 #used in Scenario 1, Scenario 2
  # Tps <- 5 #used in Scenario 3
  zeta_true <- 0.99
  S_r <- dim(de.Wmatrix)[1]
  effects <- matrix(0, S_r, Tps+1)
  effects[,1] <- mvrnorm(1, mu = rep(0, S_r), Sigma = sigma_r) #xi_1
  for (i in 2:(Tps+1)){
    mu_effects <- zeta_true * effects[,i-1]
    effects[,i] <- mvrnorm(1, mu = mu_effects, Sigma = sigma_r)
  }
  effects <- effects - mean(effects) #impose sum-to-zero constraint to entire 1~T+1
  effects_new <- effects[,Tps+1]
  effects0 <- effects[,c(1:Tps)]
  
  ## 2. Generate individual data with random effects
  N <- 4000 #number of total observations
  #assign N individuals into S regions randomly 
  region_idx <- sample(1:S_r, size = N, replace = TRUE)
  #assign N individuals into T time points randomly 
  time_idx <- sample(1:Tps, size = N, replace = TRUE)
  
  Se_true <- c(0.95,0.98) #sensitivity for pool&individual
  Sp_true <- c(0.98,0.99) #specificity for pool&individual
  
  beta_true <- c(-4, 1, -1, 1)   #prevalence = 7-7.5% (used in Scenario 1, Scenario 3)
  # beta_true <- c(-4, 2, -0.5, 2) #prevalence = 15% (used in Scenario 2)
  p <- length(beta_true)
  sigma_x <- diag(p-2)
  xx <- mvrnorm(N, mu = rep(0, p-2), Sigma = sigma_x) 
  xx2 <- rbinom(N, 1, prob = 0.5)
  data_x <- cbind(rep(1,N), xx, xx2)
  data_e <- effects0[cbind(region_idx, time_idx)] #spatial-temporal random effects
  p_true <- exp(data_x %*% beta_true + data_e)/
    (1+exp(data_x %*% beta_true + data_e))
  # mean(p_true) #prevalence
  y_t <- rbinom(N, 1, p_true) #individual true status
  p_ind <- Se_true[2]*y_t + (1-Sp_true[2])*(1-y_t) 
  y_obs <- rbinom(N, 1, p_ind) #individual test results
  
  ## Prevalence in the T+1 time point
  x_pred <- c(1,0,0,0)
  p_pred <- exp(as.numeric(x_pred %*% beta_true) + effects_new)/
    (1+exp(as.numeric(x_pred %*% beta_true) + effects_new))
  
  ## Generate pool data - under DT protocol
  size <- 4 #pool size
  J <- N/size #number of pools per state
  
  z_t <- rep(NA, J)
  data_z <- rep(NA, J)
  for (j in 1:J){
    pool_index <- (j-1)*size
    pool_data <- y_t[(pool_index+1):(pool_index+size)]
    z_t[j] <- ifelse(sum(pool_data) > 0, 1, 0)
    p_pool <- Se_true[1]*(z_t[j]) + (1-Sp_true[1])*(1-z_t[j])
    data_z[j] <- rbinom(1,1,p_pool)
  }
  
  # Part II. MH-Gibbs Sampling
  
  ## Objective function in MH algorithm for rho in spatial effects
  objective <- function(vec_xi, tau_sqrt, num_tps, rho, alpha_rho, beta_rho, W_mtx, D_mtx){
    if (rho>0 & rho<0.9999) {
      Dist_mtx_r <- (1/tau_sqrt)*(D_mtx - rho*W_mtx)
      Dist_mtx_r_inv <- solve(Dist_mtx_r)
      p1 <- (-num_tps/2)*log(det(Dist_mtx_r_inv))
      p2 <- (-1/(2*tau_sqrt))*vec_xi
      p3 <- (alpha_rho-1) * log(rho) + (beta_rho-1) * log(1-rho)
      calcu <- p1 + p2 + p3
    } else {calcu <- -Inf}
    return(calcu)
  }
  
  ## Calculate xi_sum3 (with zeta_init) function
  fun_xi_sum3 <- function(e_input, zeta, rho, W_mtx, D_mtx){
    n_tps <- dim(e_input)[2]
    matrix_s0 <- D_mtx - rho*W_mtx #without tau_sqrt
    xi_sum3 <- 0
    result1 <- 0
    xi_sum3 <- as.numeric(t(e_input[,1]) %*% matrix_s0 %*% e_input[,1])
    for (t in 2:n_tps){
      xi_vec0 <- e_input[,t-1]
      xi_vec1 <- e_input[,t]
      result1 <- t(xi_vec1 - zeta*xi_vec0) %*% matrix_s0 %*% (xi_vec1 - zeta*xi_vec0)
      xi_sum3 <- xi_sum3 + result1
    }
    return(xi_sum3)
  }  
  
  MH_gibbs <- function(z, x, y_test, beta_init, var_beta, spatial_idx,
                       W_mtx, rho_init, alpha_r, beta_r,
                       tau_sqrt, tau_a, tau_b, zeta_init, sigma_sqr_zeta,
                       temporal_idx, num_tps, 
                       accuracy_init, Se1, Se2, Sp1, Sp2,
                       var_tune, iter_warmup, iter_sampling, bf_adapt
  ){
    S <- dim(D_mtx)[1]
    W_mtx <- as.spam(W_mtx)
    D_mtx <- diag(apply(W_mtx, 1, sum)) #calculate D matrix
    D_mtx <- as.spam(D_mtx)
    Dist_mtx_init <- solve(D_mtx - rho_init*W_mtx)
    region_time_idx <- cbind(spatial_idx, temporal_idx)
    x_pred <- c(1,0,0,0)
    
    effects_init <- matrix(0, S, num_tps+1)
    effects_init[,1] <- mvrnorm(1, mu=rep(0, S),
                                Sigma=tau_sqrt*Dist_mtx_init) #starting condition
    for (t in 2:(num_tps+1)){
      effects_init[,t] <- mvrnorm(1, mu=zeta_init*effects_init[,t-1],
                                  Sigma=tau_sqrt*Dist_mtx_init)
    }
    effects_init <- effects_init - mean(effects_init) #initialize xi
    effects_init0 <- effects_init[,c(1:num_tps)]
    effects_len <- length(as.vector(t(effects_init)))
    
    post_tau_a <- tau_a + S*(num_tps+1)/2 #posterior a in tau_r
    
    n <- dim(x)[1]
    p <- dim(x)[2]
    pool_num <- length(z)
    pool_size <- n/pool_num
    Se_p <- Se1 #initial pool_Se
    Sp_p <- Sp1 #initial pool_Sp
    Se_i <- Se2 #initial ind_Se
    Sp_i <- Sp2 #initial ind_Sp
    
    pool_ind_cpp <- rep(0:(pool_num-1), each=pool_size) #pool index for each individual
    pool_ind <- rep(1:pool_num, each=pool_size) #pool index for each individual
    ## Under DT, retest individuals in the positive pools
    retest_pool_idx <- which(z > 0) #find the index of positive pools
    retest_y <- y_test[pool_ind %in% retest_pool_idx]
    #initial z tilde
    z_til <- rep(0, pool_num)
    #retested individuals
    y_til_rt <- rep(0, length(retest_pool_idx))
    #initial y tilde
    p_init <- exp(x %*% beta_init + effects_init0[region_time_idx])/
      (1+exp(x %*% beta_init + effects_init0[region_time_idx]))
    y_init <- rbinom(n, 1, p_init)
    
    prior_var_inv <- (1/var_beta) * diag(p) #prior of beta's
    prior_mu <- rep(1, p)
    
    # beta's, tau, rho, random effects, zeta, Se&Sp*2
    gibbs_samples <- matrix(0, nrow = iter_warmup + iter_sampling, ncol = p+3+4+effects_len+S)
    
    for (g in 1:(iter_warmup + iter_sampling)){
      ## Sampling individual y_i
      p_est <- exp(x %*% beta_init + effects_init0[region_time_idx])/
        (1+exp(x %*% beta_init + effects_init0[region_time_idx]))
      y_init <- sampletildey_DT(pool_ind_cpp, p_est, z, y_init, y_test,
                              c(Se_p, Se_i), c(Sp_p, Sp_i))
      
      ## Sampling Sensitivity and Specificity - did not change
      ### Pool Se and Sp 
      z_matrix <- matrix(y_init, ncol = pool_size, byrow = TRUE)
      z_til <- ifelse(rowSums(z_matrix)>0, 1, 0)
      Se_p <- rbeta(1, sum(z * z_til) + accuracy_init,
                    sum((1-z) * z_til) + accuracy_init)
      Sp_p <- rbeta(1, sum((1-z) * (1-z_til)) + accuracy_init,
                    sum(z * (1-z_til)) + accuracy_init)
      ### Individual Se and Sp
      #find the y_tilde for re-test individuals
      y_til_rt <- y_init[pool_ind %in% retest_pool_idx]
      Se_i <- rbeta(1, sum(retest_y * y_til_rt) + accuracy_init,
                    sum((1-retest_y) * y_til_rt) + accuracy_init)
      Sp_i <- rbeta(1, sum((1-retest_y) * (1-y_til_rt)) + accuracy_init,
                    sum(retest_y * (1-y_til_rt)) + accuracy_init)
      
      ## Sampling Beta's
      b_pg <- as.numeric(rep(1,n))
      c_pg <- as.numeric(x %*% beta_init + effects_init0[region_time_idx])
      omega <- rpg(n, b_pg, c_pg)
      omega_mat <- as.spam(diag(omega))
      kappa <- y_init - b_pg/2
      h_ome <- kappa/omega - effects_init0[region_time_idx]
      var_ome_inv <- solve(prior_var_inv + t(x) %*% omega_mat %*% x)
      m_ome <- prior_var_inv %*% prior_mu + t(x) %*% omega_mat %*% h_ome
      post_mu <- var_ome_inv %*% m_ome
      beta_par <- mvrnorm(1, post_mu, Sigma = var_ome_inv) #posterior beta's
      beta_init <- beta_par
      
      ## Sampling S by T random effects 
      matrix_spar0 <- D_mtx - rho_init*W_mtx #without tau_sqrt
      matrix_spar <- (1/tau_sqrt)*matrix_spar0 
      for (t in 1:(num_tps+1)){
        time_t_idx <- which(temporal_idx == t)
        temp_idx <- spatial_idx[time_t_idx]
        indicator_mat <- as.spam(model.matrix(~ -1 + factor(temp_idx, levels = seq(1:S))))
        omega0 <- omega[time_t_idx]
        omega0_mat <- as.spam(diag(omega0))
        kappa0 <- kappa[time_t_idx]
        x0 <- x[time_t_idx,]
        if (t==1){
          econd <- effects_init[,t+1] #depend on xi_{t-1} and/or {t+1}
        } else if (t==(num_tps+1)){
          econd <- effects_init[,t-1]} else {
          econd <- effects_init[,t-1] + effects_init[,t+1]
          }
        if (t != (num_tps+1)){
          m_xi <- t(indicator_mat) %*% omega0_mat %*% indicator_mat +
            matrix_spar * (1+zeta_init^2)
          b_xi <- t(indicator_mat) %*% (kappa0 - omega0_mat%*%(x0%*%beta_init)) +
            zeta_init*(matrix_spar%*%econd)
          } else {
          m_xi <- matrix_spar
          b_xi <- zeta_init*(matrix_spar%*%econd)}
        m_xi_inv <- solve(m_xi)
        mu_xi <- m_xi_inv %*% b_xi
        effects_init[,t] <-  mvrnorm(1, mu=mu_xi, Sigma=m_xi_inv)
      }
      effects_init <- effects_init - mean(effects_init)
      effects_init0 <- effects_init[,c(1:num_tps)]
      # Predicted effects in the T+1 time point
      p_pred_mcmc <- exp(as.numeric(x_pred %*% beta_init) + effects_init[,num_tps+1])/
        (1+exp(as.numeric(x_pred %*% beta_init) + effects_init[,num_tps+1]))
      
      ## Sampling zeta
      xi_sum1 <- 0
      xi_sum2 <- 0
      for (t in 2:(num_tps+1)){
        xi_vec0 <- effects_init[,t-1]
        xi_vec1 <- effects_init[,t]
        result1 <- t(xi_vec0) %*% matrix_spar %*% xi_vec0
        result2 <- t(xi_vec0) %*% matrix_spar %*% xi_vec1
        xi_sum1 <- xi_sum1 + result1
        xi_sum2 <- xi_sum2 + result2
      }
      m_zeta_inv <- 1/(xi_sum1 + (1/sigma_sqr_zeta))
      zeta_init <- rtruncnorm(1, a=-1, b=1, mean=m_zeta_inv*xi_sum2, sd=sqrt(m_zeta_inv))
      
      ## Sampling tau_sqrt
      xi_sum3 <- fun_xi_sum3(effects_init, zeta_init, rho_init, W_mtx, D_mtx)
      post_tau_b <- tau_b + xi_sum3/2
      tau_sqrt <- rinvgamma(1, shape=post_tau_a, rate=post_tau_b) #posterior tau_sqrt
      
      ## Sampling rho by MH algorithm
      u <- runif(1, 0, 1)
      rho_s <- rtruncnorm(1, a=0, b=0.9999, mean=rho_init, sd=var_tune)
      truncnorm_prop <- dtruncnorm(rho_init, a=0, b=0.9999, mean=rho_s, sd=var_tune)/
        dtruncnorm(rho_s, a=0, b=0.9999, mean=rho_init, sd=var_tune)
      xi_sum3_s <- fun_xi_sum3(effects_init, zeta_init, rho_s, W_mtx, D_mtx)
      proportion <- exp(objective(xi_sum3_s, tau_sqrt,  num_tps+1, rho_s, alpha_r, beta_r, W_mtx, D_mtx) -
                          objective(xi_sum3, tau_sqrt,  num_tps+1, rho_init, alpha_r, beta_r, W_mtx, D_mtx)) *
        truncnorm_prop
      accept_prob <- ifelse(proportion <= 1, proportion, 1)
      if (u <= accept_prob){
        rho_init <- rho_s
      } else {
        rho_init <- rho_init
      }
      
      gibbs_samples[g, 1:p] <- beta_init
      gibbs_samples[g, p+1] <- tau_sqrt
      gibbs_samples[g, p+2] <- rho_init
      gibbs_samples[g, p+3] <- zeta_init
      gibbs_samples[g, (p+4):(p+3+4)] <- c(Se_p, Se_i, Sp_p, Sp_i)
      gibbs_samples[g, (p+4+4):(p+3+4+effects_len)] <- as.vector(t(effects_init)) #save all random effects
      gibbs_samples[g, (p+4+4+effects_len):(p+3+4+effects_len+S)] <- p_pred_mcmc
    }
    return(gibbs_samples[(iter_warmup:(iter_warmup+iter_sampling)),])
  }
  
  output <- MH_gibbs(data_z, data_x, y_obs, rep(0,p), 10, region_idx, 
                     W_mtx, 0.8, 2, 1,
                     0.2, 3, 1, 0.8, 0.8,
                     time_idx, 10, 
                     1, 0.8, 0.8, 0.8, 0.8, 
                     0.2, 2500, 5000, 50)
  
  mcmcobj <- mcmc(output)
  
  est_beta_mean <- summary(mcmcobj)$statistics[,"Mean"]
  est_beta_var <- apply(mcmcobj,2, var)
  est_beta_lci <- summary(mcmcobj)$quantiles[,1]
  est_beta_uci <- summary(mcmcobj)$quantiles[,5]
  true_par <- c(beta_true,tau_sqrt_tr, rho, zeta_true, Se_true, Sp_true, 
                as.vector(t(effects)), p_pred)
  convg <- pnorm(abs(geweke.diag(mcmcobj)$z), lower.tail = F)*2 
  results <- data.frame(true_pars = true_par,
                        est_mean = est_beta_mean,
                        est_var = est_beta_var,
                        est_lci = est_beta_lci,
                        est_uci = est_beta_uci,
                        convg = convg, 
                        seed = seed)
  return(results)
}

nm_set <- 100 #randomly choose a set number
n_sim <- 500
output <- mclapply(((nm_set-1)*n_sim+1):(nm_set*n_sim), 
                   simu_func, mc.cores = nm_cores)
time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
saveRDS(output, file = paste0('VAR_DT1_', time,'.rds'))