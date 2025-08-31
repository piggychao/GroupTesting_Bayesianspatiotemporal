############################################################
#
# Bayesian Spatio-temporal Regression Models for
#   Group Testing Data Under IT Protocol
#      With Known Se&Sp 
#
#        CAR-AR(1) Model
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
Rcpp::sourceCpp('../sampletildey_IT.cpp')

nm_cores <- 28

simu_func <- function(seed){
  set.seed(seed)

  ## Part I. Data Generation
  ## 1. Generate spatial random effects
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
  spatial_rand <-  mvrnorm(1, mu = rep(0, S_r), Sigma = sigma_r)
  spatial_rand <- spatial_rand - mean(spatial_rand)

  ## 2. Generate temporal random effects 
  Tps <- 10 #used in Scenario 1, Scenario2
  # Tps <- 5 #used in Scenario 3
  sigma_r_t <- 0.25
  phi_t <- 0.8
  
  ## Generate AR1 structure
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  tempo_rand <- mvrnorm(1, mu = rep(0, Tps), 
                        Sigma = sigma_r_t/(1-phi_t^2)* ar1_cor(Tps,phi_t))
  tempo_rand <- tempo_rand - mean(tempo_rand)
  tempo_rand_new <- rnorm(1, mean = phi*tempo_rand[Tps], sd=sqrt(sigma_r_t))

  ## 2. Generate individual data with spatial random effects
  N <- 4000 #number of total observations
  #assign N individuals into S regions randomly 
  region_idx <- sample(1:S_r, size = N, replace = TRUE)
  #assign N individuals into T time points randomly 
  time_idx <- sample(1:Tps, size = N, replace = TRUE)
  
  Se <- 0.98 #known sensitivity for individual testing
  Sp <- 0.99 #known specificity for individual testing

  ## Generate individual data
  beta_true <- c(-4, 1, -1, 1)#prevalence = 7-7.5% (used in Scenario 1, Scenario 3)
  # beta_true <- c(-4, 2, -0.5, 2) #prevalence = 15% (used in Scenario 2)
  p <- length(beta_true)
  sigma_x <- diag(p-2)
  xx <- mvrnorm(N, mu = rep(0, p-2), Sigma = sigma_x) 
  xx2 <- rbinom(N, 1, prob = 0.5)
  data_x <- cbind(rep(1,N), xx, xx2)
  data_s <- spatial_rand[region_idx] #xi
  data_t <- tempo_rand[time_idx] #theta
  p_true <- exp(data_x %*% beta_true + data_s + data_t)/
    (1+exp(data_x %*% beta_true + data_s + data_t))
  #mean(p_true) #prevalence
  y_t <- rbinom(N, 1, p_true) #individual true status
  p_ind <- Se*y_t + (1-Sp)*(1-y_t) 
  y_obs <- rbinom(N, 1, p_ind) #individual test results

  ## Prevalence in next time point
  x_pred <- c(1,0,0,0)
  p_pred <- exp(as.vector(x_pred %*% beta_true) + spatial_rand + tempo_rand_new)/
    (1+exp(as.vector(x_pred %*% beta_true) + spatial_rand + tempo_rand_new))

  # Part II. MH-Gibbs Sampling
  ## Objective function in MH algorithm for rho in spatial effects
  objective <- function(xi, tau_sqrt, rho, alpha_rho, beta_rho, W_mtx, D_mtx){
    if (rho>0 & rho<0.99) {
      Dist_mtx_r <- D_mtx - rho*W_mtx
      Dist_mtx_r_inv <- solve(Dist_mtx_r)
      p1 <- (-1/2)*log(det(Dist_mtx_r_inv))
      p2 <- as.numeric((-1/(2*tau_sqrt)) * (t(xi) %*% Dist_mtx_r %*% xi))
      p3 <- (alpha_rho-1) * log(rho) + (beta_rho-1) * log(1-rho)
      calcu <- p1 + p2 + p3
    } else {calcu <- -Inf}
    return(calcu)
  }

  MH_gibbs <- function(x, y, beta_init, var_beta, spatial_idx,
                       W_mtx, rho_init, alpha_r, beta_r,
                       tau_sqrt, tau_a, tau_b, 
                       temporal_idx, num_tps, 
                       sigma_theta_init, a_temp, b_temp, phi_temp,
                       Se, Sp,
                       var_tune, iter_warmup, iter_sampling, bf_adapt
  ){
    D_mtx <- diag(apply(W_mtx, 1, sum)) #calculate D matrix
    Dist_mtx_init <- solve(D_mtx - rho_init*W_mtx)
    S <- dim(D_mtx)[1]
    G_mtx <- as.matrix(model.matrix(~ as.factor(spatial_idx) - 1))
    xi_init <-  mvrnorm(1, mu = rep(0, S), 
                        Sigma = tau_sqrt * Dist_mtx_init) #initialize xi
    post_tau_a <- tau_a + S/2 #posterior a in tau_r
    H_mtx <- as.matrix(model.matrix(~ as.factor(temporal_idx) - 1))
    Tps <- num_tps
    #variance-covariance matrix of temporal effects
    AR_mtx <- ar1_cor(Tps, phi_temp)
    sigma_t <- sigma_theta_init / (1-phi_temp^2) * AR_mtx
    sigma_t_inv <- solve(sigma_t)
    theta_init <- mvrnorm(1, mu = rep(0, Tps), Sigma = sigma_t)
    n <- length(y)
    p <- dim(x)[2]
    y_init <- y
    Se1 <- Se #Sensitivity for individual testing
    Sp1 <- Sp #Specificity for individual testing
    prior_var_inv <- (1/var_beta) * diag(p)
    prior_mu <- rep(1, p)
    post_sigma_a <- a_temp + Tps/2 #posterior a in sigma_theta
    # beta's, mean of xi, tau, rho
    gibbs_samples <- matrix(0, nrow = iter_warmup + iter_sampling, ncol = p+2+S+Tps+2+S)
    
    for (g in 1:(iter_warmup + iter_sampling)){
      #sampling individual y_i
      p_est <- exp(x%*%beta_init + xi_init[spatial_idx] + theta_init[temporal_idx])/
        (1 + exp(x%*%beta_init + xi_init[spatial_idx] + theta_init[temporal_idx]))
      ytilde <- sampletildey_IT(p_est, y_obs, Se1, Sp1)
      y_init <- ytilde
      
      b_pg <- as.numeric(rep(1,n))
      c_pg <- as.numeric(x %*% beta_init + 
                           xi_init[spatial_idx] + theta_init[temporal_idx])
      omega <- rpg(n, b_pg, c_pg)
      omega_mat <- as.spam(diag(omega))
      kappa <- y_init - b_pg/2
      h_ome <- kappa/omega - xi_init[spatial_idx] - theta_init[temporal_idx]
      var_ome_inv <- solve(prior_var_inv + t(x) %*% omega_mat %*% x)
      m_ome <- prior_var_inv %*% prior_mu + t(x) %*% omega_mat %*% h_ome
      post_mu <- var_ome_inv %*% m_ome
      beta_par <- mvrnorm(1, post_mu, Sigma = var_ome_inv) #posterior beta's
      xi_post_m_inv <- solve(t(G_mtx) %*% omega_mat %*% G_mtx + (1/tau_sqrt)*
                               (D_mtx - rho_init*W_mtx))
      xi_post_b <- t(G_mtx) %*% (kappa - omega_mat %*% 
                                   (x %*% beta_par + H_mtx %*% theta_init))
      xi_post_mu <- xi_post_m_inv %*% xi_post_b
      xi_init <- mvrnorm(1, xi_post_mu, Sigma = xi_post_m_inv) #posterior xi
      xi_init <- xi_init - mean(xi_init)
      theta_post_m_inv <- solve(t(H_mtx) %*% omega_mat %*% H_mtx + sigma_t_inv)
      theta_post_b <- t(H_mtx) %*% (kappa - omega_mat %*% 
                                      (x %*% beta_par + G_mtx %*% xi_init))
      theta_post_mu <- theta_post_m_inv %*% theta_post_b
      # posterior theta
      theta_init <- mvrnorm(1, theta_post_mu, Sigma = theta_post_m_inv)
      theta_init <- theta_init - mean(theta_init)
      ## MH algorithm for rho
      u <- runif(1, 0, 1) 
      rho_s <- rtruncnorm(1, a=0, b=0.99, mean=rho_init, sd=var_tune) #revise: truncatednorm
      truncnorm_prop <- dtruncnorm(rho_init, a=0, b=0.99, mean=rho_s, sd=var_tune)/
        dtruncnorm(rho_s, a=0, b=0.99, mean=rho_init, sd=var_tune)
      proportion <- exp(objective(xi_init, tau_sqrt, rho_s, alpha_r, beta_r,
                                  W_mtx, D_mtx) - 
                          objective(xi_init, tau_sqrt, rho_init, alpha_r, beta_r,
                                    W_mtx, D_mtx)) * truncnorm_prop 
      accept_prob <- ifelse(proportion <= 1, proportion, 1)
      if (u <= accept_prob){
        rho_init <- rho_s
      } else {
        rho_init <- rho_init
      }
  
      post_tau_b <- tau_b + (t(xi_init) %*% (D_mtx - rho_init*W_mtx) %*% xi_init)/2 
      tau_sqrt <- rinvgamma(1, shape=post_tau_a, rate=post_tau_b) #posterior tau_sqrt
      post_sigma_b <- b_temp + (1/2) * (1-phi_temp^2) * (t(theta_init) %*% solve(AR_mtx) %*% theta_init)
      var_temp <- rinvgamma(1, shape=post_sigma_a, rate=post_sigma_b) #posterior sigma_theta
      sigma_t <- var_temp / (1-phi_temp^2) * AR_mtx
      sigma_t_inv <- solve(sigma_t)
      gibbs_samples[g, 1:p] <- beta_par
      beta_init <- beta_par
      gibbs_samples[g, p+1] <- tau_sqrt
      gibbs_samples[g, p+2] <- rho_init
      gibbs_samples[g, (p+3):(p+2+S)] <- xi_init
      gibbs_samples[g, (p+2+S+1):(p+2+S+Tps)] <- theta_init
      gibbs_samples[g, p+2+S+Tps+1] <- var_temp
      
      ## Prediction in next time point
      theta_pred <- rnorm(1, phi_temp*theta_init[Tps], sqrt(var_temp))
      gibbs_samples[g, p+2+S+Tps+2] <- theta_pred
      p_pred2 <- exp(as.vector(x_pred %*% beta_par) + xi_init + theta_pred)/
        (1+exp(as.vector(x_pred %*% beta_par) + xi_init + theta_pred))
      gibbs_samples[g, (p+2+S+Tps+3):(p+2+S+Tps+2+S)] <- p_pred2
    }
    return(gibbs_samples[(iter_warmup:(iter_warmup+iter_sampling)),])
  }

  output <- MH_gibbs(data_x, y_obs, rep(0,p), 10, region_idx, 
                    W_mtx, 0.2, 2, 1,
                    0.2, 3, 3,
                    time_idx, 10,
                    0.1, 5, 2, phi_t, 0.98, 0.99,
                    0.8, 2500, 5000, 50)

  mcmcobj <- mcmc(output)
  
  est_beta_mean <- summary(mcmcobj)$statistics[,"Mean"]
  est_beta_var <- apply(mcmcobj,2, var)
  est_beta_lci <- summary(mcmcobj)$quantiles[,1]
  est_beta_uci <- summary(mcmcobj)$quantiles[,5]
  true_par <- c(beta_true,tau_sqrt_tr, rho, spatial_rand, tempo_rand, sigma_r_t, 
                 tempo_rand_new, p_pred)
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

nm_set <- 3 #randomly choose a set number
n_sim <- 500
output <- mclapply(((nm_set-1)*n_sim+1):(nm_set*n_sim), 
                   simu_func, mc.cores = nm_cores)
time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
saveRDS(output, file = paste0('CARAR_IT_', time,'.rds'))