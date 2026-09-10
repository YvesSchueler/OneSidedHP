
################################################################################
# Code for On adjusting the one-sided Hodrick-Prescott filter (2020) 
# by Elias Wolf (FU-Berlin), Frieder Mokinski (Deutsche Bundesbank), and
# Yves Schüler (Deutsche Bundesbank)

# Version date 2023/05/04
# If you encounter any bug, please mail Yves Schüler at yves.schueler (at) bundesbank.de

# This program estimates the cyclical component of the adjusted one-sided HP
# filter (HP-1s*) by estimating the two-sided HP filter (HP-2s) on an expanding sample and keeping the last observation.

# For the exact references that are cited in this program, please see
# the paper.
#
#       Input:    
#                   y     - a Tx1 vector (y_1,...,y_T)', where T is the
#                   number of observations. 
#                      
#                   lm_2  - a scalar. This is the value of the smoothing
#                   parameter usually employed for the two-sided HP filter
#                   (HP-2s). This input is optional. The default is lm_2=1,600.
#
#                   opt    - 0/1. If opt=1, the adjustment parameters are
#                   obtained through the optimization routine (computationally more intensive). This input is optional. The default is
#                   opt=0 (computationally less intensive).
#
#                   sample - 0/1. If sample = 1, the adjustment takes sample size into account. May be used for large smoothing
#                   parameters in combination with small sample sizes, for instance, if
#                   lm_2 = 400,000 and T < 100. 
#
#      Output:  
#                   ycycle_adj  - a Tx1 vector with the adjusted extracted cyclical component
#                  
#                   lm_1        - a scalar. This is the adjusted value of the smoothing parameter used as an input to the one-sided HP filter
#
#                   kappa       - a scalar. This is the mulitplicative scaling factor, with wich the cyclical component of the
#                   one-sided HP filter is rescaled
#
#       Example:
#                   ycycle_adj, lm_1, kappa = adj1s_hpfilter(y,1600)
#                   yields a Tx1 vector of the extracted cyclical component
#                   using a smoothing parameter of 650 and scaling
#                   parameter of size 1.1513 for the one-sided HP filter
################################################################################                  


adj1s_hpfilter <- function(y, lm_2, opt=TRUE, sample=FALSE){
  
  if (is.numeric(y) == FALSE){
    y = as.numeric(y)
  }
  
  if (length(y) == 0){
    stop("Time series is missing or empty")
  }

  if (length(dim(y)) > 1){
    stop("y contains more than one time series: the program is designed for one series only")
  }
      
  if (length(y) < 3){
    stop("y contains less than 3 observations: please check your data input")
  }
  
  if (any(is.na(y)) == TRUE){
    stop("y contains missing values")
  }
      
  poly_lm <- function(lm_2){
    return(-5.071743 +  0.4080780*lm_2  -0.749537*lm_2^(1/2) + 14.50001*lm_2^(1/3) -50.36123*lm_2^(1/4) + 41.44369*lm_2^(1/5))
  } 
  
  poly_kappa <- function(lm_2){
    return(0.984928 +  0.547209*(1/lm_2)  - 0.703636*(1/sqrt(lm_2)) + 3.291657*(1/(lm_2)^(1/3)) - 3.165415*(1/(lm_2)^(1/4)) + 1.759527 *(1/(lm_2)^(1/5)))
  }
  
  if (opt == FALSE){
    lm_1 <- poly_lm(lm_2)
    kappa <- poly_kappa(lm_2)
  }
    
  if (opt == TRUE){
  # finds adjustment parameters lm_1 and kappa by harmonizing PTF of HP-2s
  # and HP-1s
    
    if (sample == TRUE){
      T_eff <- length(y)  # size of vector of filter weight for observation t=T
    } else {
      T_eff <- 1000
    }
    
    ptf_cost <- function(params, lm_2, T_eff){
      
      # Parameters to be optimized
      lm_1 <- params[1]
      ratio <- params[2]
      
      I <- diag(T_eff)
      Q_t <- diff(I, differences=2)
      A_inv <- solve(I + lm_1*t(Q_t)%*%Q_t)
      
      wts_one <- A_inv[T_eff,]
      
      # FFT of the weights to obtain PTF of the one-side HP-filter
      grid <- 2*pi*seq(0,floor(T_eff/2), by=1)/T_eff
      
      k <- seq(-(T_eff-1), 0, by=1)
      hpgain <- numeric(length = length(grid))
      
      for (i in 1:length(grid)){
        hpgain[i] <- 1 - t(rep(1, length(k)))%*%(wts_one*exp(1i*grid[i]*k))
      }

      hp1_ptf <- ratio^2*abs(hpgain)^2
      
      # PTF of the two-sided HP_filter 
      hp2_ptf <- ((4*lm_2*(1-cos(grid))^2/(1+4*lm_2*(1-cos(grid))^2))^2)
      
      # Calculate loss function defined as the squared distance    
      cost <- sum((hp1_ptf-hp2_ptf)^2)
      
      return(cost)
    }

    # Initial values and optimization of the cost function
    init_vals <- c(poly_lm(lm_2), poly_kappa(lm_2))
    res <- optim(par=init_vals, fn=ptf_cost, lm_2 = lm_2, T_eff = T_eff, method='BFGS')
        
    # Collect optimal parameters
    lm_1 <- res$par[1]
    kappa <- res$par[2]
  }
        
  hp_one <- function(lm_1, kappa, y){
    T_series <- length(y)
    cycle_os <- numeric(length = T_series)
    
    for (i in 3:T_series){
      
      I <- diag(i)
      Q_t <- diff(I , differences = 2)
      A_inv <- solve(I + lm_1*t(Q_t)%*%Q_t)

      psi_t <- y[i] - A_inv[i,]%*%(y[1:i])

      cycle_os[i] <- psi_t

    }
    return(kappa*cycle_os)
  }
  
  # Apply the one-sided HP-Filter to the series
  ycycle_adj = hp_one(lm_1, kappa, y)
    
  f_out <- list("ycycle_adj" = ycycle_adj[3:length(ycycle_adj)],
                "lm_1" = lm_1,
                "kappa" = kappa)
  return(f_out)
}

