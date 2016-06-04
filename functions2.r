library(MASS)
library(doParallel)
library(foreach)
library(FNN)
library(caret)

# Definitions
i1 <- complex(imaginary = 1)

# Local polynomial regression 
lp <- function(YY, XX, xx, b) {
    K <- diag(dnorm(c(xx)/b))
    B <- solve(t(XX) %*% K %*% XX) %*% (t(XX) %*% K %*% YY)
    return(B)
}

lp_error <- function(YY, XX, betahat, xx, bw) {
    M <- YY - XX %*% betahat
    K <- dnorm(xx, sd = bw)
    return(mean(M^2*K))
}



# Gamma matrices
g1 <- function(x) matrix(c(1,1,x[1],x[2]), 2, 2)
g1inv <- function(x)  solve(g1(x))
g2 <- function(x) matrix(c(0, 1, 0, x[2]), 2, 2)
g3 <- function(x)  g1inv(x) %*% g2(x)



# Shock moments and fixed coefficients beta_{t}
estimate_shocks_and_fixed <- function(Y1, Y2, X1, X2, Z1, Z2, shocks_bw) {

    # First moments
    DY <- Y2 - Y1
    coeffs <- lp(YY = DY, 
                XX = cbind(1, X2, (-1)*Z1, Z2), 
                xx = X2 - X1, 
                b = shocks_bw)
    EU <- coeffs[1]
    EV <- coeffs[2]
    muW <- matrix(c(EU, EV), 2, 1)
    
    if (is.null(Z1) || is.null(Z2)) {
        b_Z <- NULL
        Zb <- 0
    } else {
        b_Z = matrix(coeffs[3:length(coeffs)], ncol = 1)
        Zb <- cbind( (-1)*Z1, Z2) %*% b_Z  
    }
    # Second moments
    shock_cov <- lp((DY - EU - EV*X2 - Zb)^2, 
                    cbind(1, 2*X2, X2^2), 
                    X2 - X1,
                    b = shocks_bw)

    return(list(shock_means = muW, 
                shock_variances = shock_cov[c(1,3)], 
                shock_covariance = shock_cov[2],
                fixed_coeffs = b_Z))
}

    
normal_cf <- function(muAB, chol_sigmaAB, tt) {
    tt <- matrix(tt, 1, 2)
    S <-  matrix(c(chol_sigmaAB[1], 0, chol_sigmaAB[2], chol_sigmaAB[3]),2,2)
    sigmaAB <- t(S) %*% S
    return(exp(1i * tt %*% muAB - .5 * tt %*% sigmaAB %*% t(tt)))
}


# Integrates the ratio over t  
compute_rhs <- function(Y1, Y2, X1, X2, Z1, Z2, det_bnds, modulus_bnds, reg_params, tts) { 

    # Activate parallel computing if available
    registerDoParallel(cores = detectCores())

    shocks_bw <- estimate_shocks_and_fixed_optimal_bw(Y1, Y2, X1, X2, Z1, Z2) 
    shocks_and_fixed <- estimate_shocks_and_fixed(Y1, Y2, X1, X2, Z1, Z2, shocks_bw)

    # Estimate for each (t1, t2) pair
    rhs <- foreach(i = 1:nrow(tts), .combine = rbind) %dopar% {
        tt = matrix(as.numeric(tts[i,]), 1, 2)
        mean_ratio(Y1, Y2, X1, X2, Z1, Z2, tt, shocks_and_fixed, det_bnds, modulus_bnds, reg_params)
    }
    return(rhs)
}


# Create Fourier variable grid
create_tts <- function(n_tt) {    
    t_grid <- runif(n = n_tt, min = -.5, max = .5)
    return(expand.grid(t_grid, t_grid))
}


loss <- function(rhs, tts, p) {
    
    lhs <- foreach (i = 1:nrow(tts), .combine = c) %dopar% {
        tt = matrix(as.numeric(tts[i,]), 1, 2)
        normal_cf(p[1:2], p[3:5], tt)
    }
    
    # Output the squared distance
    difference <- rhs - lhs
    return(c(Re(difference %*% Conj(difference))))
}



# Fit on training set
fit_estimates <- function(rhs, p_init, tts) {
    estimates <- foreach (col = rhs, .combine = cbind) %dopar% {
        result <- nlm(function(p) loss(as.complex(col), tts, p), p_init)$estimate
        }
    return(estimates)
}
 

# Score on test set
score_estimates <- function(rhs, estimates, tts) {
    test_errors <-  foreach (i = seq(ncol(rhs)), .combine = c) %dopar% {
        error <- loss(rhs[,i], tts, estimates[,i])
    }
    return(test_errors)
}





build_num_var <- function(Y1, Y2, X1, X2, Z1, Z2, b_Z, tt) {
    
    n_z <- length(b_Z)/2
    n_obs <- length(Y1)
    tt <- matrix(tt, 1, 2)
    
    # Subtract terms involving fixed coeffs (i.e., Z_{t}'b_{t})
    if (is.null(Z1) || is.null(Z2)) {
        Y1tilde <- Y1
        Y2tilde <- Y2
    } else {
        b_Z1 <- matrix(b_Z[1:n_z], n_z, 1)
        b_Z2 <- matrix(b_Z[(n_z+1):(2*n_z)], n_z, 1)
        Y1tilde <- Y1 - Z1 %*% b_Z1
        Y2tilde <- Y2 - Z2 %*% b_Z2
    }
    
    # Build dependent variable
    costY <- matrix(0, n_obs, 1)
    sintY <- matrix(0, n_obs, 1)
    for (i in seq(n_obs)){
        tY <- tt %*% g1inv(c(X1[i], X2[i])) %*% matrix(c(Y1tilde[i], Y2tilde[i]),2,1)
        costY[i] <- cos(tY)
        sintY[i] <- sin(tY)
    }
    return(list(costY = costY, sintY = sintY))
}
 

# Numerator
create_numerator <- function(num_var, X1, X2, Z1, Z2, det_bnds, reg_params) {
    
    # This will be a matrix 
    real_part <- fit_regression(num_var$costY, X1, X2, Z1, Z2, reg_params)
    imag_part <- fit_regression(num_var$sintY, X1, X2, Z1, Z2, reg_params)
    numerator <- real_part + i1*imag_part
    
    # Make nan of each numerator that doesn't satisfy det_bnd
    numerators <- foreach(num = iter(numerator, by = "col"), .combine = cbind) %dopar% {    
        foreach(det_bnd = det_bnds, .combine = cbind) %dopar% {
            discard_low_det_obs(num, X1, X2, det_bnd)
        }
    }
    return(numerators)
}



# Denominator variable
create_denominator <- function(shock_means, shock_variances, shock_covariance,
                               X1, X2, Z1, Z2, tt, modulus_bnds) {
    
    # Make sure dimensions are correct
    n_obs <- dim(X1)[1]
    muW <- matrix(shock_means, 2, 1)
    sigmaW <- matrix(c(shock_variances[1], shock_covariance,
                       shock_covariance, shock_variances[2]), 2, 2)

    # Compute denominator
    denom <- matrix(0, n_obs, 1) 
    for (i in seq(n_obs)) {
        denom[i] <- exp(i1* tt %*% g3(c(X1[i], X2[i])) %*% muW -
                        .5 * tt %*% g3(c(X1[i], X2[i])) %*% sigmaW %*%
                                        t(g3(c(X1[i], X2[i]))) %*% t(tt))
        }

    # Increase the modulus for each bound
    denoms <- foreach(modulus_bnd = modulus_bnds, .combine = cbind) %dopar% { 
        increase_modulus(denom, modulus_bnd)
    }
    return(denoms)
}


# Right-hand size variable for given Fourier variable t
mean_ratio <- function(Y1, Y2, X1, X2, Z1, Z2, tt, 
                       shocks_and_fixed, det_bnds, modulus_bnds, reg_params) {
    # Numerator variable
    num_var <- build_num_var(Y1, Y2, X1, X2, Z1, Z2, 
                             shocks_and_fixed$fixed_coeffs, tt)

    # Creates several numerators and denominators according to tuning parameters
    numerators <- create_numerator(num_var, X1, X2, Z1, Z2, det_bnds, reg_params) 
    denominators <- create_denominator(shocks_and_fixed$shock_means, 
                                     shocks_and_fixed$shock_variances, 
                                     shocks_and_fixed$shock_covariance,
                                      X1, X2, Z1, Z2, tt, modulus_bnds)  
   
    # One ratio for each (reg_param, det_bnd, modulus_bnd) combination 
    ratios <- foreach(n = numerators, .combine = cbind) %:%
                foreach(d = denominators, .combine = cbind) %dopar% {
                    out <- n/d
              }
    
    # Mean ratios
    mean_ratios <- foreach(r = ratios, .combine = c) %dopar% { 
        mean(r[is.finite(r)])
    }     
    return(mean_ratios)
}

# Increases |z| keeping arg(z) constant
increase_modulus <- function(z, modulus_bnd) {
    r <- abs(z)
    z[r < modulus_bnd] <- z[r < modulus_bnd]/r[r < modulus_bnd]*modulus_bnd
    return(z)
}

# Indexes observations with lower Gamma1 determinant
discard_low_det_obs <- function(v, X1, X2, det_bnd) {
    idx <- apply(cbind(X1,X2), 1, function(x) det(g1(x)) < det_bnd)
    v[idx] <- NaN
    return(v)
}


# Auxiliary function that generates data satisfying model assumptions
create_data <- function(n_obs, with_Z = TRUE) {

    # Draw coefficients for first period
    mAB <- c(1,2)
    covAB <- matrix(c(2,1,1,2),nrow = 2)
    AB <- mvrnorm(n = n_obs, mu = mAB, Sigma = covAB)
    A <- AB[,1,drop = F]; B <- AB[,2,drop = F]

    # Draw shocks in second period
    mUV <- c(.5,.5)
    covUV <- matrix(c(1,0,0,1), nrow = 2)
    W <- mvrnorm(n = n_obs, mu = mUV, Sigma = covUV)
    U <- W[,1]; V <- W[,2]
   
    # Draw regressors
    # X1 correlated with A,B
    X1 <-  .2*A + .5*B + 
        .5*A^2 - .2*B^2 + 
rnorm(n = 1,mean = 0, sd = sqrt(5))
    
    # De-mean X1
    X1 <- X1 - mean(X1) 
    X_sd <- sd(X1)
    
    # X2 is not correlated with A,B, but has same mean (zero) and std. dev.
    X2 <- matrix(rnorm(n = n_obs, mean = 0, sd = X_sd), n_obs, 1)
    X <- cbind(X1, X2)


    # Draw Z1 and Z2
    if (with_Z) {
        Z11 = rnorm(n  = n_obs, mean = 3, sd = 2)
        Z21 = rnorm(n  = n_obs, mean = 4, sd = 2)
        Z12 = rnorm(n  = n_obs, mean = 3, sd = 2)
        Z22 = rnorm(n  = n_obs, mean = 4, sd = 2)
    } else {
        Z11 = 0
        Z21 = 0
        Z12 = 0
        Z22 = 0
    }

    # Compute dependent variable
    Y1 <- A + B*X1 + 1*Z11 + 2*Z21
    Y2 <- (A + U) + (B + V)*X2 + 1*Z12 + 2*Z22
    Y <- cbind(Y1, Y2)

    # Bind everything in a list to return
    dataset <- list()
    dataset$X1 <- X1
    dataset$X2 <- X2
    dataset$Y1 <- Y1
    dataset$Y2 <- Y2
    dataset$A <- A
    dataset$B <- B
    dataset$W <- W
    dataset$Z1 <- cbind(Z11, Z21)
    dataset$Z2 <- cbind(Z12, Z22)
    dataset$beta_1 <- c(1, 2)
    dataset$beta_2 <- c(1, 2)


    return(dataset)
}




estimate_shocks_and_fixed_optimal_bw <- function(Y1, Y2, X1, X2, Z1, Z2, 
                            bw_lower = .05, bw_upper = .5, n_bws = 20, n_folds = 10) {
    
    n_obs <- dim(Y1)[1]
    shocks_bws <- seq(bw_lower, bw_upper, length.out = n_bws)
    folds = createFolds(seq(n_obs), k = n_folds) 
    
    cv_mse <- foreach(shocks_bw = shocks_bws, .combine = cbind) %dopar% {
        bw_mse <- foreach(fold = folds, .combine = c) %dopar% {
            shocks_and_fixed <- estimate_shocks_and_fixed(Y1 = Y1[-fold,,drop=F], 
                                                         Y2 = Y2[-fold,,drop=F], 
                                                         X1 = X1[-fold,,drop=F], 
                                                         X2 = X2[-fold,,drop=F], 
                                                         Z1 = Z1[-fold,,drop=F], 
                                                         Z2 = Z2[-fold,,drop=F], 
                                                         shocks_bw = shocks_bw)
            errors <- rep(0, n_folds)
            for (i in fold) {
                errors[i] <-  lp_error(YY = Y2[i] - Y1[i],
                                    XX = cbind(1, X2[-i,], (-1)*Z1[-i,], Z2[-i,]),
                                    xx = X2[i] - X1[i],
                                    betahat = c(shocks_and_fixed$shock_means,
                                                shocks_and_fixed$fixed_coeffs),
                                    bw = shocks_bw)
            }
            mean(errors, na.rm = TRUE)
        }    
    }
    
    msem <- apply(cv_mse, 2, function(x) mean(x, na.rm = TRUE))
    msesd <- apply(cv_mse, 2, function(x) sd(x, na.rm = TRUE))
    return(shocks_bws[which.min(as.numeric(msem + msesd))])
}


fit_regression <- function(Y, X1, X2, Z1 = NULL, Z2 = NULL, reg_params) {
        estimates <- foreach(k = reg_params, .combine = cbind) %dopar% {
        model <- tryCatch( knn.reg(train = cbind(X1, X2, Z1, Z2), y = Y, k  = k), 
                          error = function(e) { 
                              warning(sprintf("Skipped KNN for k = %f", k))
                              return(NULL)
                          })
        model$pred
    }
    return(as.matrix(estimates))
}


fhhps <- function(Y1, Y2, X1, X2, Z1, Z2,
                  mu_init = c(0,0), 
                  sigma_init = matrix(c(1,0,0,1), 2, 2),
                  det_bnds = seq(.05, .1, .01), 
                  modulus_bnds = seq(.05, .1, .01),
                  reg_params = c(1,3,4,5),
                  n_tt = 10)  {
    
    p_init <- c(mu_init, chol(sigma_init)[c(1,3,4)])
    
    tts_train <- create_tts(n_tt)
    rhs_train <- compute_rhs(Y1, Y2, X1, X2, Z1, Z2, det_bnds, modulus_bnds, reg_params, tts_train)
    estimates <- fit_estimates(rhs_train, p_init, tts_train)
 
    tts_test <- create_tts(n_tt)
    rhs_test <- compute_rhs(Y1, Y2, X1, X2, Z1, Z2, det_bnds, modulus_bnds, reg_params, tts_test)
    scores <- score_estimates(rhs_test, estimates, tts_test)
    table <- build_table(modulus_bnds, det_bnds, reg_params, n_tt, scores, estimates)

    return(table)
    
}



build_table <- function(modulus_bnds, det_bnds, reg_params, n_tt, scores, estimates) {
    
    parameters <- expand.grid(modulus_bnds, det_bnds, reg_params)
    means <- t(estimates[1:2,])
    covariances <- foreach (i = seq(ncol(estimates)), .combine = rbind) %do% {
        L <- matrix(c(estimates[3,i],estimates[4,i],0,estimates[5,i]), 2, 2)
        Sigma <- t(L) %*% L
        c(Sigma[1,1], Sigma[1,2], Sigma[2,2])
    }
    n_tts <- rep(n_tt, dim(covariances)[1])
    table <- cbind(parameters, scores, n_tts, means, covariances)
    colnames(table) <- c("modulus_bnd", "det_bnd", "reg_params",  "score", "n_tt", 
                         "EA", "EB", "VA", "CAB", "VB")
    return(table)
}

