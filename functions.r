library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(plyr)
library(caret)
library(doParallel)
library(foreach)
library(mgcv)
library(gam)

#Get C++ functions (takes about a few moments to compile)
#sourceCpp("cpp_functions.cpp")

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



estimate_conditional_first_moments <- function(Y1, Y2, X1, X2, Z1, Z2, lb) {

    n_obs <- dim(Y1)[1]
    XZ <- cbind(X1, X2, Z1, Z2)
    n_regs <- dim(XZ)[2]

    # Local regression
    EY1 <- fit_model(Y1, X1, X2, Z1, Z2, at_xx = FALSE)
    EY2 <- fit_model(Y2, X1, X2, Z1, Z2, at_xx = FALSE)
    EDY2 <- fit_model(Y2 - Y1, X1, X2, Z1, Z2, at_xx = TRUE)
               
    # Invert matrices to get results
    cmoments <- matrix(0, 2, n_obs)
    for (i in seq(n_obs)) {
        b <- matrix(c(EY1[i], EY2[i] - EDY2[i]), 2, 1)
        A <- matrix(c(1, 1, X1[i], X2[i]), 2, 2) 
    
        if (rcond(A) > lb) { 
	    cmoments[,i] <- Re(solve(A) %*% b)	
        } else {
	    cmoments[,i] <- NaN
        }
    }
    return(cmoments)
}


estimate_conditional_second_moments <- 
    function(Y1, Y2, X1, X2, Z1, Z2, EA1_x, EB1_x, lb) {
    
    n_obs <- dim(Y1)[1]
    cmoments <- matrix(0, 3, n_obs)
    
    # Local Regressions
    EY1 <- fit_model(Y1, X1, X2, Z1, Z2)
    EDY2 <- fit_model(Y2 - Y1, X1, X2, Z1, Z2, at_xx = TRUE) 
    EY1sq <- fit_model(Y1^2, X1, X2, Z1, Z2) 
    EY2sq <- fit_model(Y2^2, X1, X2, Z1, Z2)
    EY1Y2 <- fit_model(Y1*Y2, X1, X2, Z1, Z2)
    EDY2sq <- fit_model((Y2 - Y1)^2, X1, X2, Z1, Z2, at_xx = TRUE)

    EA1_B1X2_x = EA1_x - EB1_x*X2

    
    for (i in seq(n_obs)) {

        # Invert
        b <- matrix(c(EY1sq[i], 
                    EY2sq[i] - EDY2sq[i] - 2*EA1_B1X2_x[i]*EDY2[i],
                    EY1Y2[i] - EY1[i]*EDY2[i]), 3, 1)
        A <- matrix(c(1, X1[i]^2, 2*X1[i],
                    1, X2[i]^2, 2*X2[i],
                    1, X1[i]*X2[i], X1[i] + X2[i]), 
                    byrow = TRUE, 3, 3)
        
        # Include in answer only if acceptable rcond
        if (rcond(A) > lb) { 
            cmoments[,i] <- Re(solve(A) %*% b)	
        } else {
            cmoments[,i] <- NaN
        }
    }
	
    return(cmoments)
}


outlier_to_nan <- function(x, q_low, q_high) {
    output <- x
    output[(x < quantile(x, q_low, na.rm = TRUE))] <- NA
    output[(x > quantile(x, q_high, na.rm = TRUE))] <- NA
    return(output)
}
    

compute_unconditional_moments <- 
    function(EA1_x, EB1_x, EA1sq_x, EB1sq_x, EA1B1_x) {
        
    means <- c(mean(EA1_x, na.rm = TRUE), mean(EB1_x, na.rm = TRUE))
    vars <- c(mean(EA1sq_x, na.rm = TRUE) - means[1]^2,
            mean(EB1sq_x, na.rm = TRUE) - means[2]^2)
    covar <- c(mean(EA1B1_x, na.rm = TRUE) - means[1]*means[2])
    return(list(random_coeff_means = means, 
                random_coeff_variances = vars, 
                random_coeff_covariance = covar))
}


# Estimates RCs under Normality assumptions
estimate_main <- function(Y1,Y2, X1, X2, Z1, Z2, 
                    mean_rcond_bnd, cov_rcond_bnd,
                    bw_lower, bw_upper, n_bws,
                    q1_low, q1_high, q2_low, q2_high) {

    
    # Shocks moments and fixed coeffs
    shocks_bw <- estimate_shocks_and_fixed_optimal_bw(Y1, Y2, X1, X2, Z1, Z2, bw_lower, bw_upper, n_bws)
    shocks_and_fixed <-  estimate_shocks_and_fixed(Y1, Y2, X1, X2, Z1, Z2, shocks_bw)
    b_Z <- shocks_and_fixed$fixed_coeffs
    n_z <- length(b_Z)/2

    # Subtract terms involving fixed coeffs (i.e., Z_{t}'b_{t})
    if (is.null(Z1) || is.null(Z2)) {
        Y1tilde <- Y1
        Y2tilde <- Y2
    } else {
        Y1tilde <- Y1 - Z1 %*% b_Z[c(1,n_z), 1, drop = FALSE]    
        Y2tilde <- Y2 -  Z2 %*% b_Z[c(n_z+1, 2*n_z), 1, drop = FALSE]
    }
    
    # Conditional first moments 
    conditional_means <- 
        estimate_conditional_first_moments(Y1tilde, Y2tilde, X1, X2, Z1, Z2, 
                                           lb = mean_rcond_bnd)

    # Trim outliers 
    EA1_x <- outlier_to_nan(x = conditional_means[1,], q1_low, q1_high)
    EB1_x <- outlier_to_nan(x = conditional_means[2,], q1_low, q1_high)

    # Conditional second moments
    conditional_second_moments <- 
        estimate_conditional_second_moments(Y1tilde, Y2tilde, X1, X2, Z1, Z2, 
            EA1_x = EA1_x, EB1_x = EB1_x, lb = cov_rcond_bnd)

   
    # Trim outliers
    EA1sq_x <- outlier_to_nan(x = conditional_second_moments[1,], q2_low, q2_high)
    EB1sq_x <- outlier_to_nan(x = conditional_second_moments[2,], q2_low, q2_high)
    EA1B1_x <- outlier_to_nan(x = conditional_second_moments[3,], q2_low, q2_high)
    
    # Compute centered, unconditional moments
    rc_moments <- compute_unconditional_moments(EA1_x, EB1_x, 
                                               EA1sq_x, EB1sq_x, EA1B1_x) 

    return(c(shocks_and_fixed, rc_moments))
}


# Main wrapper function
fhhps <- function(Y1, Y2, X1, X2, Z1, Z2,
                 mean_rcond_bnd = .1, cov_rcond_bnd = .1, 
                 bw_lower = .1, bw_upper = 1.5, n_bws = 20, 
                 q1_low = 0.01, q1_high = .99, q2_low = 0.00, q2_high = .98) {
     

    # Estimate moments
    moments <- estimate_main(Y1,Y2, X1, X2, Z1, Z2, 
                    mean_rcond_bnd, cov_rcond_bnd,
                    bw_lower, bw_upper, n_bws,
                    q1_low, q1_high, q2_low, q2_high) 
    
    return(moments)
}





#### Auxiliary function generates data satisfying model assumptions

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
                                         bw_lower, bw_upper, n_bws) {
    
    n_obs <- dim(Y1)[1]
    shocks_bws <- seq(bw_lower, bw_upper, length.out = n_bws)
    
    errors <- matrix(0, n_obs, n_bws)

    for (i in seq(n_obs)) {
        for (j in seq(n_bws)) { 
            shocks_and_fixed <- estimate_shocks_and_fixed(Y1 = Y1[-i,,drop=F], 
                                                         Y2 = Y2[-i,,drop=F], 
                                                         X1 = X1[-i,,drop=F], 
                                                         X2 = X2[-i,,drop=F], 
                                                         Z1 = Z1[-i,,drop=F], 
                                                         Z2 = Z2[-i,,drop=F], 
                                                         shocks_bw = shocks_bws[j])
            
            errors[i, j] <-  lp_error(YY = Y2[i] - Y1[i],
                                    XX = cbind(1, X2[-i,], (-1)*Z1[-i,], Z2[-i,]),
                                    xx = X2[i] - X1[i],
                                    betahat = c(shocks_and_fixed$shock_means,
                                                shocks_and_fixed$fixed_coeffs),
                                    bw = shocks_bws[j])
        }
    }
    mse <- apply(errors, 2, function(x) mean(x, na.rm = TRUE))
    return(shocks_bws[which.min(mse)])
}




fit_model <- function(Y, X1, X2, Z1, Z2, at_xx = FALSE) 
{
    if (is.null(Z1) || is.null(Z2)) {
        dataset <- data.frame(Y = Y, X1 = X1, X2 = X2)
        fmla <- formula(Y ~ s(X1) + s(X2))
    } else {
        dataset <- data.frame(Y = Y, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2)
        regressors <- foreach(col=colnames(dataset), .combine=c) %:% 
        when(col != "Y") %do% 
            sprintf("s(%s)", col) 
        fmla <- as.formula(paste(c("Y ~ ", paste(regressors, collapse = " + "))))
    }
    model <- gam(fmla, data = dataset)
    estimates <- loo_predict(model, dataset, at_xx)
    return(as.matrix(estimates))
}


loo_predict <- function(obj, dataset, at_xx) {
    registerDoParallel(cores = detectCores())
    print("Using this number of cores")
    print(getDoParWorkers())
    yhat <- foreach(i = 1:nrow(dataset), .combine = rbind) %dopar% {
        
        newdata = dataset[i,]
        if (at_xx) newdata["X1"] = newdata["X2"] 
        predict(update(object = obj, 
                       formula. = obj$formula,
                       data = dataset[-i, ]), 
                newdata = newdata, 
                type = "response")
    }
    return(data.frame(result = yhat[, 1], row.names = NULL))
}

