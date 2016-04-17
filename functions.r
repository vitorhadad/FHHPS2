library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(plyr)
library(caret)

#Get C++ functions (takes about a few moments to compile)
sourceCpp("cpp_functions.cpp")

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

    shock_variances <- shock_cov[c(1,2)]
    shock_covariances <- shock_cov[c(3)]

    return(list(shock_means = muW, 
                shock_variances = shock_cov[c(1,2)], 
                shock_covariance = shock_cov[3],
                fixed_coeffs = b_Z))
}


estimate_conditional_first_moments <- function(X, Y, Z, lb, bw1, bw2) {

    n_obs <- dim(Y)[1]
    XZ <- cbind(X, Z)
    n_regs <- dim(XZ)[2]

    results <- matrix(0, 2, n_obs)

    for (i in seq(n_obs)) {
        
        Y1 <- matrix(Y[-i,1], n_obs-1, 1)
        Y2 <- matrix(Y[-i,2], n_obs-1, 1)
        W <- matrix(cbind(XZ[-i,]), n_obs-1, n_regs)
        w <- matrix(c(XZ[i,]), 1, n_regs)
        w_xx <- matrix(c(X[i,2], X[i,2], Z[i,]), 1, n_regs)
         
        EY1 <- kreg_R(Y1, W, w, bw1)
        EY2 <- kreg_R(Y2, W, w, bw1)
        EDY2 <- kreg_R(Y2 - Y1, W, w_xx, bw2)
               

        b <- matrix(c(EY1, EY2 - EDY2), 2, 1)
        A <- matrix(c(1, 1, X[i,1], X[i,2]), 2, 2)

        results[,i] <- Re(solve(A) %*% b)	
    }
    
    if (rcond(A) > lb) { 
	results[,i] <- Re(solve(A) %*% b)	
    } else {
	results[,i] <- NaN
    }

    return(results)
}



estimate_conditional_second_moments <- 
    function(X, Y, Z, EA1_x, EB1_x, lb, bw1, bw2) {

    n_obs <- dim(Y)[1]
    XZ <- cbind(X, Z)
    n_regs <- dim(XZ)[2]
    
    results <- matrix(0, 3, n_obs)
    rconds <- rep(0, n_obs)
    for (i in seq(n_obs)) {

	# Create regressor matrices
        W <- matrix(cbind(XZ[-i,]), n_obs-1, n_regs)
        w <- matrix(c(XZ[i,]), 1, n_regs)
        w_xx <- matrix(c(X[i,2], X[i,2], Z[i,]), 1, n_regs)

	Y1 <- matrix(Y[-i,1], n_obs-1, 1)
	Y2 <- matrix(Y[-i,2], n_obs-1, 1)
	DY2 <- matrix(Y[-i,2] - Y[-i,1], n_obs-1, 1)

	# Conditional means
	EY1 <- kreg_R(Y1, W, w, bw1)	
	EDY2 <- kreg_R(DY2, W, w_xx, bw2)
        EY1sq <- kreg_R(Y1^2, W, w, bw1)	
        EY2sq <- kreg_R(Y2^2, W, w, bw1)
        EY1Y2 <- kreg_R(Y1*Y2, W, w, bw1)
        EDY2sq <- kreg_R(DY2^2, W, w_xx, bw2)

        EA1_B1X2_x = EA1_x[i] - EB1_x[i]*X[i,2]

        # Invert
        b <- matrix(c(EY1sq, 
                    EY2sq - EDY2sq - 2*EA1_B1X2_x*EDY2,
                    EY1Y2 - EY1*EDY2), 3, 1)
        A <- matrix(c(1, X[i,1]^2, 2*X[i,1],
                    1, X[i,2]^2, 2*X[i,2],
                    1, X[i,1]*X[i,2], X[i,1] + X[i,2]), 
                    byrow = TRUE, 3, 3)
        
        rconds[i] <- rcond(A)
        # Include in answer only if acceptable rcond
        if (rcond(A) > lb) { 
                results[,i] <- Re(solve(A) %*% b)	
        } else {
                results[,i] <- NaN
        }
    }
	
    return(results)
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
estimate_main <- function(X, Z, Y, 
                    mean_rcond_bnd, cov_rcond_bnd,
                    shocks_bw,
                    mean_bw1, mean_bw2, cov_bw1, cov_bw2,
                    q1_low, q1_high, q2_low, q2_high) {

    n_obs <- dim(Y)[1]
    if (is.null(Z)) {
        Z1 <- NULL
        Z2 <- NULL
    } else {
        n_z <- dim(Z)[2]/2
        Z1 <- Z[,seq(1, n_z)]
        Z2 <- Z[,seq(n_z+1, 2*n_z)]
    }
    
    # Shocks moments and fixed coeffs
    shocks_and_fixed <- estimate_shocks_and_fixed(Y, X, Z1, Z2, shocks_bw)
    b_Z = shocks_and_fixed$fixed_coeffs

    # Subtract terms involving fixed coeffs (i.e., Z_{t}'b_{t})
    if (is.null(Z)) {
        Ytilde <- Y
    } else {
        Z1b1 <- Z1 %*% b_Z[c(1,n_z), 1, drop = FALSE]
        Z2b2 <- Z2 %*% b_Z[c(n_z+1, 2*n_z), 1, drop = FALSE]
        Ytilde <- Y - cbind(Z1b1, Z2b2)
    }
    
    # Conditional first moments 
    conditional_means <- 
        estimate_conditional_first_moments(X, Ytilde, Z,
            lb = mean_rcond_bnd, 
            bw1 = mean_bw1, 
            bw2 = mean_bw2)

    # Trim outliers 
    EA1_x <- outlier_to_nan(x = conditional_means[1,], q1_low, q1_high)
    EB1_x <- outlier_to_nan(x = conditional_means[2,], q1_low, q1_high)

    # Conditional second moments
    conditional_second_moments <- estimate_conditional_second_moments(X, Ytilde, Z, 
            EA1_x = EA1_x, EB1_x = EB1_x,
            lb = cov_rcond_bnd, 
            bw1 = cov_bw1, 
            bw2 = cov_bw2)
   
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
                 shocks_bw = .1,
                 mean_bw1 = .1, mean_bw2 = .1, 
                 cov_bw1 = .1, cov_bw2 = .1,
                 q1_low = 0.01, q1_high = .99,
                 q2_low = 0.00, q2_high = .98) {
     
    # Ensure appropriate matrix form
    n_obs <- length(Y1)
    Y <- matrix(cbind(Y1, Y2), n_obs, 2)    
    X <- matrix(cbind(X1, X2), n_obs, 2)
    
    if (is.null(Z1) || is.null(Z2)) {
        Z <- NULL
    } else {
        Z <- matrix(cbind(Z1, Z2), n_obs, 2*dim(Z1)[2])
    }
    
    # Estimate moments
    moments <- estimate_main(X, Z, Y, mean_rcond_bnd, cov_rcond_bnd, shocks_bw, 
                mean_bw1, mean_bw2, cov_bw1, cov_bw2, 
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

estimate_shocks_and_fixed_CV <- function(Y1,Y2,X1,X2,Z1,Z2, 
                                         shocks_bws, n_folds = 5) {
    
    n_obs <- dim(Y1)[1]
    folds <- createFolds(y = seq(n_obs), k = n_folds)

    errors <- matrix(0, length(shocks_bws), n_folds)
    colnames(errors) <- names(folds)
    rownames(errors) <- shocks_bws

    for (i in seq_along(shocks_bws)) {
    bw <- shocks_bws[i]
        
        for (j in seq_along(folds)) {
            idx <- folds[[j]]
        
            Y1_test <- Y1[idx,,drop = F]
            Y2_test <- Y2[idx,,drop = F]
            X1_test <- X1[idx,,drop = F]
            X2_test <- X2[idx,,drop = F]
            Z1_test <- Z1[idx,,drop = F]
            Z2_test <- Z2[idx,,drop = F]
            Y1_train <- Y1[-idx,,drop = F]
            X1_train <- X1[-idx,,drop = F]
            Y2_train <- Y2[-idx,,drop = F]
            X2_train <- X2[-idx,,drop = F]
            Z1_train <- Z1[-idx,,drop = F]
            Z2_train <- Z2[-idx,,drop = F]

            shocks_and_fixed <- estimate_shocks_and_fixed(Y1 = Y1_train, 
                                                         Y2 = Y2_train, 
                                                         X1 = X1_train, 
                                                         X2 = X2_train, 
                                                         Z1 = Z1_train, 
                                                         Z2 = Z2_train, 
                                                         bw)
            
           errors[i, j] <-  lp_error(YY = Y2_test - Y1_test,
                                    XX = cbind(1, X2_test, (-1)*Z1_test, Z2_test),
                                    xx = X2_test - X1_test,
                                    betahat <- c(shocks_and_fixed$shock_means,
                                                shocks_and_fixed$fixed_coeffs),
                                     bw = bw)
        }
    }
    return(errors)
}


