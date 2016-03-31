rm(list = ls())
source("functions.r")


# Parameters
n_obs <- 1000
q1_low <- 0
q1_high <- 1
mean_rcond_bnd <- .1
cov_rcond_bnd <- .1
bw1 <- .1
bw2 <- .1

ds <- create_data(n_obs)
Y <- ds$Y
X <- ds$X
Z <- NULL
n_obs <- dim(Y)[1]

moments <- fhhps(X1 = X[,1], X2 = X[,2],
                Z1 = Z[,1], Z2 = Z[,2],
                Y1 = Y[,1], Y2 = Y[,2],
                mean_rcond_bnd, cov_rcond_bnd, 
                mean_bw1 = bw1, mean_bw2 = bw2, 
                cov_bw1 = bw1, cov_bw2 = bw2, 
                q1_low = q1_low, q1_high = q1_high, 
                q2_low = q1_low, q2_high = q1_high) 

print(moments)
