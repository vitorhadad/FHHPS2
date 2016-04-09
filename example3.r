rm(list = ls())
source("functions.r")


# Parameters
n_obs <- 5000
q1_low <- 0.01
q1_high <- 0.99
mean_rcond_bnd <- .1
cov_rcond_bnd <- .1
bw0 <- .1
bw1 <- .1
bw2 <- .1
lb <- .1

ds <- create_data(n_obs, X_mean = 0)
Z1 <- ds$Z1
Z2 <- ds$Z2

umoments <- fhhps(X1 = ds$X[,1], X2 = ds$X[,2],
                Z1 = ds$Z1, Z2 = ds$Z2,
                Y1 = ds$Y[,1], Y2 = ds$Y[,2],
                mean_rcond_bnd, cov_rcond_bnd, 
                mean_bw1 = bw1, mean_bw2 = bw2, 
                cov_bw1 = bw1, cov_bw2 = bw2, 
                q1_low = q1_low, q1_high = q1_high, 
                q2_low = q1_low, q2_high = q1_high) 

print(umoments)
