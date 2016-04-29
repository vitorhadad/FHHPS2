rm(list = ls())
source("cv_functions.r")

n_obs = 50

# Generate Data 
dataset <- create_data(n_obs = n_obs, with_Z = TRUE)
X1 = dataset$X1
X2 = dataset$X2
Z1 = dataset$Z1
Z2 = dataset$Z2
Y1 = dataset$Y1
Y2 = dataset$Y2
lb = 0

bw_lower = .1
bw_upper = .2
n_bws = 1

# Estimate
output <- fhhps(X1 = dataset$X1, X2 = dataset$X2,
                Z1 = dataset$Z1, Z2 = dataset$Z2,
                Y1 = dataset$Y1, Y2 = dataset$Y2,
                bw_lower = .1, bw_upper = .2, n_bws = 2)


