rm(list = ls())
source("functions2.r")


n_obs <- 200
n_tt <- 3

# Generate Data 
dataset <- create_data(n_obs = n_obs, with_Z = TRUE)
    
Y1 <- dataset$Y1
Y2 <- dataset$Y2
X1 <- dataset$X1
X2 <- dataset$X2
Z1 <- dataset$Z1
Z2 <- dataset$Z2

# Estimate
t <- proc.time()
output <- fhhps(X1 = dataset$X1, X2 = dataset$X2,
                Z1 = dataset$Z1, Z2 = dataset$Z2,
                Y1 = dataset$Y1, Y2 = dataset$Y2, 
                n_tt = n_tt)
t <- proc.time() - t

