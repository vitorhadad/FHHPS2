rm(list = ls())
source("functions2.r")


# Generate Data 
n_obs_grid <- c(500, 1000, 2000)
det_bnds <- c(.01, .05, .1, .2)
modulus_bnds <- c(.01, .05, .1, .2)
n_tts <- c(7, 10, 15)
k <- 1
filename = "sim_results_with_cf.txt"

# Parameter grid
params <- expand.grid(n_obs_grid, det_bnds, modulus_bnds, n_tts)
# Shuffled
params <- params[sample(nrow(params)),]

while (TRUE) {

n_obs <- params[k,1]
det_bnd <- params[k,2]
modulus_bnd <- params[k,3]
n_tt <- params[k,4]

# Generate Data 
data <- create_data(n_obs = n_obs, with_Z = TRUE)
    
# Estimate
t <- proc.time()
output <- fhhps(X1 = dataset$X1, X2 = dataset$X2,
                Z1 = dataset$Z1, Z2 = dataset$Z2,
                Y1 = dataset$Y1, Y2 = dataset$Y2, 
                det_bnd = det_bnd, modulus_bnd = modulus_bnd, 
                n_tt = n_tt)
t <- proc.time() - t

write.table(rbind(c(n_obs, det_bnd, modulus_bnd, n_tt, unlist(output), t[3])), 
            file = filename, append = TRUE,
            row.names = FALSE, col.names = if (file.exists(filename)) FALSE else TRUE,
            quote = FALSE) 

k = k + 1
}
