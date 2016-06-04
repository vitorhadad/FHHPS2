rm(list = ls())
source("functions2.r")


# Generate Data 
n_obs_grid <- c(200, 1000, 2000)
n_tts <- c(10, 20, 30)

n_obs <- sample(n_obs_grid, 1)
n_tt <- sample(n_tts, 1)
det_bnds_grid <- seq(.01, .2, .01)
modulus_bnds_grid <- seq(.01, .2, .01)
mu_init = c(0,0)
sigma_init = matrix(c(1, 0, 0, 1), 2, 2)

det_bnd <- sample(det_bnds_grid, 1)
modulus_bnd <- sample(modulus_bnds_grid, 1)
reg_params <- c(1, 3, 5, 10, 20, 50)
filename = sprintf("sims.txt", n_obs)

# Generate Data 
dataset <- create_data(n_obs = n_obs, with_Z = TRUE)
   
Y1 <- dataset$Y1
X1 <- dataset$X1
Y2 <- dataset$Y2
X2 <- dataset$X2
Z1 <- dataset$Z1
Z2 <- dataset$Z2

table <- fhhps(Y1, Y2, X1, X2, Z1, Z2,
            mu_init = c(0,0), 
            sigma_init = matrix(c(1,0,0,1), 2,2),
            det_bnds = det_bnd,
            modulus_bnds = modulus_bnd,
            reg_params = reg_params,
            n_tt = n_tt)

print(table) 

for (line in seq(nrow(table))) {
    write.table(rbind(table[line,,drop=F]), 
            file = filename, append = TRUE,
            row.names = FALSE, col.names = if (file.exists(filename)) FALSE else TRUE,
            quote = FALSE) 
}


# Submit new job    
system("qsub job.pbs")




