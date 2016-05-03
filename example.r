rm(list = ls())
source("functions.r")

n_obs = 5000

# Generate Data 
dataset <- create_data(n_obs = n_obs, with_Z = TRUE)
X1 = dataset$X1
X2 = dataset$X2
Z1 = dataset$Z1
Z2 = dataset$Z2
Y1 = dataset$Y1
Y2 = dataset$Y2
bw_lower = .1
bw_upper = 1.2
n_bws = 9
lb = 0.1
q1_low = 0.01
q2_low = 0.00
q1_high = .99
q1_high = 0.98


# Generate Data 
data <- create_data(n_obs = n_obs, with_Z = TRUE)
    
# Estimate
output <- fhhps(X1 = dataset$X1, X2 = dataset$X2,
                Z1 = dataset$Z1, Z2 = dataset$Z2,
                Y1 = dataset$Y1, Y2 = dataset$Y2, 
                mean_rcond_bnd = lb, cov_rcond_bnd = lb, 
                bw_lower = .1, bw_upper = 1.2, n_bws = 9)

write.table(rbind(c(n_obs, lb, unlist(output))), 
            file = "sim_results.txt", append = TRUE,
            row.names = FALSE, col.names = FALSE,
            quote = FALSE) 


