rm(list = ls())
source("functions.r")

# Generate Data 
data <- create_data(n_obs = 1000)

# Simulation
n_iter = 1
for (k in seq(n_iter)) {

    output <- fhhps(X1 = data$X1, X2 = data$X2,
                   Z1 = data$Z1, Z2 = data$Z2,
                   Y1 = data$Y1, Y2 = data$Y2)
    
    write.table(rbind(unlist(output)), 
                file = "sim_results.txt", append = TRUE,
                row.names = FALSE, col.names = FALSE) 
}
