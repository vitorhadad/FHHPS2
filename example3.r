rm(list = ls())
source("functions.r")

n_sims = 1000
n_obs = 1000

print("Simulation Example")

for (k in seq(n_sims)) {

    # Generate Data 
    data <- create_data(n_obs = n_obs, with_Z = TRUE)
    
    # Estimate
    output <- fhhps(X1 = data$X1, X2 = data$X2,
                   Z1 = data$Z1, Z2 = data$Z2,
                   Y1 = data$Y1, Y2 = data$Y2)
    
    sprintf("Simulation %d of %d", k, n_sims)

    # Output
    write.table(rbind(unlist(output)), 
                file = "sim_results.txt", append = ifelse(k == 1, FALSE, TRUE),
                row.names = FALSE, col.names = ifelse(k == 1, TRUE, FALSE),
                quote = FALSE) 
}


