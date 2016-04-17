rm(list = ls())
source("functions.r")

n_sims = 1000
n_obs = 500

# Generate Data 
data <- create_data(n_obs = n_obs, with_Z = TRUE)

X1 = data$X1
X2 = data$X2
Z1 = data$Z1
Z2 = data$Z2
Y1 = data$Y1
Y2 = data$Y2
shocks_bws = c(.01, .05, .1, .3, .5)
n_folds = 5

errors = estimate_shocks_and_fixed_CV(Y1, Y2, X1, X2, Z1, Z2, shocks_bws = shocks_bws)


# Estimate
#output <- fhhps(X1 = data$X1, X2 = data$X2,
#                Z1 = data$Z1, Z2 = data$Z2,
#                Y1 = data$Y1, Y2 = data$Y2)

    # Output
#write.table(rbind(unlist(output)), file = "sim_results.txt")

