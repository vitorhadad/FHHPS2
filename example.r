
source("functions.r")
    
r = c(.05,.1,.3,.5)
pgrid = expand.grid(r,r,r)

n_obs <- 1000



for (i in seq(dim(pgrid)[1])) {
    
    p = as.numeric(pgrid[i,])

    for (k in range(10)) {
        
        ds <- create_data(n_obs)
        moments <- fhhps(ds$Y[,1], ds$Y[,2], ds$X[,1], ds$X[,2], Z1 = NULL, Z2 = NULL,
                    mean_rcond_bnd = p[1], cov_rcond_bnd = p[1], 
                    mean_bw1 =  p[2], mean_bw2 = p[2], cov_bw1 = p[3], cov_bw2 = p[3]) 
        
        output <- rbind(c(n_obs, p, moments))
        
        write.table(output, file = "example_sims.txt", append = TRUE,
                    row.names = FALSE, col.names = FALSE)    
    }
    
}

