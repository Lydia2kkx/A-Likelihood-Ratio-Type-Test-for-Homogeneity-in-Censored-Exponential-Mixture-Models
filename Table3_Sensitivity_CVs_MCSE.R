
n_values <- c(20, 50, 100)
M_values <- c(10000, 50000, 100000, 500000)

K <- 1000
block_size <- 50
alpha <- 0.05

# Store estimated critical values
CVs <- array(NA, dim = c(length(n_values), K, length(M_values)))

# Read simulation results
for (i in seq_along(n_values)) {
  
  n <- n_values[i]
  
  for (block in 1:(K/block_size)) {
    
    filename <- paste0("TS_n=", n, "_block", block, ".rds")
    filepath <- file.path("CV_simulations", filename)
    
    TS_block <- readRDS(filepath)
    
    for (j in 1:block_size) {
      
      ind <- (block-1)*block_size + j
      
      for (k in seq_along(M_values)) {
        
        CVs[i, ind, k] <- quantile(
          TS_block[j, 1:M_values[k]],
          probs = 1-alpha,
          names = FALSE
        )
      }
    }
    
    rm(TS_block)
  }
}

# Calculate MCSEs and sensitivity
results <- NULL

for (i in seq_along(n_values)) {
  
  CV_ref <- CVs[i, , length(M_values)]
  
  for (k in seq_along(M_values)) {
    
    CV <- CVs[i, , k]
    
    result <- data.frame(
      n = n_values[i],
      M = M_values[k],
      Mean_CV = mean(CV),
      MCSE = sd(CV),
      Mean_Difference = mean(CV - CV_ref),
      Mean_Abs_Difference = mean(abs(CV - CV_ref))
    )
    
    results <- rbind(results, result)
  }
}

print(results)

save(CVs, results, file = "CV_MCSE_sensitivity_Type2.RData")

write.csv(results,
          file = "CV_MCSE_sensitivity_Type2.csv",
          row.names = FALSE)
