library(foreach)
library(doParallel)

n_values <- c(20, 50, 100)
M <- 500000
K <- 1000
block_size <- 50


TSs_fun <- function(n, ind) {
  set.seed(1000*n + ind)
  r <- as.integer(0.8*n)
  k <- 1:(r-1)
  TSs <- numeric(M)
  for (j in 1:M) {
    obs <- sort(rexp(n))
    S <- cumsum(obs[1:r])
    T_sum <- S[r] + (n-r)*obs[r]
    
    Hs <- -k*log(k) - (r-k)*log(r-k) + k*log(S[k]) + (r-k)*log(T_sum-S[k])
    TS <- -r*log(r) + r*log(T_sum) - min(Hs)
    
    TSs[j] <- max(TS, 0)
  }
  return(TSs)
}

dir.create("CV_simulations", showWarnings = FALSE)

jobs <- expand.grid(
  n = n_values,
  block = 1:(K/block_size)
)

registerDoParallel(48)

tictoc::tic()

files <- foreach(i = 1:nrow(jobs), .combine = c) %dopar% {
  n <- jobs$n[i]
  block <- jobs$block[i]
  
  indices <- ((block-1)*block_size+1):(block*block_size)
  
  TS_block <- matrix(NA, nrow = block_size, ncol = M)
  
  for (j in 1:block_size) {
    TS_block[j, ] <- TSs_fun(n, indices[j])
  }
  
  rownames(TS_block) <- paste0("run_", indices)
  
  filename <- paste0("TS_n=", n, "_block", block, ".rds")
  filepath <- file.path("CV_simulations", filename)
  
  saveRDS(TS_block, filepath, compress = FALSE)
  
  filename
}

exectime <- tictoc::toc(quiet = TRUE)

print(files)
print(exectime$toc - exectime$tic)

save.image("Sensitivity_CVs_simulation.RData")