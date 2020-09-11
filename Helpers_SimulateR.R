## Contains the helper functions

runSimulation <- function(alpha,s,x,beta,seed){
  ## This function runs the simulation
  ## Inputs
  ## alpha: vector containing alpha parameter for each gene
  ## s: matrix containing scale values. Rows:genes and columns:samples
  ## x: design matrix
  ## beta: matrix containing GLM cofficients. Rows: genes columns samples 
  ## seed: Random number seed
  ## Output
  ## k: The counts matrix

  # Setting the random number seed
  set.seed(100)
  
  # Generating the scaled mean expression based on the GLM
  # log2(q) = xs %*% betas
  q <- t(2^(x %*% t(beta))) # Unscaled mean expression
  mu <- s * q # scaled mean expression
  
  # Generating the counts based on negative binomial distribution
  nrows <- nrow(mu)
  ncols <- ncol(mu)
  k <- matrix(nrow = nrows, ncol = ncols) # The counts matrix
  size <- 1 / alpha
  for (r in 1:nrows){
    for (c in 1:ncols){
      k[r,c] <- rnbinom(n = 1, 
                        size = size[r],
                        mu = mu[r,c])
    }
  }
  
  # Returning the counts
  return(k)
}

genDesignMat <- function(n.indv, n.biop, plot.heatmap){
  ## This function generates the design matrix
  ## Inputs
  ## n.indv: Number of indiviuals
  ## n.biop: Number of biopsies per indiviuals
  ## plot.heatmap: If TRUE plots the matrix as a heatmap
  ## Output
  ## x: Design matrix
  
  ## TO DO: Input validations
  
  x <- matrix(nrow = n.indv * n.biop, ncol = n.indv + 1)
  x[,] <- 0
  x[,1] <- 1 # intercept indication
  x[((n.biop*n.indv/2)+1):(n.biop*n.indv),2] <- 1 # Treatment indication
  
  ## Filing the individual indications
  for (i in 1:(n.indv - 1)){
    startIdx <- (i-1)* n.biop + n.biop + 1
    stopIdx <- startIdx + n.biop - 1
    x[startIdx:stopIdx,i+2] <- 1
  }
  
  # OPTIONAL: Plotting the design matrix as a heatmap
  if (plot.heatmap){
  heatmap(x,Rowv = NA, Colv = NA, revC = T, scale = "none",
          xlab = "Beta", ylab = "Indiviuals")
  }
  
  return(x)
}
