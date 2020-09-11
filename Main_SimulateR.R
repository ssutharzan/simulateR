#### Created by Sreeskandarajan Sutharzan (University of Michigan)##
#### on 09/24/19 ###################################################

#### Setting up the environment
rm(list = ls())
setwd("") # Change this to your working directory
source("code/Helpers_SimulateR.R") # The location of the helper function


### The parameter choices
## Setting the seed
set.seed(10)

## Setting the parameters
n.indv <- 20 #Number of individuals
n.biop <- 4 #Number of biopsies
n.genes <- 5000 #Number of genes
n.genes.up <- round(n.genes*0.05) #Number of upregulated genes
n.genes.down <- round(n.genes*0.05) #Number of downregulated genes
n.genes.noreg <- n.genes-(n.genes.up+n.genes.down)
beta.intercept.params <- c(2,1) #Beta intercept mean and SD
beta.up.params <- c(1,1) #Upregulated beta mean and SD
beta.down.params <- c(1,1) #Downregulated beta mean and SD
beta.noreg.params <- c(0,0) #Nonregulated beta mean and SD
beta.indv.params.mean <- c(0,0.1) #Indiviual beta mean mean and SD
beta.indv.params.sd <- c(0,0.05) #Indiviual beta variance mean and SD
alpha.params <- c(0,0.01) # Alpha mean and SD
s <- 1 #Scale

## Generating the design matrix
x <- genDesignMat(n.indv = n.indv, n.biop = n.biop, plot.heatmap = T)

## alpha: vector containing alpha parameter for each gene
alpha <- abs(rnorm(n.genes,alpha.params[1],alpha.params[2]))
hist(alpha)

# beta.indv.params: matrix containing beta parameters for the indiviuals 
# Rows:genes and columns: 1-mean 2-sd
beta.indv.params <- cbind(rnorm
                          (n.genes,beta.indv.params.mean[1],
                            beta.indv.params.mean[2]),
                          sqrt(abs(rnorm(n.genes,beta.indv.params.sd[1],
                                         beta.indv.params.sd[2]))))

# generating the beta cofficients for the indiviuals
# beta.indv: matrix containing beta coefficients for the indiviuals 
# Rows:genes and columns: beta coffiecients
beta.indv <- t(sapply(1:n.genes,
                      function(indv) 
                        rnorm(
                          n.indv-1,beta.indv.params[indv,1],
                          beta.indv.params[indv,2])))
hist(beta.indv)

# beta.major: matrix containing the major beta cofficients
# Rows:genes and columns: beta coffiecients
beta.intercept <- rnorm(n.genes,beta.intercept.params[1],
                        beta.intercept.params[2])
hist(beta.intercept)


beta.up <- abs(rnorm(n.genes.up,beta.up.params[1],
                     beta.up.params[2]))
hist(beta.up)
beta.down <- -abs(rnorm(n.genes.down,beta.down.params[1],
                        beta.down.params[2]))
hist(beta.down)
beta.noreg <- rnorm(n.genes.noreg,beta.noreg.params[1],
                    beta.noreg.params[2])
hist(beta.noreg)
beta.treat <- c(beta.up,beta.down,beta.noreg)
hist(beta.treat)

# beta: matrix containing GLM cofficients. Rows: genes columns cofficients 
beta <- cbind(beta.intercept,beta.treat,beta.indv)

## Running the simulation
seed <- 100 # Random number seed

k <- runSimulation(alpha = alpha, s = s,
                   x = x, beta = beta, seed = seed)
hist(log(k),100)

## Writing the outputs
# Count data
gene.names <- sapply(1:n.genes, function(x) paste0("gene",x))
rownames(k) <- gene.names
write.csv(k,file="data/simT_counts_new.csv",
          row.names = T)

# Col data
patient <- as.factor(
  unlist(lapply(as.list(1:n.indv), function(x) rep(paste0("patient",x),n.biop))))
control <- rep("control",n.biop*n.indv/2)
treated <- rep("case",n.biop*n.indv/2)
treatment <- as.factor(c(control,treated))
coldat.names <- sapply(1:(n.indv*n.biop), function(x) paste0("sample",x))
col.dat <- data.frame(patient,treatment)
row.names(col.dat) <- coldat.names
write.csv(col.dat,file="data/simT_coldat.csv",
          row.names = T)