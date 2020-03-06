#### pSeedling 
## A function to calculate probability of reaching seedling stage, from seed stage. 
## Ruth Kelly 15/06/2017,  based on Caswell 2001 -2nd Ed. page 124-125
## Code adapted from 'lifeTimeRepEvents' function in the package Mage. 

pSeedling <- function(matU, seed1 = 1, seedling1 = 1){
  
  if(missing(matU)){stop('matU missing')}
  if(seed1 > seedling1){stop('seed1 must not be greater than seedling1')}
  #Function to determine probability of reaching reproduction, age at maturity and reproductive lifespan (Code adapted from H. Caswell's matlab code):
  
  ### calculate width/length of U matrix
  uDim <- dim(matU)[1]
  
  ### surv is the proportion surviving in eash column/stage of the U matrix
  surv <- colSums(matU)

  #### make vector of 0 for seedstage and 1's for plant stages 
  plant_stages <- rep(0, times = uDim)
  plant_stages[seedling1:uDim] <- 1
  
  #### Probability of survival to first seedling stage
  
  ### create a new absorption stage by representing the non-seeds stage as
  ### zeros in the transition matrix
  Uprime <- matU
  Uprime[,seedling1:uDim] <- 0

  ### create Mprime, 2 row matrix, same ncols as the U matrix 
  Mprime <- matrix(0,2,uDim)
  
for (p in 1:uDim) {
    if (plant_stages[p]==1) Mprime[2,p] = 1 else
      Mprime[1,p] = 1-surv[p]
  }


Bprime  <- Mprime%*%(ginv(diag(uDim)-Uprime))


### this set to count from the first stage (i.e. youngest seed stage)
prop_seedling <- Bprime[2,seed1]

return(prop_seedling)  
## 
 
}


