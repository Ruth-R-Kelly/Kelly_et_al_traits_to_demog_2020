## "exceptionalLife" function by Kevin Healy. 

## Editted by Ruth Kelly 06/06/2018 such that 99% lifespan is now 99.9% 


### Calculate from this how long things actually live by seeing how many are 
## still alive at time point x.  
## The function 'exceptionalLife', by Kevin Healy does this.

### I've modified the function so the arbitary maximum life is 10000, instead of 100

exceptionalLife<-function(matU,startLife=1){
  popVector=rep(0,dim(matU)[1])
  popVector[startLife]=10000
  lifespanLeftover=matrix(0,10000,1)
  for (n in 1:10000){
    lifespanLeftover[n]=sum(popVector)
    popVector=matU%*%popVector
  }
  Lexcept.50=min(which(lifespanLeftover<5000))
  if(Lexcept.50==Inf) {Lexcept.50=NA}
  Lexcept.95=min(which(lifespanLeftover<500))
  if(Lexcept.95==Inf) {Lexcept.95=NA}
  Lexcept.99=min(which(lifespanLeftover<100))
  if(Lexcept.99==Inf) {Lexcept.99=NA}
  Lexcept.999=min(which(lifespanLeftover<10))
  if(Lexcept.999==Inf) {Lexcept.999=NA}
  
  return(list(Lexcept.50 = Lexcept.50, Lexcept.95 = Lexcept.95, Lexcept.99 = Lexcept.99, Lexcept.999 = Lexcept.999))
}
