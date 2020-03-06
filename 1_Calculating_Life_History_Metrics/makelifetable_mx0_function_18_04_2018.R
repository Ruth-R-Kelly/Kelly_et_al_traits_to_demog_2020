### function to calculate mx allowing adapted from Owen Jones et al.'s function in Mage 
### to include reproduction at age 0 (i.e. in the first year) . This is important
## for calculations with Annual plants

### adapted by Ruth Kelly 18 04 2018




makeLifeTable_mx0 <- function (matU, matF = NULL, matC = NULL, startLife = 1, nSteps = 1000) 
{
  matDim = ncol(matU)
  matUtemp = matU
  survivorship = array(NA, dim = c(nSteps, matDim))
  
  for (o in 1:nSteps) {
    survivorship[o, ] = colSums(matUtemp %*% matU)
    matUtemp = matUtemp %*% matU
  }
  
  lx = survivorship[, startLife]
  lx = c(1, lx[1:(length(lx) - 1)])
  out = data.frame(x = 0:(length(lx) - 1), lx = lx)
  if (!missing(matF)) {
    if (sum(matF, na.rm = T) == 0) {
      warning("matF contains only 0 values")
    }
    
    
    ageFertility = array(0, dim = c(nSteps, matDim))
    fertMatrix = array(0, dim = c(nSteps, matDim))
    matUtemp2 = matU
    e = matrix(rep(1, matDim))
    for (q in 1:nSteps) {
      fertMatrix = matF %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% 
                                                                 matUtemp2)))))
      ageFertility[q, ] = colSums(fertMatrix)
      matUtemp2 = matUtemp2 %*% matU
    }
    mx = ageFertility[, startLife]
    
    ### this next two lines are added to calculate reproduction in year 0.  
    ### I think the first is probably unnecessary and it will always resolve to 
    ## matF  i.e the other other terms resolve to 1. And fertility could be 
    ## directly calculated in year nought as the sum of the startLife column in matF. 
    ## So this here is an overkill just incase I've missed something approach.
    mx_time0 <- matF %*% (matU%^%0)* (as.numeric((ginv(diag(t(e) %*% 
                                                              (matU%^%0)))))) 
    mx_time0_sl <- colSums(mx_time0) [startLife]
    mx = c(mx_time0_sl, mx[1:(length(mx) - 1)])
    out$mx = mx
  }
  if (!missing(matC)) {
    if (sum(matC, na.rm = T) == 0) {
      warning("matC contains only 0 values")
    }
    ageClonality = array(0, dim = c(nSteps, matDim))
    clonMatrix = array(0, dim = c(nSteps, matDim))
    matUtemp2 = matU
    e = matrix(rep(1, matDim))
    for (q in 1:nSteps) {
      clonMatrix = matC %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% 
                                                                 matUtemp2)))))
      ageClonality[q, ] = colSums(clonMatrix)
      matUtemp2 = matUtemp2 %*% matU
    }
    cx = ageClonality[, startLife]
    cx_time0 <- matC %*% (matU%^%0)* (as.numeric((ginv(diag(t(e) %*% 
                                                              (matU%^%0))))))
    cx_time0_sl <- colSums(cx_time0) [startLife]
    cx = c(cx_time0_sl, mx[1:(length(mx) - 1)])
    out$cx = cx
  }
  return(out)
}