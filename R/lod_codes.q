#' Estimates the parameter (nu or k) of the generalized Gamma distribution
#' @param par Starting value for the parameter
#' @param skx Skewness of the data for which the parameter is being estimated
#' @return Returns the estimated parameter nu or k 
#' @references
#' Stacy EW, Mihrama GA (1965) Parameter Estimation for a Generalized Gamma Distribution. Volume 7 issue 3 Pages 349-368. Based on 3rd equation of equation 21.
nu<-function(par,skx){
  tetra.gamma<-D(expression(trigamma(aa)),"aa")
  nu0<-exp(par)
  aa<-nu0
  val<-((eval(tetra.gamma)/(trigamma(nu0)^(3/2)))+abs(skx))^2
  return(val)
}


#' Estimates the parameters of the generalized Gamma distribution
#' @param x Vector of exposure observations
#' @references
#' Stacy EW, Mihrama GA (1965) Parameter Estimation for a Generalized Gamma Distribution. Volume 7 issue 3 Pages 349-368. Based on equation 21.
get.gengamma<-function(x){
  y<-log(x)
  # Compute moments - mean, variance and skewness
  muy<-mean(y)
  vary<-var(y)
  sky<-skewness(y)
  # starting value
  par1<-log(10)
  # Compute parameter estimate of the Generalized Gamma distribution
  # based on equation 21 in Stacey and Mihram (1965)
  nu0<-exp(optim(par1,fn=nu,skx=sky,method="BFGS",control=list(maxit=3000,parscale=abs(par1)))$par)
  p0<-sign(sky)*sqrt(trigamma(nu0))/sqrt(vary)
  a0<-exp(muy-(digamma(nu0)/p0))
  return(c(a0,p0,nu0))
}

#' Estimates the exposure value for observation below the limit of detection (LOD) assuming 
#' generalized gamma distribution for the exposure.
#' @param pars Parameters of the Generalized Gamma distribution - a vector of length 3 giving the 
#' parameters (scale,d,k) or equivalently (a,nu,p)
#' @param lod The value of the limit of detection
#' @references
#' Arunajadai SG., Rauh VA (2012) Handling covariates subject to limits of detection in regression. Volume 19 issue 3 Pages 369-391.
get.exposure<-function(pars,lod){
  fn1<-function(x){
    x*dgengamma(x,pars[1],pars[2],pars[3])
  }
  return(integrate(fn1,0,lod)$value/pgengamma(lod,pars[1],pars[2],pars[3]))
}
