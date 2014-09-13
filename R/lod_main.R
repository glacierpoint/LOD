#' Estimates the exposure value for observation below the limit of detection (LOD) assuming 
#' generalized gamma distribution for the exposure.
#' @import
#' @param data Data with two columns - first column is named exposure and second column is named status. 
#' Status takes on two values 0-detect and 1-non-detect i.e. below the limit of detection. The exposure values 
#' for the non-detects is the limit of detection
#' @param n.impute The required number of imputed data sets of the exposure data 
#' @param control List of control parameters that will be passed on to optim. See optim for details.
#' @return Returns a list with the following elements
#' \itemize{
#' \item convergence - Indicator if the optimization algorithm to estimate parameters of the 
#' generalized gamma distribution converged. 0 indicates convergence.
#' \item parameters - The parameter estimates of the generalized gamma distribution
#' \item se - The standard error estimates of the parameters of the generalized gamma distribution
#' \item imputed.data - A matrix of imputed values with number of columns equal to n.impute and number 
#' of rows equal to number of rows of data.
#' }
#' @examples
#' \dontrun{
#' 
#' my.lod<-lod(data10)
#' names(my.lod)         
#' 
#' plot(my.lod)
#' my.lod[my.lod$status==1,]
#' }
#' @references
#' Arunajadai SG., Rauh VA (2012) Handling covariates subject to limits of detection in regression. Volume 19 issue 3 Pages 369-391.
lod<-function(data,n.impute=10,control=list(maxit=2000)) UseMethod("lod")

lod.default<-function(data,n.impute=10,control=list(maxit=2000)){
  values<-data$exposure
  lod<-data$status
  detects<-values[which(lod==0)]
  non.detects<-values[which(lod==1)]
  cat("Computing starting values \n")
  pars<-log(unlist(get.gengamma(detects)))
  cat("     Starting values computed \n")
  fn<-function(pars,detects,non.detects){
    # compute the negative log-likelihood function
    pars<-exp(pars)
    lhood<-c(dgengamma(detects,pars[1],pars[2],pars[3]),pgengamma(non.detects,pars[1],pars[2],pars[3]))
    nllhood<-(-sum(log(lhood)))
    return(nllhood)
  }
  cat("Estimating parameters of Generalized Gamma Distribution \n")
  opt1<-try(optim(pars,fn=fn,detects=detects,non.detects=non.detects,hessian=TRUE,control=control),silent=FALSE)
  cat("     Parameters of Generalized Gamma Distribution estimated \n")
  par0<-opt1$par
  vcov0<-solve(opt1$hessian)
  par1<-exp(opt1$par)
  ses1<-deltamethod(list(~x1,~x2,~x3),mean=par0,cov=vcov0,ses=TRUE)
  names(par1)<-names(ses1)<-c("scale or a ","d or p ","k or nu")
  #  Get index of non.detects to store imputed values
  lod.ind<-split(1:nrow(data),data$status)[[2]]
  lod.ind<-split(lod.ind,non.detects)
  unique.lod<-as.numeric(unique(names(lod.ind)))
  impute.data<-matrix(0,ncol=n.impute,nrow=length(values))
  n<-length(unique.lod)
  count<-1
  cat("Imputing data \n")
  while(count<=n.impute){
    # get impuated parameters
    imp.pars0<-rmvnorm(1,par0,vcov0)
    imp.pars1<-exp(imp.pars0)
    exposure.lod<-values
    for(i in 1:n){
      mu<-try(get.exposure(imp.pars1,unique.lod[i]),silent=TRUE)
      if(class(mu)=="try-error" | is.nan(mu) | !is.finite(mu)){
        ind<-1
        break
      }else{
        ind<-0
        exposure.lod[lod.ind[[i]]]<-mu
      }
    }
    if(ind==1){
      ind<-0
      next
    }else{
      impute.data[,count]<-exposure.lod
      cat("    impute: ",count,"\n")
      count<-count+1
    }
  }
  cat("Data imputed \n")
  colnames(impute.data)<-paste("impute",1:n.impute,sep="")
  impute.data<-data.frame(data,impute.data)
  list1<-list(convergence=opt1$convergence,parameters=par1,se=ses1,imputed.data=impute.data)
  class(list1)<-"lod"
  return(list1)
}

#' Plots  the obsrved vs therotical quantiles based on the generalized gamma distribution 
#' generalized gamma distribution for the exposure.
#' @import
#' @examples
#' \dontrun{
#' 
#' my.lod<-lod(data10)
#' names(my.lod)         
#' 
#' plot(my.lod)
#' my.lod[my.lod$status==1,]
#' }
#' @references
#' Arunajadai SG., Rauh VA (2012) Handling covariates subject to limits of detection in regression. Volume 19 issue 3 Pages 369-391.
plot.lod<-function(obj){
  data1<-obj$imputed.data
  pars<-obj$parameters
  F1<-ecdf(data1$exposure)
  F1U<-ecdf(data1$exposure[data1$status==0])
  a1<-F1(data1$exposure)
  a2<-qgengamma(a1,pars[1],pars[2],pars[3])
  a3<-F1U(a2)
  a4<-F1U(data1$exposure)
  plot(a4,a3,xlab="Observed",ylab="Theoretical",pch=19,main="Exposures that are detectable")
  abline(0,1,col=2,lwd=2)
}


#' Simulated Exposure data -10\% below LOD
#'
#' A dataset containing simulated exposure for 
#' for 250 subjects with two different limits of detection with 10\% of the data below LOD
#'
#' @format A data frame with 259 rows and 2 variables - exposure and status
#' 
#' The variables are as follows:
#' \itemize{
#'   \item exposure Exposure measurements
#'   \item Status 0=detect  1= below limit of detection. Exposure values for status=1 is  the limit of detection for that measurement.
#' }
#' @name data10
NULL

#' Simulated Exposure data -25\% below LOD
#'
#' A dataset containing simulated exposure for 
#' for 250 subjects with two different limits of detection with 25\% of the data below LOD
#'
#' @format A data frame with 259 rows and 2 variables - exposure and status
#' 
#' The variables are as follows:
#' \itemize{
#'   \item exposure Exposure measurements
#'   \item Status 0=detect  1= below limit of detection. Exposure values for status=1 is  the limit of detection for that measurement.
#' }
#' @name data25
NULL

#' Simulated Exposure data -50\% below LOD
#'
#' A dataset containing simulated exposure for 
#' for 250 subjects with two different limits of detection with 50\% of the data below LOD
#'
#' @format A data frame with 259 rows and 2 variables - exposure and status
#' 
#' The variables are as follows:
#' \itemize{
#'   \item exposure Exposure measurements
#'   \item Status 0=detect  1= below limit of detection. Exposure values for status=1 is  the limit of detection for that measurement.
#' }
#' @name data50
NULL

#' Simulated Exposure data -75\% below LOD
#'
#' A dataset containing simulated exposure for 
#' for 250 subjects with two different limits of detection with 75\% of the data below LOD
#'
#' @format A data frame with 259 rows and 2 variables - exposure and status
#' 
#' The variables are as follows:
#' \itemize{
#'   \item exposure Exposure measurements
#'   \item Status 0=detect  1= below limit of detection. Exposure values for status=1 is  the limit of detection for that measurement.
#' }
#' @name data75
NULL

#' Simulated Exposure data -90\% below LOD
#'
#' A dataset containing simulated exposure for 
#' for 250 subjects with two different limits of detection with 90\% of the data below LOD
#'
#' @format A data frame with 259 rows and 2 variables - exposure and status
#' 
#' The variables are as follows:
#' \itemize{
#'   \item exposure Exposure measurements
#'   \item Status 0=detect  1= below limit of detection. Exposure values for status=1 is  the limit of detection for that measurement.
#' }
#' @name data90
NULL

#' Simulated Exposure data -95\% below LOD
#'
#' A dataset containing simulated exposure for 
#' for 250 subjects with two different limits of detection with 95\% of the data below LOD
#'
#' @format A data frame with 259 rows and 2 variables - exposure and status
#' 
#' The variables are as follows:
#' \itemize{
#'   \item exposure Exposure measurements
#'   \item Status 0=detect  1= below limit of detection. Exposure values for status=1 is  the limit of detection for that measurement.
#' }
#' @name data95
NULL
