#####Samples from an MVN Distribution
mvnsamp<-function(sig.mn,prec){
  n = length(sig.mn)
  prec = forceSymmetric(prec)
  L = t(chol(prec))
  v = forwardsolve(L,sig.mn)
  mu = backsolve(t(L),v)
  z = rnorm(n,0,1)
  y = backsolve(t(L),z)
  x = y + mu
  return(x)
}