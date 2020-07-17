zsamp <-function(i,lambda.i,chi.i,psi.i){
  res = rgig(n = 1, lambda = lambda.i, chi = chi.i[i],psi = psi.i )
  return(res)
}