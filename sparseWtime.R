

sparseWtime<-function(g.s,tm){
  W.base = sparseW(g.s)
  W.time = Matrix(matrix(0,g.s*g.s*tm,g.s*g.s*tm),sparse = TRUE)
  for(i in 1:tm){
    W.time[((i-1)*g.s*g.s+1):(i*g.s*g.s),((i-1)*g.s*g.s+1):(i*g.s*g.s)] = W.base
  }
  for(i in 1:(g.s^2)){
    W.time[i,(i+g.s^2)] = 1
  }
  for(i in (g.s^2 +1):((tm-1)*g.s^2)){
    W.time[i,i+g.s^2] = 1
    W.time[i,i-g.s^2] = 1
  }
  for(i in ((tm-1)*g.s^2 + 1):(tm*g.s*g.s)){
    W.time[i,i-g.s^2] = 1
  }
  return(W.time)
}
