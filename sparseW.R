
###Makes a spam adjacent matrix
sparseW <-function(g.s){
  S = g.s*g.s
  ###Everything Else
  int.seq = seq(from = (g.s+1),to = (S-g.s),by = 1)
  ind = int.seq%%g.s==0
  int.seq = int.seq[!ind]
  ind = int.seq%%g.s ==1 
  int.seq = int.seq[!ind]
  col.ind = unlist(lapply(as.matrix(int.seq),colfunc,gd = g.s))
  row.ind = unlist(lapply(as.matrix(int.seq),rowfunc))
  W = sparseMatrix(i = row.ind,j = col.ind,x = rep(1,length(col.ind)),dims = c(S,S))
  # for(i in int.seq){
  #   W[i,i-1] = 1
  #   W[i,i+1] = 1
  #   W[i,i-g.s-1] = 1
  #   W[i,i-g.s] = 1
  #   W[i,i-g.s+1] = 1
  #   W[i,i+g.s-1] = 1
  #   W[i,i+g.s] = 1
  #   W[i,i+g.s+1] = 1
  #   #print(i)
  # }
  ###Top Left Corner
  W[1,2] = 1
  W[1,1+g.s] = 1
  #W[1,2 + g.s] = 1
  ###Top Right Corner
  W[g.s,g.s-1] = 1
  #W[g.s,2*g.s-1] = 1
  W[g.s,2*g.s] = 1
  ###Bottom Left Corner
  W[((g.s-1)*g.s + 1),((g.s-2)*g.s + 1) ] = 1
  #W[((g.s-1)*g.s + 1),((g.s-2)*g.s + 2) ] = 1
  W[((g.s-1)*g.s + 1),((g.s-1)*g.s + 2) ] = 1
  ###Bottom Right Corner
  W[S,S-1] = 1
  W[S, S -g.s] = 1
  #W[S, (S - g.s -1)] = 1
  ###Top Row
  for(i in 2:(g.s-1)){
    W[i,i-1] = 1
    W[i,i+1] = 1
    #W[i,i-1 + g.s] = 1
    W[i,i + g.s] = 1
    #W[i,i + g.s +1] = 1
  }
  ###Bottom Row
  for(i in (g.s*(g.s-1)+ 2):(S-1)){
    W[i,i-1] = 1
    W[i,i+1] = 1
    #W[i,i-1 - g.s] = 1
    W[i,i - g.s] = 1
    #W[i,i - g.s +1] = 1
  }
  ###Left Column
  for(i in 1:(g.s-2)){
    W[(i*g.s + 1),((i-1)*g.s+1)] = 1
    #W[(i*g.s + 1),((i-1)*g.s+2)] = 1
    W[(i*g.s + 1),(i*g.s + 2)] = 1
    W[(i*g.s + 1),((i+1)*g.s+1)] = 1
    #W[(i*g.s + 1),((i+1)*g.s+2)] = 1
  }
  ###Right Column
  for(i in 2:(g.s-1)){
    W[(i*g.s ),((i-1)*g.s )] = 1
    #W[(i*g.s),((i-1)*g.s -1)] = 1
    W[(i*g.s),(i*g.s-1)] = 1
    W[(i*g.s),((i+1)*g.s)] = 1
    #W[(i*g.s),((i+1)*g.s -1)] = 1
  }
  return(W)
}