
###Format the Data
main.df = big.pm.lst.wint[[1]]
for(i in 2:length(big.pm.lst.wint)){
  main.df = rbind(main.df,big.pm.lst.wint[[i]])
}
date.list = sort(unique(main.df$Date))
t = length(date.list)
main.df$time <-NA
for(i in 1:length(date.list)){
  ind = main.df$Date == date.list[i]
  main.df$time[ind] = i
}
###Divide the data by year
main.df$year <-NA
###year 2010-2011
t.ind = c(60:149)
length(t.ind)
date.list[t.ind]
year.ind = is.element(main.df$time,t.ind)
main.df$year[year.ind] = 1
###year 2011-2012
t.ind = c(150:239)
length(t.ind)
date.list[t.ind]
year.ind = is.element(main.df$time,t.ind)
main.df$year[year.ind] = 2
###year 2012-2013
t.ind = c(241:330)
length(t.ind)
date.list[t.ind]
year.ind = is.element(main.df$time,t.ind)
main.df$year[year.ind] = 3
###year 2013-2014
t.ind = c(331:420)
length(t.ind)
date.list[t.ind]
year.ind = is.element(main.df$time,t.ind)
main.df$year[year.ind] = 4
ind.rm = is.na(main.df$year)
main.df = main.df[!ind.rm,]

###Add in time
date.list = sort(unique(main.df$Date))
t = length(date.list)
t.season = 90
main.df$time <-NA
for(i in 1:length(date.list)){
  ind = main.df$Date == date.list[i]
  main.df$time[ind] = i
}

###Remove Stations without much data
S.total = length(unique(main.df$stnID))
station.list = unique(main.df$stnID)
S.data = c()
S.range = c()
for(i in 1:S.total){
  ind = main.df$stnID == station.list[i]
  S.data = c(S.data, sum(ind))
  S.range = cbind(S.range,range(main.df$time[ind]))
  plot(main.df$time[ind],main = paste(i))
}
# Take out stations with big data gaps, plus two that're far away from the rest
#133
stations.rm = c(1,5,20:29)
stations.rm = c(stations.rm,which(S.range[1,] > 20),which(S.range[2,]<340) )
for(i in stations.rm){
  station.id = station.list[i]
  ind = main.df$stnID == station.id
  main.df = main.df[!ind,]
}
station.list = unique(main.df$stnID)

###Parameters
g.s = 15
Q = g.s*g.s
n.knots.c = 5
dg.c = 3
n.basis.c = n.knots.c + dg.c
n.phi = n.basis.c*n.basis.c
dg.re = 3
n.knots.re = 12
n.basis.re = dg.re + n.knots.re
Q.re = n.basis.re*n.basis.re*n.basis.re
S = length(station.list)
yr = 4
P = 12
r = 0.5

###Observations Locations
obs.loc.df = unique(cbind(main.df$lat,main.df$lon))
x.obs.loc = obs.loc.df[,2]
y.obs.loc = obs.loc.df[,1]
map('state')
points(x.obs.loc, y.obs.loc)

###Grid Cell Locations
x.grid = rep(seq(from = min(x.obs.loc)-1, to = max(x.obs.loc)+1, length.out = g.s),g.s)
y.grid = c()
y.off.set = (max(y.obs.loc) +1 - (min(y.obs.loc)-1))/g.s
for(i in 1:g.s){
  y.grid = c(y.grid,rep( min(y.obs.loc)-1 + y.off.set*2*(i-1), g.s))
}
points(x.grid,y.grid,col = 'red')
W = sparseW(g.s)
D = sparseMatrix(i = seq(1,Q),j = seq(1,Q),x = colSums(W))
#M = (D-W)%*%(D-W)
M = (D-W)
I = sparseMatrix(i = seq(1,Q),j = seq(1,Q),x = 0.0001)

###Knot locations (for coefficients)
x.knots = seq(min(x.grid),max(x.grid),length.out = n.knots.c)
y.knots = seq(min(y.grid),max(y.grid),length.out = n.knots.c)
B.x = bs(x = x.obs.loc,knots = x.knots, degree = dg.c)
B.y = bs(x = y.obs.loc, knots = y.knots, degree = dg.c)
B = matrix(NA,S,n.phi)
ct = 1
for(i in 1:n.basis.c){
  for(j in 1:n.basis.c){
    B[,ct] = B.x[,i]*B.y[,j]
    ct = ct +1
  }
}
B.single = Matrix(B,sparse = TRUE)
B.Mat = B.single
ct = 1
for(i in 1:t.season){
  for(j in 1:yr){
    if(ct >1){
      B.Mat = Matrix(rbind(B.Mat,B.single),sparse = TRUE)
    }
    ct = ct +1
  }
}
W = sparseW(g.s)
D = sparseMatrix(i = seq(1,Q),j = seq(1,Q),x = colSums(W))
M = (D - W)
W.phi = sparseW(n.basis.c)
D.phi = sparseMatrix(i = seq(1,n.phi),j = seq(1,n.phi),x = colSums(W.phi))
M.phi = (D.phi - W.phi)

###Knot locations
x.knots.re = seq(min(x.grid),max(x.grid),length.out = n.knots.re)
y.knots.re = seq(min(y.grid),max(x.grid),length.out = n.knots.re)
t.knots.re = seq(1,t.season,length.out = n.knots.re)
t.list = c()
for(i in 1:t.season){
  t.list = c(t.list,rep(i,S))
}
B.x.re = bs(x = rep(x.obs.loc,t.season),knots = x.knots.re, degree = dg.re)
B.y.re = bs(x = rep(y.obs.loc,t.season), knots = y.knots.re, degree = dg.re)
B.t.re = bs(x = t.list, knots = t.knots.re, degree = dg.re)

B = matrix(NA,S*t.season,Q.re)
ct = 1
for(i in 1:n.basis.re){
  for(j in 1:n.basis.re){
    for(k in 1:n.basis.re){
      B[,ct] = B.x.re[,i]*B.y.re[,j]*B.t.re[,k]
      ct = ct +1
    }
  }
}
B.Mat.re= matrix(0,S*t.season*yr,Q.re*yr)
for(j in 1:yr){
  B.Mat.re[((j-1)*S*t.season + 1):(j*S*t.season),((j-1)*Q.re+1):(j*Q.re)] = B
}
B.Mat.re = Matrix(B.Mat.re, sparse = TRUE)
W.psi = sparseWtime(n.basis.re,n.basis.re)
D.psi = sparseMatrix(i = seq(1,Q.re), j = seq(1,Q.re), x = colSums(W.psi))
M.psi = (D.psi- W.psi)

###Assign locations to grid cells
grid.assign = c()
for(i in 1:S){
  loc.i = (x.obs.loc[i] - x.grid)^2 + (y.obs.loc[i] - y.grid)^2
  ind = which(loc.i == min(loc.i))
  grid.assign = c(grid.assign, ind)
}
grid.assign.overall = rep(grid.assign,yr*t.season)

###Determing which location time pairs are observed
obs.df = array(NA, c(S, t.season, yr))
eta.df = array(NA, c(S,t.season,yr))
X.array = array(99,c(S,t.season,yr,P))
for(s in 1:S){
  station.id = station.list[s]
  for(i in 1:t.season){
    for(j in 1:yr){
      ind = which(main.df$stnID == station.id & main.df$time == i + (j-1)*t.season & main.df$year == j )
      if(length(ind)>0){
        obs.df[s,i,j] = 1
        eta.df[s,i,j] = main.df$PM[ind]
        X.array[s,i,j,1] = main.df$wndspd[ind]
        X.array[s,i,j,2] = main.df$air.afternoon[ind]
        X.array[s,i,j,3] = main.df$air.night[ind]
        X.array[s,i,j,4] = main.df$dswrf[ind]
        X.array[s,i,j,5] = main.df$tcdc[ind]
        X.array[s,i,j,6] = main.df$prcp[ind]
        X.array[s,i,j,7] = main.df$rh[ind]
        X.array[s,i,j,8] = main.df$hpbl.night[ind]
        X.array[s,i,j,9] = main.df$hpbl.afternoon[ind]
        X.array[s,i,j,10] = main.df$lftx[ind]
        X.array[s,i,j,11] = main.df$lts[ind]
        X.array[s,i,j,12] = main.df$tke[ind]
      }
    }
  }
  print(s)
}
#save.image(paste0(path,"PreImage.Rdata"))

###Create the true function
theta = (1-2*r)/(r*(1-r))
lambda = sqrt(2/(r*(1-r)))

###Format Data for Loop
Y = c()
for(j in 1:yr){
  for(i in 1:t.season){
    Y = c(Y,eta.df[,i,j])
  }
}
Y.ind = is.na(Y)
Y = Y[!Y.ind]
B.Mat.re = B.Mat.re[!Y.ind,]
n = length(Y)
B.Mat = B.Mat[!Y.ind,]
ind.count = c()
for(j in 1:yr){
  for(i in 1:t.season){
    ind.count = c(ind.count, sum(!is.na(eta.df[,i,j])))
  }
}
ind.yr.list = list()
start.ind = 1
ct = 1
for(j in 1:yr){
  ind.t.list = list()
  for(i in 1:t.season){
    stop.ind = start.ind + ind.count[ct] -1
    ind.t.list[[i]] = seq(from = start.ind, to = stop.ind)
    start.ind = start.ind + ind.count[ct]
    ct = ct +1
  }
  ind.yr.list[[j]] = ind.t.list
}
X.list = list()
for(k in 1:P){
  x.vec = c()
  for(j in 1:yr){
    for(i in 1:t.season){
      x.vec = c(x.vec,X.array[,i,j,k])
    }
  }
  x.vec = x.vec[!Y.ind]
  X.list[[k]] = sparseMatrix(i = seq(1,n),j = seq(1,n),x = x.vec,dims = c(n,n))
}
grid.assign.overall = grid.assign.overall[!Y.ind]

# G.col = c()
# for(j in 1:yr){
#   for(i in 1:t.season){
#     ind.t = obs.df[,i,j] == 1
#     id.na = is.na(ind.t)
#     ind.t[id.na] = FALSE
#     G.col = c(G.col, grid.assign[ind.t] + (i-1)*Q + (j-1)*Q*t.season)
#   }
# }
# G.mat = sparseMatrix(i = seq(1,n), j = G.col, x = rep(1,n),dims = c(n,Q*t.season*yr))

###MCMC Initialization
G = 20000
burn = 10000
phi.g = matrix(0,n.phi,P)
tau.g = rep(1,P)
alpha.tau = 1
beta.tau = 1
upsilon.g = 1
sig.upsilon = 1000
psi.g = rep(0,Q.re*yr)
psi.array.g = array(0,c(Q.re,yr))
gamma.g = rep(1,yr)
alpha.gamma = 1
beta.gamma = 1
Z.g = rep(1,n)
#Z.g = c()
#for(i in 1:t){
#  for(j in 1:yr){
#    ind.z = !is.na(eta.true[,i,j])
#    Z.g = c(Z.g, Z.true[ind.z,i,j])
#  }
#}
Z.inv.g = 1/Z.g
Z.inv.g.D = sparseMatrix(i = seq(1,n),j = seq(1,n),x = Z.inv.g)
eta.g = rep(0,n)
###Storage Devices
phi.array = array(99,c(n.phi,P,G))
tau.array = array(99,c(P,G))
upsilon.array = array(99,c(1,G))
psi.array = array(99,c(Q.re*yr, G))
zeta.array = array(99, c(yr,G))
gamma.array = array(99, c(yr, G))
Z.array = array(99, c(n,G))

st.time = Sys.time()
for(g in 1:G){
  ###Sample Phi and Tau
  for(i in 1:P){
    eta.g = upsilon.g*rep(1,n) + B.Mat.re%*%psi.g + theta*Z.g
    for(j in 1:P){
      if(i!=j){
        eta.g = eta.g + X.list[[j]]%*%B.Mat%*%as.matrix(phi.g[,j])
      }
    }
      phi.mean = (1/lambda^2)*t(X.list[[i]]%*%B.Mat)%*%(Z.inv.g*(Y - eta.g))
      phi.sig = (1/lambda^2)*(t(X.list[[i]]%*%B.Mat)%*%(Z.inv.g.D%*%X.list[[i]]%*%B.Mat)) + (1/tau.g[i])*M.phi
      phi.g[,i] = mvnsamp(phi.mean,phi.sig)
      phi.array[,i,g] = phi.g[,i]
      tau.g[i] = rinvgamma(1,n.phi/2 + alpha.tau,rate = 1/2*sum(phi.g[,i]*(M.phi%*%phi.g[,i])) + beta.tau)
      tau.array[i,g] = tau.g[i]                     
  }
  ###Sample Upsilon
  eta.g = B.Mat.re%*%psi.g + theta*Z.g
  for(j in 1:P){
    eta.g = eta.g + X.list[[j]]%*%B.Mat%*%as.matrix(phi.g[,j])
  }
  upsilon.sig = 1/((1/lambda^2)*sum(Z.inv.g) + 1/sig.upsilon)
  upsilon.mean = upsilon.sig*((1/lambda^2)*sum(Z.inv.g.D%*%(Y - eta.g)))
  upsilon.g = rnorm(1,upsilon.mean,sqrt(upsilon.sig))
  upsilon.array[1,g] = upsilon.g
  ###Sample the psi's
  eta.g = upsilon.g + theta*Z.g
  for(j in 1:P){
    eta.g = eta.g + X.list[[j]]%*%B.Mat%*%as.matrix(phi.g[,j])
  }
  for(j in 1:yr){
    ind.yr = unlist(ind.yr.list[[j]])
    psi.mean = (1/lambda^2)*t(B.Mat.re[ind.yr,((j-1)*Q.re + 1):(j*Q.re)])%*%Z.inv.g.D[ind.yr,ind.yr]%*%(Y[ind.yr] - eta.g[ind.yr]) 
    psi.sig = (1/lambda^2)*t(B.Mat.re[ind.yr,((j-1)*Q.re + 1):(j*Q.re)])%*%Z.inv.g.D[ind.yr,ind.yr]%*%B.Mat.re[ind.yr,((j-1)*Q.re + 1):(j*Q.re)] + (1 /gamma.g[j])*M.psi
    psi.array.g[,j] = mvnsamp(psi.mean,psi.sig)
  }
  psi.g = c()
  for(j in 1:yr){
    psi.g = c(psi.g,psi.array.g[,j])
  }
  pm = mean(psi.g)
  psi.g = psi.g - pm
  psi.array.g = psi.array.g - pm
  ###Sample the Gamma's
  for(j in 1:yr){
    gamma.power = Q.re/2 + alpha.gamma
    gamma.exp = psi.array.g[,j]%*%M.psi%*%psi.array.g[,j]
    gamma.exp = (1/2)*gamma.exp + beta.gamma
    gamma.g[j] = rinvgamma(1,gamma.power,rate = as.numeric(gamma.exp))
  }
  ###Sample the Z's
  eta.g = upsilon.g*rep(1,n) + B.Mat.re%*%psi.g
  for(j in 1:P){
    eta.g = eta.g + X.list[[j]]%*%B.Mat%*%as.matrix(phi.g[,j])
  }
  Z.chi = (Y - eta.g)^2/(lambda^2)
  Z.g = apply(matrix(seq(1,n),n,1),1,zsamp,lambda.i = 1/2, chi.i = Z.chi, psi.i = 2 + theta^2/(lambda^2))
  Z.inv.g = 1/Z.g
  Z.inv.g.D = sparseMatrix(i = seq(1,n),j = seq(1,n),x = Z.inv.g)
  psi.array[,g] = psi.g
  gamma.array[,g] = gamma.g
  Z.array[,g] = Z.g
  print(g)
}
sp.time = Sys.time()

#save.image(paste0(path,"All cov Winter50_3D_CAR_Dg3_AND_Dg3_SVC_15_knots"))

###Check convergence
plot(phi.array[1,1,1:G])
plot(phi.array[20,1,1:G])
plot(phi.array[10,1,1:G])
plot(phi.array[30,1,1:G])

plot(tau.array[1,1:G])
plot(tau.array[2,1:G])

plot(upsilon.array[1,1:G])

plot(psi.array[1,1:G])
plot(psi.array[2000,1:G])
plot(psi.array[400,1:G])
plot(psi.array[300,1:G])

plot(gamma.array[1,1:G])
plot(gamma.array[2,1:G])
plot(gamma.array[3,1:G])
plot(gamma.array[4,1:G])

###Knot locations
B.x.d = bs(x = x.grid,knots = x.knots, degree = dg.c)
B.y.d = bs(x = y.grid, knots = y.knots, degree = dg.c)
B.d = matrix(NA,Q,n.phi)
ct = 1
for(i in 1:n.basis.c){
  for(j in 1:n.basis.c){
    B.d[,ct] = B.x.d[,i]*B.y.d[,j]
    ct = ct +1
  }
}
phi.array.dis = array(NA,c(Q,P,G))
for(g in 1:G){
  phi.array.dis[,,g] = B.d%*%phi.array[,,g]
}

phi.hat = apply(phi.array.dis[,,burn:G],c(1,2),mean)
phi.lower = apply(phi.array.dis[,,burn:G],c(1,2),quantile,prob = 0.025)
phi.upper = apply(phi.array.dis[,,burn:G],c(1,2),quantile,prob = 0.975)
phi.sig = matrix(0,Q,P)
ind = phi.lower >= 0
phi.sig[ind] = 1
ind = phi.upper<=0
phi.sig[ind] = -1

###Make the Maps
library(rgdal)
library(raster)
library(RColorBrewer)
#path = "C:/Users/Travis/Stella Self Dropbox/Stella Watson/Synced Files/Research/GMRF with Unknown Structures/Quantile Regression Spatial Smoother Paper/Spline Redo/Data Apps/Winter50/3d car dg 3 AND dg 3 SVC/"
m1 = readOGR(paste('C:/Users/scwatson/Downloads/states_21basic'))
ind = is.element(m1$SUB_REGION,c("New England", "Middle Atlantic", "South Atlantic", 'East North Central','East South Central' ))
m1 = m1[ind,]
plot.names = c("Windspeed", "Daytime Air Temp.", "Nighttime Air Temp.", "DSWRF", "Percent Cloud Cover", "Precipitation", "Relative Humidity", "Nighttime HPBL", "Daytime HPBL", "LFTX", "LTS", "TKE")
#windows(width = 20,height = 18)
#par(mfrow = c(2,4))
###Make the Effect Surface Maps
for(p in 1:P){
  base = raster(xmn = min(x.grid),xmx = max(x.grid),ymn = min(y.grid),ymx = max(y.grid),nrow = g.s, ncol = g.s)
  base$val <- NA
  base$val.neg <- NA
  base$val.pos <- NA
  base$val.sig<- NA
  ind.cells = extract(base,y = cbind(x.grid,y.grid),cellnumbers = TRUE)[,1]
  ind.neg = phi.hat[,p] <= 0
  ind.pos = phi.hat[,p] >= 0
  phi.hat.neg = phi.hat[,p]
  phi.hat.neg[ind.pos] = NA
  phi.hat.pos = phi.hat[,p]
  phi.hat.pos[ind.neg] = NA
  base$val.neg[ind.cells] = phi.hat.neg
  base$val.pos[ind.cells] = phi.hat.pos
  base$val[ind.cells] = phi.hat[,p]
  ind.sig = phi.sig[,p] != 0
  phi.plot.sig= rep(NA,dim(phi.hat)[1])
  phi.plot.sig[ind.sig] = 1
  base$val.sig[ind.cells] = phi.plot.sig
  p4 = as.matrix(rasterToPoints(mask(base$val.sig,m1)))
  #p4 = p4[!is.na(p4[,3]),]
  base = disaggregate(base,10)
  windows(width = 10)
  pdf(paste0(path,"RCWinter50 ",plot.names[p],'.pdf'))
  plot(m1,main = plot.names[p])
  p1 = mask(base$val.neg,m1)
  p2 = mask(base$val.pos,m1)
  p3 = mask(base$val,m1)
  plot(p1,add = TRUE,zlim = c(min(phi.hat),max(phi.hat)),col = brewer.pal(9,"Blues")[seq(from = 9, to = 1)],legend = FALSE,breaks = seq(from = min(phi.hat),to = 0, length.out = 10))
  plot(p2,add = TRUE,zlim = c(min(phi.hat),max(phi.hat)),col = brewer.pal(9,"Reds"),legend = FALSE,breaks = c(seq(from = 0, to = max(phi.hat),length.out = 9)))
  #plot(p1,add = TRUE,zlim = c(min(phi.hat),max(phi.hat)),col = rainbow(10,start = .4, end = .6,rev = TRUE),legend = FALSE)
  #plot(p2,add = TRUE,zlim = c(min(phi.hat),max(phi.hat)),col = rainbow(10,start = 0, end = .4,rev = TRUE),legend = FALSE)
  points(p4[,1],p4[,2],pch = 19)
  plot(m1,add = TRUE)
  dev.off()
  dev.off()
  #plot(m1)
  #plot(p3,add = TRUE)
  #plot(m1,add = TRUE)
  #dev.off()
}
windows(width = 3)
pdf(paste0(path,"RWinter50 Legend.pdf"))
#image(z = matrix(c(seq(from = min(phi.hat),to = 0, length.out = 10),seq(from = 0, to = max(phi.hat),length.out = 10)),1,20), col = rainbow(20,start = 0,end = 0.6,rev = TRUE ),axes = FALSE)
image(z = t(matrix(c(seq(from = min(phi.hat),to = 0, length.out = 9),seq(from = 0, to = max(phi.hat),length.out = 9)),1,18)), col = c(brewer.pal(9,"Blues")[seq(from = 9, to = 1)],brewer.pal(9,"Reds")),axes = FALSE,breaks = c(seq(from = min(phi.hat),to = 0, length.out = 10),seq(from = 0, to = max(phi.hat),length.out = 9)))
axis(labels = c(round(min(phi.hat),2),0,round(max(phi.hat),2)),at = c(0,0.5,1),side = 1)
dev.off()
dev.off()

###Make the Significance Maps
for(p in 1:P){
  base = raster(xmn = min(x.grid),xmx = max(x.grid),ymn = min(y.grid),ymx = max(y.grid),nrow = g.s, ncol = g.s)
  base$val <- NA
  base$val.neg <- NA
  base$val.pos <- NA
  ind.cells = extract(base,y = cbind(x.grid,y.grid),cellnumbers = TRUE)[,1]
  ind.neg = phi.sig[,p] == -1
  ind.pos = phi.sig[,p] == 1
  phi.hat.neg = rep(NA,dim(phi.hat)[1])
  phi.hat.neg[ind.neg] = -1
  phi.hat.pos = rep(NA,dim(phi.hat)[1])
  phi.hat.pos[ind.pos] = 1
  base$val.neg[ind.cells] = phi.hat.neg
  base$val.pos[ind.cells] = phi.hat.pos
  base$val[ind.cells] = phi.sig[,p]
  base = disaggregate(base,10)
  windows(width = 10)
  pdf(paste0(path,"Sig Winter50 ",plot.names[p],'.pdf'))
  plot(m1,main = plot.names[p])
  p1 = mask(base$val.neg,m1)
  p2 = mask(base$val.pos,m1)
  p3 = mask(base$val,m1)
  #plot(p3,add = TRUE,zlim = c(-1,1),col = rainbow(3,start = 0,end = 0.6))
  plot(p1,add = TRUE,zlim = c(-1,1),col = rainbow(10,start = .4, end = .6,rev = TRUE),legend = FALSE)
  plot(p2,add = TRUE,zlim = c(-1,1),col = rainbow(10,start = 0, end = .4,rev = TRUE),legend = FALSE)
  plot(m1,add = TRUE)
  dev.off()
  dev.off()
  
  res = list(phi.array, tau.array, upsilon.array, psi.array, gamma.array)
  return(res)
}


