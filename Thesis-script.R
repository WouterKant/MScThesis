############ FUNCTION DEFINITIONS #############

### SPLIT SAFE MWU
{
library(rootSolve)

ucalc = function(a,b){
  na = length(a)
  nb = length(b)
  a = setNames(a,rep('a',na))
  b = setNames(b,rep('b',nb))
  ties = isTRUE(intersect(a,b) != 0)
  if(isTRUE(ties)){
    for(i in unique(r)){
      ti = sum(r==i)
      tiescorrection = tiescorrection + (ti^3 - ti)/((na+nb)*(na+nb-1))
    }
  }
  r = rank(c(a,b),ties.method = "average")
  u = sum(r[names(r) == 'a']) - na*(na+1)/2
  return(u)
}

lcalc = function(delta,na,nb){
  mean = (delta+0.5)*na*nb
  expect = function(l){
    denomH1 = dwilcox(0:(na*nb),na,nb)*exp(l*(0:(na*nb)/(na*nb)-0.5))
    sum((denomH1*0:(na*nb))/sum(denomH1)) - mean
  }
  return(multiroot(expect,0)$root)
}

SplitMWUDelta = function(a,b,prior = 'uniform'){
  evalgrid = seq(-0.5,0.5,0.1)
  if(prior == 'uniform'){
    Wd = rep(1/length(evalgrid),length(evalgrid))
  }
  if(prior == 'normal'){
    Wd = dnorm(evalgrid,0,0.5)
  }
  if(prior == 'invnorm'){
    Wd = 1/dnorm(evalgrid,0,0.1)
  }
  Wd = Wd/sum(Wd)
  Evals = numeric(length(a))
  stopt = length(a)
  check = FALSE
  tWd = numeric(length(Wd))
  for(i in 1:length(a)){
    u = ucalc(a[i],b[i])
    na = nb = 1
    for(j in 1:length(Wd)){
      l = lcalc(evalgrid[j],na,nb)
      tWd[j] = Wd[j]*exp(l*(u/(na*nb)-0.5))/sum(dwilcox(0:(na*nb),na,nb)*exp(l*(0:(na*nb)/(na*nb)-0.5)))
    }
    Evals[i] = sum(tWd)
    Wd = tWd/sum(tWd)
    if(prod(Evals[1:i])>20 && check == FALSE){
      check = TRUE
      stopt = i
    }
  }
  return(list(Evals,stopt,prod(Evals)))
}


SplitMWUDeltaStop = function(a,b,prior = 'uniform'){
  evalgrid = seq(-0.5,0.5,0.1)
  if(prior == 'uniform'){
    Wd = rep(1/length(evalgrid),length(evalgrid))
  }
  if(prior == 'normal'){
    Wd = dnorm(evalgrid,0,0.5)
  }
  if(prior == 'invnorm'){
    Wd = 1/dnorm(evalgrid,0,0.5)
  }
  Wd = Wd/sum(Wd)
  Evals = numeric(length(a))
  stopt = length(a)
  check = FALSE
  tWd = numeric(length(Wd))
  i = 0
  while(i < length(a) && prod(Evals[1:i])<20){
    i = i+1
    u = ucalc(a[i],b[i])
    na = nb = 1
    for(j in 1:length(Wd)){
      l = lcalc(evalgrid[j],na,nb)
      tWd[j] = Wd[j]*exp(l*(u/(na*nb)-0.5))/sum(dwilcox(0:(na*nb),na,nb)*exp(l*(0:(na*nb)/(na*nb)-0.5)))
    }
    Evals[i] = sum(tWd)
    Wd = tWd/sum(tWd)
    if(prod(Evals[1:i])>20 && check == FALSE){
      check = TRUE
      stopt = i
    }
  }
  return(list(Evals,stopt,prod(Evals[1:i])))
}

}

### SEQUENTIAL RANK TEST
{
MurielsPrior = function(data,prior = 'divided'){
  N = nrow(data)
  E = numeric(N-1)
  E[1] = 1
  #theta = 1/thetagrid
  #Wtheta = Wtheta/sum(Wtheta)
  a = seq(-100,100)
  if(prior == 'divided'){
    wgrid = (1+0.05)^a
    thetagrid = numeric(length(wgrid)-1)
    for(i in 1:(length(wgrid)-1)){
      thetagrid[i] = (wgrid[i+1]+wgrid[i])/2
    }
    Wtheta = rep(1/length(thetagrid),length(thetagrid))
  }
  if(prior == 'normal'){
    thetagrid = seq(0.01,3,0.01)
    Wtheta = dnorm(thetagrid,mean = 1)
  }
  if(prior == 'uniform'){
    thetagrid = seq(0.01,3,0.01)
    Wtheta = rep(1/length(thetagrid),length(thetagrid))
  }
  Wtheta = Wtheta/sum(Wtheta)
  for(n in 1:(N-1)){
    sortdata = matrix(matrix(data[1:n,],ncol = 2)[order(matrix(data[1:n,],ncol = 2)[,1]),],ncol =2)
    G = matrix(data = 0, nrow = n,ncol = 2)
    for(j in 1:n){
      G[j,1] = sum(sortdata[j:n,2] == 1)
      G[j,2] = sum(sortdata[j:n,2] == 2)
    }
    k = rank(c(sortdata[,1],data[n+1,1]))[n+1]
    ind = c(data[n+1,2]==1,data[n+1,2]==2)
    if(k == n+1){
      Ttheta = sapply(thetagrid,function(theta){
        thetas = c(1,theta)
        return(prod(G[1:(k-1),]%*%thetas/(c(ind%*%thetas)+ G[1:(k-1),]%*%thetas))* (n+1))
      })
    }else{
      Ttheta = sapply(thetagrid,  function(theta){
        thetas = c(1,theta)
        return(c(ind%*%thetas)/(c(ind%*%thetas)+G[k,]%*%thetas) * prod(G[1:(k-1),]%*%thetas/(c(ind%*%thetas)+ G[1:(k-1),]%*%thetas))* (n+1))
      })
    }
    Ttheta
    Twtheta = Wtheta * Ttheta
    E[n] = sum(Twtheta)
    Wtheta = Twtheta / sum(Twtheta)
  }
  return(list(c(1,E),prod(E), Wtheta))
}

MurielsPriorStop = function(data,prior = 'divided'){
  N = nrow(data)
  E = numeric(N-1)
  E[1] = 1
  #theta = 1/thetagrid
  #Wtheta = Wtheta/sum(Wtheta)
  a = seq(-100,100)
  if(prior == 'divided'){
    wgrid = (1+0.05)^a
    thetagrid = numeric(length(wgrid)-1)
    for(i in 1:(length(wgrid)-1)){
      thetagrid[i] = (wgrid[i+1]+wgrid[i])/2
    }
    Wtheta = rep(1/length(thetagrid),length(thetagrid))
  }
  if(prior == 'normal'){
    thetagrid = seq(0.01,3,0.01)
    Wtheta = dnorm(thetagrid,mean = 1)
  }
  if(prior == 'uniform'){
    thetagrid = seq(0.01,3,0.01)
    Wtheta = rep(1/length(thetagrid),length(thetagrid))
  }
  Wtheta = Wtheta/sum(Wtheta)
  n = 0
  while(n <(N-1) && prod(E[1:n]) < 20){
    n = n+1
    sortdata = matrix(matrix(data[1:n,],ncol = 2)[order(matrix(data[1:n,],ncol = 2)[,1]),],ncol =2)
    G = matrix(data = 0, nrow = n,ncol = 2)
    for(j in 1:n){
      G[j,1] = sum(sortdata[j:n,2] == 1)
      G[j,2] = sum(sortdata[j:n,2] == 2)
    }
    k = rank(c(sortdata[,1],data[n+1,1]))[n+1]
    ind = c(data[n+1,2]==1,data[n+1,2]==2)
    if(k == n+1){
      Ttheta = sapply(thetagrid,function(theta){
        thetas = c(1,theta)
        return(prod(G[1:(k-1),]%*%thetas/(c(ind%*%thetas)+ G[1:(k-1),]%*%thetas))* (n+1))
      })
    }else{
      Ttheta = sapply(thetagrid,  function(theta){
        thetas = c(1,theta)
        return(c(ind%*%thetas)/(c(ind%*%thetas)+G[k,]%*%thetas) * prod(G[1:(k-1),]%*%thetas/(c(ind%*%thetas)+ G[1:(k-1),]%*%thetas))* (n+1))
      })
    }
    Ttheta
    Twtheta = Wtheta * Ttheta
    E[n] = sum(Twtheta)
    Wtheta = Twtheta / sum(Twtheta)
  }
  return(list(c(1,E),prod(E), n,Wtheta))
}
}

#HELPER FUNCTION FOR CALCULATION [1]
logcosh = function(x){
  s = sign(x) * x
  p = exp(-2 * s)
  return(s + log1p(p) - log(2))
}

# LOGISTIC ALTERNATIVE GENERATIVE DISTRIBUTION
{
#data generator
rlogalt = function(n, delta = 0){
  unif = runif(n)
  return(-log(1/unif-1)+delta)
}

#check if it works
{grid = seq(-10,10)
  plot(1/(1+exp(-grid)))
  hist(rlogalt(1000,0))
  plot(sort(rlogalt(1000,0)),1:1000)
}
}

# EFRON-DE LA PEÑA E-VALUE FUNCTIONS
{
Ecorr= function(lambda,z){
  if(z>=0){
    ind = lambda>=0
    return((-ind+1)*exp(abs(lambda)*-z - 0.5*lambda^2*(-z)^2)  +  ind*(2 - exp(abs(lambda)*(-z) - 0.5*lambda^2*(-z)^2)))
  }
  if(z<0){
    ind = lambda>=0
    return((ind)*exp(abs(lambda)*z - 0.5*lambda^2*(z)^2)  +  (-ind+1)*(2 - exp(abs(lambda)*(z) - 0.5*lambda^2*(z)^2)))
  }
}

Ecosh = function(lambda,z){
  return(exp(lambda*z - logcosh(lambda*z)))
}

WEcosh = function(lambda,z,W){
  return((exp(lambda*z - logcosh(lambda*z)))*(1-W)+W)
}

Epena = function(lambda,z){
  return(exp(lambda*z - 0.5*lambda^2*z^2))
}

Eram = function(lambda,z){
  if(z<0){
    ind = lambda>=0
    return((-ind+1)*exp(abs(lambda)*-z - 0.5*lambda^2*(-z)^2)  +  ind*(2 - exp(abs(lambda)*(-z) - 0.5*lambda^2*(-z)^2)))
  }
  if(z>=0){
    ind = lambda>=0
    return((ind)*exp(abs(lambda)*z - 0.5*lambda^2*(z)^2)  +  (-ind+1)*(2 - exp(abs(lambda)*(z) - 0.5*lambda^2*(z)^2)))
  }
}

WEram = function(lambda,z,W = 0){
  if(z<0){
    ind = lambda>=0
    return( ((-ind+1)*exp(abs(lambda)*-z - 0.5*lambda^2*(-z)^2)  +  ind*(2 - exp(abs(lambda)*(-z) - 0.5*lambda^2*(-z)^2)))*(1-W)+W)
  }
  if(z>=0){
    ind = lambda>=0
    return(  ((ind)*exp(abs(lambda)*z - 0.5*lambda^2*(z)^2)  +  (-ind+1)*(2 - exp(abs(lambda)*(z) - 0.5*lambda^2*(z)^2)))*(1-W)+W)
  }
}

}

# EFRON-DE LA PEÑA EXPERIMENT FUNCTIONS WITH OPTIONAL STOPPING
{
Efronspowertest = function(z, prior = 'unif',shapeinvg = 1,rate = 1,cascale = 1,grid = 'sup'){
  lgrid = seq(0.1,100,0.1)
  if(prior == 'cauchy'){
    Wl = dcauchy(lgrid, scale = cascale)
  }
  else if(prior == 'invgamma'){
    Wl = dinvgamma(lgrid,shapeinvg,rate)
  }
  else{Wl = rep(1,length(lgrid))}
  Wl = Wl/sum(Wl)
  Sr = 1
  i = 0
  while(i< length(z) && Sr<20){
    i = i+1
    tWl = Wl * Epena(lgrid,z[i])
    Sr = Sr * sum(tWl)
    Wl = tWl/sum(tWl)
  }
  return(c(i,Sr>20,Sr))
}

Efronspowertestexp = function(z, prior = 'unif',shapeinvg = 1,rate = 1,cascale = 1,grid = 'exp', type = 'NewE', twosided = FALSE){
  if(grid == 'exp'){
    a = seq(-2000,2000)
    lgrid = (1+0.01)^a
    evalgrid = numeric(length(lgrid)-1)
    cgrid = numeric(length(lgrid)-1)
    for(i in 1:(length(lgrid)-1)){
      evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
      cgrid[i] = lgrid[i+1]-lgrid[i]
    }
  }
  if(prior == 'cauchy'){
    Wl = dcauchy(evalgrid, scale = 1)*cgrid
  }
  else if(prior == 'invgamma'){
    Wl = dinvgamma(evalgrid,shapeinvg,rate)*cgrid
  }
  else{Wl = rep(1,length(evalgrid))*cgrid}
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'CorrectedE'){
    while(i< length(z) && Sr<20){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z) && Sr<20){
      i = i+1
      tWl = Wl * Ecosh(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z) && Sr<20){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

Efronspowertestmaxl = function(z, type = 'NewE',twosided = FALSE,W = 0){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'CorrectedE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z) & Sr<20){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],W)
      Erec[i] = WEcosh(evalgrid[which.max(Wl)],z[i],W)
      if(i == 1){
        Erec[1] = 1
      }
      Sr = Sr * Erec[i]
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
      if(i == 1){
        Sr = 1
        Erec[1] = 1
      }
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

RamdasTest = function(z,sigmethod = 'second',Eval ='Efron',W = 0){
  sigma = 0
  a = 0
  E = numeric(length(z))
  i = 0
  while(i<length(z) & prod(E[1:i])<20){
    i = i+1
    if(sigmethod == 'first'){
      a = 1/sigma
    }
    if(sigmethod == 'second'){
      a = mean(z[1:(i-1)])/(mean(z[1:(i-1)])^2 + sigma)
    }
    if(Eval == 'Efron'){
      E[i] = Epena(a,z[i])
    }
    if(Eval == 'CorrectedE'){
      E[i] = Ecorr(a,z[i])
    }
    if(Eval == 'NewE'){
      E[i] = WEcosh(a,z[i],W)
    }
    if(Eval == 'Ram'){
      E[i] = Eram(a,z[i])
    }
    if(i == 1 | i == 2){
      E[i] = 1
    }
    sigma = sum((z[1:i]-mean(z[1:i]))^2)/i
  }
  return(list(E,prod(E[1:i]),i,prod(E[1:i])>20))
}

RamdasTestV = function(z,sigmethod = 'second',Eval ='Efron',startV = 0){
  sigma = 0
  a = 0
  E = numeric(length(z))
  i = 0
  Vind = seq(0,1,0.01)
  V = c(startV, rep(0.25/100,100))
  while(i<length(z) & prod(E[1:i])<20){
    i = i+1
    if(sigmethod == 'first'){
      a = 1/sigma
    }
    if(sigmethod == 'second'){
      a = mean(z[1:(i-1)])/(mean(z[1:(i-1)])^2 + sigma)
    }
    if(Eval == 'Efron'){
      E[i] = Epena(a,z[i])
    }
    if(Eval == 'CorrectedE'){
      E[i] = Ecorr(a,z[i])
    }
    if(Eval == 'NewE'){
      E[i] = Ecosh(a,z[i],sum(V*Vind))
      V = V * ((1-Vind)*E[i]+Vind)
      V = V/sum(V)
    }
    if(i == 1 | i == 2){
      E[i] = 1
    }
    sigma = sum((z[1:i]-mean(z[1:i]))^2)/i
  }
  return(list(E,prod(E[1:i]),i,prod(E[1:i])>20))
}

EfronsScaleInvNorm = function(z, normsd = 1,grid = 'exp', type = 'CorrectedE',twosided = FALSE){
  evalgrid = seq(0.1,100,0.1)
  if(grid == 'exp'){
    a = seq(-1000,1000)
    lgrid = (1+0.01)^a
    evalgrid = numeric(length(lgrid)-1)
    cgrid = numeric(length(lgrid)-1)
    for(i in 1:(length(lgrid)-1)){
      evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
      cgrid[i] = lgrid[i+1]-lgrid[i]
    }
  }
  Wl = dnorm(evalgrid, sd = normsd)*cgrid
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  z = z/z[1]
  if(type == 'CorrectedE'){
    while(i< length(z) && (Sr<20|is.nan(Sr))){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z) && (Sr<20|is.nan(Sr))){
      i = i+1
      tWl = Wl * Ecosh(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z) && (Sr<20|is.nan(Sr))){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'Ram'){
    while(i< length(z) && (Sr<20|is.nan(Sr))){
      i = i+1
      tWl = Wl * Eram(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
      if(i == 1){
        Sr = 1
        Erec[1] = 1
      }
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec,prod(Erec[1:i])>20))
}

EfronsGExp = function(z, type = 'CorrectedE', twosided = FALSE){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'CorrectedE'){
    while(i< length(z) && Sr<20){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z) && Sr<20){
      i = i+1
      tWl = Wl * Ecosh(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z)&& Sr<20){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
      if(is.nan(Sr)){
        Sr = 0
      }
    }
  }
  if(type == 'Ram'){
    while(i< length(z)&& Sr<20){
      i = i+1
      tWl = Wl * SRamCorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
      if(i == 1){
        Sr = 1
        Erec[1] = 1
      }
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec,prod(Erec[1:i])>20))
}

EfronsGExpW = function(z, type = 'CorrectedE', twosided = FALSE,W){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'NewE'){
    while(i< length(z) && Sr<20){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],W)
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'Ram'){
    while(i< length(z) && (Sr<20|is.nan(Sr))){
      i = i+1
      tWl = Wl * WEram(evalgrid,z[i],W)
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

EfronsGExpV = function(z, type = 'NewE', twosided = FALSE,startV = 0.5){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  n = length(z)
  Vind = seq(0,1,0.01)
  V = c(startV, rep(0.25/100,100))
  if(type == 'CorrectedE'){
    while(i< length(z) && Sr<20){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z) && Sr<20){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],sum(V*Vind))
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
      V = V * ((1-Vind)*Sr+Vind)
      V = V/sum(V)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z) && (Sr<20|is.nan(Sr))){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec,V))
}

EfronsGExpSigma = function(z, type = 'NewE', twosided = FALSE){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'NewE'){
    while(i< length(z) && Sr<20){
      i = i+1
      tWl = Wl * Snewsigma(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z)&& Sr<20){
      i = i+1
      tWl = Wl * Ssigma(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
      if(is.nan(Sr)){
        Sr = 0
      }
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec,prod(Erec[1:i])>20))
}

Ebernoulli = function(z){
  y = (z>=0)*2-1
  theta = 0.5
  i = 0
  Erec = numeric(length(z))
  while(i< length(z) && prod(Erec[1:i])<20){
    i = i+1
    if(y[i] == 1){
      Erec[i] = theta/0.5
    }else{
      Erec[i] = (1-theta)/0.5
    }
    theta = (sum(y[1:i]>0)+0.5)/(i+1)
  }
  return(list(i,prod(Erec[1:i]),Erec,prod(Erec[1:i])>20))
}
}

# EFRON-DE LA PEÑA EXPERIMENT FUNCTIONS WITHOUT OPTIONAL STOPPING FOR EXPERIMENT, CALLED LOG BECAUSE THEY ARE TO BE PLOTTED ON LOGARITHMIC SCALE
{
Efronspowertestexplog = function(z, prior = 'unif',shapeinvg = 1,rate = 1,cascale = 1,grid = 'exp', type = 'CorrectedE',twosided = FALSE){
  evalgrid = seq(0.1,100,0.1)
  if(grid == 'exp'){
    a = seq(-1000,1000)
    lgrid = (1+0.01)^a
    evalgrid = numeric(length(lgrid)-1)
    cgrid = numeric(length(lgrid)-1)
    for(i in 1:(length(lgrid)-1)){
      evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
      cgrid[i] = lgrid[i+1]-lgrid[i]
    }
  }
  if(prior == 'cauchy'){
    Wl = dcauchy(evalgrid, scale = 1)*cgrid
  }
  else if(prior == 'invgamma'){
    Wl = dinvgamma(evalgrid,shapeinvg,rate)*cgrid
  }
  else{Wl = rep(1,length(evalgrid))*cgrid}
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  i = 0
  if(type == 'CorrectedE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      #      if()
      tWl = Wl * Ecosh(evalgrid,z[i])
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl))
}

Efronspowertestmaxllog = function(z, type = 'NewE',twosided = FALSE,h = 0){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'CorrectedE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],h)
      Erec[i] = WEcosh(evalgrid[which.max(Wl)],z[i],h)
      if(i == 1){
        Erec[1] = 1
      }
      Sr = Sr * Erec[i]
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
      if(i == 1){
        Sr = 1
        Erec[1] = 1
      }
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

RamdasTestlog = function(z,sigmethod = 'second',Eval ='Efron'){
  sigma = 0
  a = 0
  E = numeric(length(z))
  i = 0
  while(i<length(z)){
    i = i+1
    a = mean(z[1:(i-1)])/(mean(z[1:(i-1)])^2 + sigma)
    if(Eval == 'Efron'){
      E[i] = Epena(a,z[i])
    }
    if(Eval == 'CorrectedE'){
      E[i] = Ecorr(a,z[i])
    }
    if(Eval == 'NewE'){
      E[i] = Ecosh(a,z[i])
    }
    if(Eval == 'Ram'){
      E[i] = Eram(a,z[i])
    }
    if(i == 1 | i == 2){
      E[i] = 1
    }
    sigma = sum((z[1:i]-mean(z[1:i]))^2)/i
  }
  return(list(E,prod(E[1:i]),i))
}

RamdasTestWlog = function(z,sigmethod = 'second',Eval ='Efron',W){
  sigma = 0
  a = 0
  E = numeric(length(z))
  i = 0
  while(i<length(z)){
    i = i+1
    a = mean(z[1:(i-1)])/(mean(z[1:(i-1)])^2 + sigma)
    if(Eval == 'Efron'){
      E[i] = Epena(a,z[i])
    }
    if(Eval == 'CorrectedE'){
      E[i] = Ecorr(a,z[i])
    }
    if(Eval == 'NewE'){
      E[i] = WEcosh(a,z[i],W)
    }
    if(Eval == 'Ram'){
      E[i] = WEram(a,z[i],W)
    }
    if(i == 1 | i == 2){
      E[i] = 1
    }
    sigma = sum((z[1:i]-mean(z[1:i]))^2)/i
  }
  return(list(E,prod(E[1:i]),i))
}

RamdasTestVlog = function(z,sigmethod = 'second',Eval ='Efron',startV = 0){
  sigma = 0
  a = 0
  E = numeric(length(z))
  i = 0
  Vind = seq(0,1,0.01)
  V = c(startV, rep(0.25/100,100))
  while(i<length(z)){
    i = i+1
    if(sigmethod == 'first'){
      a = 1/sigma
    }
    if(sigmethod == 'second'){
      a = mean(z[1:(i-1)])/(mean(z[1:(i-1)])^2 + sigma)
    }
    if(Eval == 'Efron'){
      E[i] = Epena(a,z[i])
    }
    if(Eval == 'CorrectedE'){
      E[i] = Ecorr(a,z[i])
    }
    if(Eval == 'NewE'){
      E[i] = WEcosh(a,z[i],sum(V*Vind))
      V = V * ((1-Vind)*E[i]+Vind)
      V = V/sum(V)
    }
    if(i == 1 | i == 2){
      E[i] = 1
    }
    sigma = sum((z[1:i]-mean(z[1:i]))^2)/i
  }
  return(list(E,prod(E[1:i]),i,prod(E[1:i])>20))
}

EfronsScaleInvNormlog = function(z, normsd = 1,grid = 'exp', type = 'CorrectedE',twosided = FALSE){
  evalgrid = seq(0.1,100,0.1)
  if(grid == 'exp'){
    a = seq(-1000,1000)
    lgrid = (1+0.01)^a
    evalgrid = numeric(length(lgrid)-1)
    cgrid = numeric(length(lgrid)-1)
    for(i in 1:(length(lgrid)-1)){
      evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
      cgrid[i] = lgrid[i+1]-lgrid[i]
    }
  }
  Wl = dnorm(evalgrid, sd = normsd)*cgrid
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  z = z/z[1]
  if(type == 'CorrectedE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Ecosh(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'Ram'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Eram(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

EfronsScaleInvNormWlog = function(z, normsd = 1,grid = 'exp', type = 'CorrectedE',twosided = FALSE,W){
  evalgrid = seq(0.1,100,0.1)
  if(grid == 'exp'){
    a = seq(-1000,1000)
    lgrid = (1+0.01)^a
    evalgrid = numeric(length(lgrid)-1)
    cgrid = numeric(length(lgrid)-1)
    for(i in 1:(length(lgrid)-1)){
      evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
      cgrid[i] = lgrid[i+1]-lgrid[i]
    }
  }
  Wl = dnorm(evalgrid, sd = normsd)*cgrid
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  z = z/z[1]
  if(type == 'CorrectedE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],W)
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'Ram'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEram(evalgrid,z[i],W)
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

EfronsGExplog = function(z, type = 'CorrectedE', twosided = FALSE){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'CorrectedE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Ecorr(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Ecosh(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'OrigE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Epena(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
      if(i == 1){
        Sr = 1
        Erec[1] = 1
      }
    }
  }
  if(type == 'Ram'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * Eram(evalgrid,z[i])
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
      if(i == 1){
        Sr = 1
        Erec[1] = 1
      }
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

EfronsGExpWlog = function(z, type = 'NewE', twosided = FALSE,W){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],W)
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'Ram'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEram(evalgrid,z[i],W)
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = tWl/sum(tWl)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

EfronsGExpWlognolik = function(z, type = 'NewE', twosided = FALSE,W){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],W)
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      tWl2 = Wl * Ecosh(evalgrid,z[i])
      Wl = tWl2/sum(tWl2)
    }
  }
  if(type == 'Ram'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEram(evalgrid,z[i],W)
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      tWl2 = Wl * Eram(evalgrid,z[i])
      Wl = tWl2/sum(tWl2)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

EfronsGExpWlog3 = function(z, type = 'NewE', twosided = FALSE,W){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],W)
      tWl3 = Wl * Ecosh(evalgrid,z[i])
      Erec[i] = sum(tWl3)
      Sr = Sr * sum(tWl3)
      Wl = tWl/sum(tWl)
    }
  }
  if(type == 'Ram'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEram(evalgrid,z[i],W)
      tWl3 = Wl * Eram(evalgrid,z[i])
      Erec[i] = sum(tWl3)
      Sr = Sr * sum(tWl3)
      Wl = tWl/sum(tWl)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec))
}

EfronsGExpVlog = function(z, type = 'CorrectedE', twosided = FALSE,startV = 0.5){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  n = length(z)
  Vind = seq(0,1,0.01)
  V = c(startV, rep(0.25/100,100))
  V = V/sum(V)
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],sum(V*Vind))
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = (tWl)/sum(tWl)
      V = V * ((1-Vind)*Sr+Vind)
      V = V/sum(V)
    }
  }
  if(type == 'Ram'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEram(evalgrid,z[i],sum(V*Vind))
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = (tWl)/sum(tWl)
      V = V * ((1-Vind)*Sr+Vind)
      V = V/sum(V)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec,sum(V*Vind),V ))
}

EfronsGExpVlog2point = function(z, type = 'CorrectedE', twosided = FALSE){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  n = length(z)
  Vind = c(0,0.2)
  V = c(0.5,0.5)
  V = V/sum(V)
  if(type == 'NewE'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEcosh(evalgrid,z[i],sum(V*Vind))
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = (tWl)/sum(tWl)
      V = V * ((1-Vind)*Sr+Vind)
      V = V/sum(V)
    }
  }
  if(type == 'Ram'){
    while(i< length(z)){
      i = i+1
      tWl = Wl * WEram(evalgrid,z[i],sum(V*Vind))
      Erec[i] = sum(tWl)
      Sr = Sr * sum(tWl)
      Wl = (tWl)/sum(tWl)
      V = V * ((1-Vind)*Sr+Vind)
      V = V/sum(V)
    }
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec,sum(V*Vind),V ))
}

EfronsGExpSigmalog = function(z, twosided = FALSE){
  a = seq(-1000,1000)
  lgrid = (1+0.01)^a
  evalgrid = numeric(length(lgrid)-1)
  for(i in 1:(length(lgrid)-1)){
    evalgrid[i] = (lgrid[i+1]+lgrid[i])/2
  }
  Wl = rep(1/length(evalgrid),length(evalgrid))
  if(twosided == TRUE){
    Wl = c(rev(Wl),Wl)
    evalgrid = c(-rev(evalgrid),evalgrid)
  }
  Wl = Wl/sum(Wl)
  Sr = 1
  Erec = numeric(length(z))
  i = 0
  while(i< length(z)){
    i = i+1
    tWl = Wl * Ecoshsigma(evalgrid,z[i])
    Erec[i] = sum(tWl)
    Sr = Sr * sum(tWl)
    Wl = tWl/sum(tWl)
  }
  return(list(i,Sr>20,evalgrid[which.max(Wl)],Sr,Wl,Erec,prod(Erec[1:i])>20))
}

Ebernoullilog = function(z){
  y = (z>=0)*2-1
  theta = 0.5
  i = 0
  Erec = numeric(length(z))
  while(i< length(z)){
    i = i+1
    if(y[i] == 1){
      Erec[i] = theta/0.5
    }else{
      Erec[i] = (1-theta)/0.5
    }
    theta = (sum(y[1:i]>0)+0.5)/(i+1)
  }
  return(list(i,prod(Erec[1:i]),Erec,prod(Erec[1:i])>20))
}

}

# STOPPING TIMES MWU UNDER DIFFERENT DISTRIBUTIONS
{
{
  # MWU TEST & SMWU TEST
  thetagrid = seq(0.2,3,0.1)
  B = 1000
  ngrid = seq(1,1000,1)
  reg = numeric(B)
  alpha = 0.05
  regmwu = numeric(B)
  regmwurecrec = numeric(length(thetagrid))
  safmwu = numeric(B)
  safmwurecrec = numeric(length(thetagrid))
  for(i in 1:length(thetagrid)){
    rm = TRUE
    sm = TRUE
    q = 1
    while(q < max(ngrid) && rm == T){
      for(j in 1:B){
        x = rnorm(ngrid[q],thetagrid[i])
        y = rnorm(ngrid[q])
        if(rm == T){regmwu[j] = wilcox.test(x,y, paired = FALSE)$p.value}
      }
      if(rm == T){regt = mean(regmwu<alpha)
      if(regt >=0.8){
        regmwurecrec[i] = ngrid[q]
        rm = F
      }
      }
      print(ngrid[q])
      q = q+1
    }
    print(c('HEYOO', thetagrid[i]))
  }
}
regmwurecrec

#T DISTRIBUTION
{
  # MWU TEST & SMWU TEST
  thetagrid = seq(0.2,3,0.1)
  B = 1000
  ngrid = seq(1,1000,1)
  reg = numeric(B)
  alpha = 0.05
  regmwu = numeric(B)
  regmwurecrecT = numeric(length(thetagrid))
  safmwu = numeric(B)
  safmwurecrec = numeric(length(thetagrid))
  for(i in 1:length(thetagrid)){
    rm = TRUE
    sm = TRUE
    q = 1
    while(q < max(ngrid) && rm == T){
      for(j in 1:B){
        x = rt(ngrid[q],3,thetagrid[i])
        y = rt(ngrid[q],3)
        if(rm == T){regmwu[j] = wilcox.test(x,y, paired = FALSE)$p.value}
      }
      if(rm == T){regt = mean(regmwu<alpha)
      if(regt >=0.8){
        regmwurecrecT[i] = ngrid[q]
        rm = F
      }
      }
      print(ngrid[q])
      q = q+1
    }
    print(c('HEYOO', thetagrid[i]))
  }
}
regmwurecrecT

#CAUCHY DISTRIUTION
{
  # MWU TEST & SMWU TEST
  thetagrid = seq(0.2,3,0.1)
  B = 1000
  ngrid = seq(1,1000,1)
  reg = numeric(B)
  alpha = 0.05
  regmwu = numeric(B)
  regmwurecrecCauch = numeric(length(thetagrid))
  safmwu = numeric(B)
  safmwurecrec = numeric(length(thetagrid))
  for(i in 1:length(thetagrid)){
    rm = TRUE
    sm = TRUE
    q = 1
    while(q < max(ngrid) && rm == T){
      for(j in 1:B){
        x = rt(ngrid[q],3,thetagrid[i])
        y = rt(ngrid[q],3)
        if(rm == T){regmwu[j] = wilcox.test(x,y, paired = FALSE)$p.value}
      }
      if(rm == T){regt = mean(regmwu<alpha)
      if(regt >=0.8){
        regmwurecrecCauch[i] = ngrid[q]
        rm = F
      }
      }
      print(ngrid[q])
      q = q+1
    }
    print(c('HEYOO', thetagrid[i]))
  }
}
regmwurecrecCauch

#LOGALT DISTRIBUTION
{
  # MWU TEST & SMWU TEST
  thetagrid = seq(0.2,3,0.1)
  B = 1000
  ngrid = seq(1,1000,1)
  reg = numeric(B)
  alpha = 0.05
  regmwu = numeric(B)
  regmwurecrecLogAlt = numeric(length(thetagrid))
  safmwu = numeric(B)
  safmwurecrec = numeric(length(thetagrid))
  for(i in 1:length(thetagrid)){
    rm = TRUE
    sm = TRUE
    q = 1
    while(q < max(ngrid) && rm == T){
      for(j in 1:B){
        x = rt(ngrid[q],3,thetagrid[i])
        y = rt(ngrid[q],3)
        if(rm == T){regmwu[j] = wilcox.test(x,y, paired = FALSE)$p.value}
      }
      if(rm == T){regt = mean(regmwu<alpha)
      if(regt >=0.8){
        regmwurecrecLogALt[i] = ngrid[q]
        rm = F
      }
      }
      print(ngrid[q])
      q = q+1
    }
    print(c('HEYOO', thetagrid[i]))
  }
}
regmwurecrecLogAlt
}

########### SECTION 1 #############

# FIGURE 1.1
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    Bern = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rnorm(N,1)
    Bernt = Ebernoullilog(z)[[3]]
    for(n in 1:N){
      Bern[n,b] = prod(Bernt[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(Bern)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Normal(mean = 1, variance = 1)',lwd = 2)
    #legend(0,20,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 10)
  }
}

# FIGURE 1.2
{
  {
    thetagrid = seq(0.20,2.5,0.10)
    library(invgamma)
    library(rootSolve)
    B = 1000
    mBern = numeric(length(thetagrid))
    wBern = numeric(length(thetagrid))
    Bernrec = numeric(B)
  }
  
  for(j in 1:length(thetagrid)){
    for(b in 1:B){
      a = rnorm(100000,thetagrid[j])
      c = rnorm(100000)
      z = a-c
      Bernrec[b] = Ebernoulli(z)[[1]]
      print(b)
    }
    
    wBern[j] = quantile(Bernrec,0.8)
    Bernrec[Bernrec>quantile(Bernrec,0.8)] = quantile(Bernrec,0.8)
    mBern[j] = mean(Bernrec)
    
    print(thetagrid[j])
  }
  {
    plot(thetagrid, mBern/regmwurecrec[1:24],xlim = c(0.25,2.5),ylim = c(1,4),col = 'black',xlab = 'Delta',ylab = 'Comparative average stopping time',type = 'l',lwd = 2,main = 'X ~ Normal(variance = 1)')
    lines(thetagrid,wBern/regmwurecrec[1:24],xlim = c(0.25,2.5),ylim = c(1,4),col = 'blue',lwd = 2)
    legend(1.5,4,c('E_Bern (Average)', 'E_Bern(Worst Case)'),text.col = c('black','blue'),text.width = 0.7)
  }
  
}

########### SECTION 2 #############
# TABLE 2.1
{
library(invgamma)
B = 1000
J= 5
cscale = numeric(J)
uscale = numeric(J)
gscal1 = numeric(J)
gscal2 = numeric(J)
gscal3 = numeric(J)
testscal1 = numeric(J)
testscal2 = numeric(J)

scalegrid = c(0.01,0.1,1,10,100)

for(j in 1:5){ # for each scale size
  cauchyrec = numeric(B)
  unifrec = numeric(B)
  invgrec1 = numeric(B)
  invgrec2 = numeric(B)
  invgrec3 = numeric(B)
  testrec1 = numeric(B)
  testrec2 = numeric(B)
  for(b in 1:B){
    z = rnorm(1000,1)*scalegrid[j]
    cauchyrec[b] = Efronspowertestexp(z,'cauchy')[[1]]
    unifrec[b] = Efronspowertestexp(z)[[1]]
    invgrec2[b] = Efronspowertestexp(z, prior = 'invgamma', 1)[[1]]
    invgrec3[b] = Efronspowertestexp(z, prior = 'invgamma', 10)[[1]]
    testrec1[b] = Efronspowertestexp(z,'cauchy',cascale = 0.5)[[1]]
    testrec2[b] = Efronspowertestexp(z,'cauchy',cascale = 2)[[1]]
    print(b)
  }
  cscale[j] = mean(cauchyrec)
  uscale[j] = mean(unifrec)
  gscal1[j] = mean(invgrec1)
  gscal2[j] = mean(invgrec2)
  gscal3[j] = mean(invgrec3)
  testscal1[j] = mean(testrec1)
  testscal2[j] = mean(testrec2)
}
cscale
uscale
gscal1
gscal2
gscal3
testscal1
testscal2
}

# FIGURE 2.1
{
plot(seq(-8,8,0.01),Epena(1,seq(-8,8,0.01)),type = 'l',lwd = 3,ylim = c(0,2),ylab ='E',xlab = 'z_i', main = 'E-variables for a singular data point with lambda = 1')

lines(seq(-8,8,0.01), sapply(seq(-8,8,0.01), function(z){
  lamba = 1
  if(z<0){
    ind = lambda>=0
    return((-ind+1)*exp(abs(lambda)*-z - 0.5*lambda^2*(-z)^2)  +  ind*(2 - exp(abs(lambda)*(-z) - 0.5*lambda^2*(-z)^2)))
  }
  if(z>=0){
    ind = lambda>=0
    return((ind)*exp(abs(lambda)*z - 0.5*lambda^2*(z)^2)  +  (-ind+1)*(2 - exp(abs(lambda)*(z) - 0.5*lambda^2*(z)^2)))
  }
}),type = 'l',lwd = 3,ylim = c(0,2),ylab ='E',xlab = 'z_i', main = 'E_Ram(lambda = 1)',col = 'red')

lines(seq(-8,8,0.01),Ecosh(1,seq(-8,8,0.01)),type = 'l',lwd = 3,ylim = c(0,2),ylab ='E',xlab = 'z_i', main = 'E_cosh(lambda = 1)',col = 'green')

lines(seq(-8,8,0.01), sapply(seq(-8,8,0.01), function(z){
  lambda = 1
  if(z>=0){
    ind = lambda>=0
    return((-ind+1)*exp(abs(lambda)*-z - 0.5*lambda^2*(-z)^2)  +  ind*(2 - exp(abs(lambda)*(-z) - 0.5*lambda^2*(-z)^2)))
  }
  if(z<0){
    ind = lambda>=0
    return((ind)*exp(abs(lambda)*z - 0.5*lambda^2*(z)^2)  +  (-ind+1)*(2 - exp(abs(lambda)*(z) - 0.5*lambda^2*(z)^2)))
  }
}),type = 'l',lwd = 2,ylim = c(0,2),ylab ='E_corr',xlab = 'z_i', main = 'E_corr(lambda = 1)',col = 'blue')

legend(-8,2,c('E_Peña','E_Ram','E_Corr','E_Cosh'), text.col = c('black','red', 'green', 'blue'))
max(Epena(1,seq(-5,5,0.01)))
}

############ SECTION 3 #############

#FIGURE 3.1 AND 3.3
{
{
  library(invgamma)
  library(rootSolve)
  B = 1000
  N = 50
  GrunwaldEfron = matrix(0,nrow = N,ncol = B)
  GrunwaldRam = matrix(0,nrow = N,ncol = B)
  GrunwaldCorr = matrix(0,nrow = N,ncol = B)
  GrunwaldCosh = matrix(0,nrow = N,ncol = B)
}
for(b in 1:B){
  z = rnorm(N,1)
  GrunwaldEfront = EfronsGExplog(z,type = 'OrigE',twosided = T)[[6]]
  GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
  GrunwaldCorrt = EfronsGExplog(z,type = 'CorrectedE',twosided = T)[[6]]
  GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
  for(n in 1:N){
    GrunwaldEfron[n,b] = prod(GrunwaldEfront[1:n])
    GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
    GrunwaldCorr[n,b] = prod(GrunwaldCorrt[1:n])
    GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
  }
  print(b)
}
{
  cols = rainbow(4)
  plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Normal(mean = 1, variance = 1)',ylim = c(0,15),lwd = 4)
  lines(rowSums(log(GrunwaldRam)/B),col = 'blue' ,lwd = 2)
  lines(rowSums(log(GrunwaldCorr)/B),col = 'red', lwd = 2)
  lines(rowSums(log(GrunwaldEfron)/B),col = 'green', lwd = 2)
  EfronBadNormGrun = recordPlot()
  legend(0,15,c('Gr.E_Peña','Gr.E_Ram','Gr.E_Corr','Gr.E_Cosh'),text.col = c('green','blue','red','black'),text.width = 12)
}
#with Ram method
{
  library(invgamma)
  library(rootSolve)
  B = 1000
  N = 50
  RamdasEfron = matrix(0,nrow = N,ncol = B)
  RamdasRam = matrix(0,nrow = N,ncol = B)
  RamdasCorr = matrix(0,nrow = N,ncol = B)
  RamdasCosh = matrix(0,nrow = N,ncol = B)
}
for(b in 1:B){
  z = rnorm(N,1)
  RamdasEfront = RamdasTestlog(z,Eval = 'Efron')[[1]]
  RamdasRamt = RamdasTestlog(z,Eval = 'Ram')[[1]]
  RamdasCorrt = RamdasTestlog(z,Eval = 'CorrectedE')[[1]]
  RamdasCosht = RamdasTestlog(z,Eval = 'NewE')[[1]]
  for(n in 1:N){
    RamdasEfron[n,b] = prod(RamdasEfront[1:n])
    RamdasRam[n,b] = prod(RamdasRamt[1:n])
    RamdasCorr[n,b] = prod(RamdasCorrt[1:n])
    RamdasCosh[n,b] = prod(RamdasCosht[1:n])
  }
  print(b)
}
{
  plot(rowSums(log(RamdasCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Normal(mean = 1, variance = 1)',lwd = 2)
  lines(rowSums(log(RamdasRam)/B),col = 'blue' ,lwd = 2)
  lines(rowSums(log(RamdasCorr)/B),col = 'red', lwd = 2)
  lines(rowSums(log(RamdasEfron)/B),col = 'green', lwd = 2)
  EfronbadNormRam = recordPlot()
  legend(0,14,c('Ram.E_Peña','Ram.E_Ram','Ram.E_Corr','Ram.E_Cosh'),text.col = c('green','blue','red','black'),text.width = 12)
}
}

#FIGURE 3.2
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    GrunwaldEfron = matrix(0,nrow = N,ncol = B)
    GrunwaldRam = matrix(0,nrow = N,ncol = B)
    GrunwaldCorr = matrix(0,nrow = N,ncol = B)
    GrunwaldCosh = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rcauchy(N,2)
    GrunwaldEfront = EfronsGExplog(z,type = 'OrigE',twosided = T)[[6]]
    GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
    GrunwaldCorrt = EfronsGExplog(z,type = 'CorrectedE',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
    for(n in 1:N){
      GrunwaldEfron[n,b] = prod(GrunwaldEfront[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCorr[n,b] = prod(GrunwaldCorrt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Cauchy(location = 2, dispersion = 1)',ylim = c(0,15),lwd = 4)
    lines(rowSums(log(GrunwaldRam)/B),col = 'blue' ,lwd = 2)
    lines(rowSums(log(GrunwaldCorr)/B),col = 'red', lwd = 2)
    lines(rowSums(log(GrunwaldEfron)/B),col = 'green', lwd = 2)
    EfronbadCauchyGrun = recordPlot()
    legend(0,15,c('Gr.E_Peña','Gr.-E_Ram','Gr.E_Corr','Gr.E_Cosh'),text.col = c('green','blue','red','black'),text.width = 12)
  }
}

#FIGURE 3.4
{
{
  library(invgamma)
  library(rootSolve)
  B = 1000
  N = 50
  RamRam = matrix(0,nrow = N,ncol = B)
  RamCosh = matrix(0,nrow = N,ncol = B)
  GrunwaldRam = matrix(0,nrow = N,ncol = B)
  GrunwaldCosh = matrix(0,nrow = N,ncol = B)
  NormRam = matrix(0,nrow = N,ncol = B)
  NormCosh = matrix(0,nrow = N,ncol = B)
}
for(b in 1:B){
  z = rnorm(N,1)
  RamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
  RamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
  GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
  GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
  NormRamt = EfronsScaleInvNormlog(z,1,type = 'Ram',twosided = TRUE)[[6]]
  NormCosht = EfronsScaleInvNormlog(z,1,type = 'NewE',twosided = TRUE)[[6]]
  for(n in 1:N){
    RamRam[n,b] = prod(RamRamt[1:n])
    RamCosh[n,b] = prod(RamCosht[1:n])
    GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
    GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
    NormRam[n,b] = prod(NormRamt[1:n])
    NormCosh[n,b] = prod(NormCosht[1:n])
  }
  print(b)
}
{
  plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X ~ Normal(mean = 1, variance = 1)',ylim = c(0,20),lwd = 2)
  lines(rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
  lines(rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
  lines(rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
  lines(rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
  lines(rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
  Norm1 = recordPlot()
  legend(0,20,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 10)
}
}

#FIGURE 3.5
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    RamRam = matrix(0,nrow = N,ncol = B)
    RamCosh = matrix(0,nrow = N,ncol = B)
    GrunwaldRam = matrix(0,nrow = N,ncol = B)
    GrunwaldCosh = matrix(0,nrow = N,ncol = B)
    NormRam = matrix(0,nrow = N,ncol = B)
    NormCosh = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rnorm(N,2)
    RamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
    NormRamt = EfronsScaleInvNormlog(z,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(z,1,type = 'NewE',twosided = TRUE)[[6]]
    for(n in 1:N){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X ~ Normal(mean = 2, variance = 1)',ylim = c(0,30),lwd = 2)
    lines(rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
    Norm2 = recordPlot()
    legend(0,31,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 10)
  }
}

#FIGURE 3.6
{
{
  thetagrid = seq(0.20,3,0.10)
  library(invgamma)
  library(rootSolve)
  B = 1000
  mRamdasRam = numeric(length(thetagrid))
  mRamdasCosh = numeric(length(thetagrid))
  mGrunwaldRam = numeric(length(thetagrid))
  mGrunwaldCosh = numeric(length(thetagrid))
  mNormRam = numeric(length(thetagrid))
  mNormCosh = numeric(length(thetagrid))
  wRamdasRam = numeric(length(thetagrid))
  wRamdasCosh = numeric(length(thetagrid))
  wGrunwaldRam = numeric(length(thetagrid))
  wGrunwaldCosh = numeric(length(thetagrid))
  wNormRam = numeric(length(thetagrid))
  wNormCosh = numeric(length(thetagrid))
  RamdasRamrec = numeric(B)
  RamdasCoshrec = numeric(B)
  GrunwaldRamrec = numeric(B)
  GrunwaldCoshrec = numeric(B)
  NormRamrec = numeric(B)
  NormCoshrec = numeric(B)
}

for(j in 1:length(thetagrid)){
  for(b in 1:B){
    a = rnorm(100000,thetagrid[j])
    c = rnorm(100000)
    z = a-c
    RamdasRamrec[b] = RamdasTest(z,Eval = 'Ram')[[3]]
    RamdasCoshrec[b] = RamdasTest(z,Eval = 'NewE')[[3]]
    GrunwaldRamrec[b] = EfronsGExp(z,type = 'Ram',twosided = T)[[1]]
    GrunwaldCoshrec[b] = EfronsGExp(z,type = 'NewE',twosided = T)[[1]]
    NormRamrec[b] = EfronsScaleInvNorm(z,1,type = 'Ram',twosided = TRUE)[[1]]
    NormCoshrec[b] = EfronsScaleInvNorm(z,1,type = 'NewE',twosided = TRUE)[[1]]
    print(b)
  }
  wRamdasRam[j] = quantile(RamdasRamrec,0.8)
  wRamdasCosh[j] = quantile(RamdasCoshrec,0.8)
  wGrunwaldRam[j] = quantile(GrunwaldRamrec,0.8)
  wGrunwaldCosh[j] = quantile(GrunwaldCoshrec,0.8)
  wNormRam[j] = quantile(NormRamrec,0.8)
  wNormCosh[j] = quantile(NormCoshrec,0.8)
  RamdasRamrec[RamdasRamrec>quantile(RamdasRamrec,0.8)] = quantile(RamdasRamrec,0.8)
  mRamdasRam[j] = mean(RamdasRamrec)
  RamdasCoshrec[RamdasCoshrec>quantile(RamdasCoshrec,0.8)] = quantile(RamdasCoshrec,0.8)
  mRamdasCosh[j] = mean(RamdasCoshrec)
  GrunwaldRamrec[GrunwaldRamrec>quantile(GrunwaldRamrec,0.8)] = quantile(GrunwaldRamrec,0.8)
  mGrunwaldRam[j] = mean(GrunwaldRamrec)
  GrunwaldCoshrec[GrunwaldCoshrec>quantile(GrunwaldCoshrec,0.8)] = quantile(GrunwaldCoshrec,0.8)
  mGrunwaldCosh[j] = mean(GrunwaldCoshrec)
  NormRamrec[NormRamrec>quantile(NormRamrec,0.8)] = quantile(NormRamrec,0.8)
  mNormRam[j] = mean(NormRamrec)
  NormCoshrec[NormCoshrec>quantile(NormCoshrec,0.8)] = quantile(NormCoshrec,0.8)
  mNormCosh[j] = mean(NormCoshrec)
  print(thetagrid[j])
}
{
  plot(thetagrid,mGrunwaldCosh/regmwurecrec,xlim = c(0.25,3),ylim = c(0,4),col = 'blue',xlab = 'Delta',ylab = 'Comparative average stopping time',type = 'l',lwd = 3,main = 'X ~ Normal(variance = 1)')
  lines(thetagrid,mGrunwaldRam/regmwurecrec,col = 'green',lwd =2)
  lines(thetagrid,mRamdasRam/regmwurecrec,col = 'black' ,lwd = 2)
  lines(thetagrid,mRamdasCosh/regmwurecrec,col = 'red',lwd = 2)
  lines(thetagrid,mNormRam/regmwurecrec,col = 'orange',lwd = 2)
  lines(thetagrid,mNormCosh/regmwurecrec, col = 'purple',lwd = 2)
  legend(1.75,4,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 1)
}

NormMStop = recordPlot()
{
  plot(thetagrid,wGrunwaldCosh/regmwurecrec,xlim = c(0.25,3),ylim = c(0,5),col = 'blue',xlab = 'Delta',ylab = 'Comparative worst case stopping time',type = 'l',lwd = 3,main = 'X ~ Normal(variance = 1)')
  lines(thetagrid,wGrunwaldRam/regmwurecrec,col = 'green',lwd =2)
  lines(thetagrid,wRamdasRam/regmwurecrec,col = 'black' ,lwd = 2)
  lines(thetagrid,wRamdasCosh/regmwurecrec,col = 'red',lwd = 2)
  lines(thetagrid,wNormRam/regmwurecrec,col = 'orange',lwd = 2)
  lines(thetagrid,wNormCosh/regmwurecrec, col = 'purple',lwd = 2)
  legend(1.75,5,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 1)
}

NormWStop = recordPlot()
}

#FIGURE 3.7
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    MaxLCosh = matrix(0,nrow = N,ncol = B)
    MaxLRam = matrix(0,nrow = N,ncol = B)
    GrunwaldRam = matrix(0,nrow = N,ncol = B)
    GrunwaldCosh = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rnorm(N,1)
    MaxLCosht = Efronspowertestmaxllog(z, twosided = T)[[6]]
    MaxLRamt = Efronspowertestmaxllog(z,twosided = T)[[6]]
    GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
    for(n in 1:N){
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
    }
    print(b)
  }
  {
    cols = rainbow(4)
    plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'green',main = 'X ~ Normal(mean = 1, variance = 1)',ylim = c(0,15),lwd = 4)
    lines(rowSums(log(GrunwaldRam)/B),col = 'blue' ,lwd = 2)
    lines(rowSums(log(MaxLCosh)/B),col = 'black' ,lwd = 2)
    lines(rowSums(log(MaxLRam)/B),col = 'red' ,lwd = 2)
    MaxLBad = recordPlot()
    legend(0,15,c('Max.L.Gr.E_Cosh', 'Max.L.Gr.E_Ram', 'Gr.E_Ram','Gr.E_Cosh'),text.col = c('black','red','blue','green'),text.width = 12)
  }
}

## SEMI-HEAVY TAILED DISTRIBUTIONS ##

#FIGURE 3.8
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    RamRam = matrix(0,nrow = N,ncol = B)
    RamCosh = matrix(0,nrow = N,ncol = B)
    GrunwaldRam = matrix(0,nrow = N,ncol = B)
    GrunwaldCosh = matrix(0,nrow = N,ncol = B)
    NormRam = matrix(0,nrow = N,ncol = B)
    NormCosh = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rt(N,3,1)
    RamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
    NormRamt = EfronsScaleInvNormlog(z,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(z,1,type = 'NewE',twosided = TRUE)[[6]]
    for(n in 1:N){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X ~ Student`s t(ncp = 1, df = 3)',ylim = c(-2,20),lwd = 2)
    lines(rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
    T1 = recordPlot()
    legend(0,20,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 10)
  }
}
#FIGURE 3.9
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    RamRam = matrix(0,nrow = N,ncol = B)
    RamCosh = matrix(0,nrow = N,ncol = B)
    GrunwaldRam = matrix(0,nrow = N,ncol = B)
    GrunwaldCosh = matrix(0,nrow = N,ncol = B)
    NormRam = matrix(0,nrow = N,ncol = B)
    NormCosh = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rt(N,3,2)
    RamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
    NormRamt = EfronsScaleInvNormlog(z,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(z,1,type = 'NewE',twosided = TRUE)[[6]]
    for(n in 1:N){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X ~ Student`s t(ncp = 2, df = 3)',ylim = c(0,30),lwd = 2)
    lines(rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
    T2 = recordPlot()
    legend(0,31,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 10)
  }
}

#FIGURE 3.10
{
  {
    thetagrid = seq(0.4,3,0.10)
    library(invgamma)
    library(rootSolve)
    B = 200
    mRamdasRam = numeric(length(thetagrid))
    mRamdasCosh = numeric(length(thetagrid))
    mGrunwaldRam = numeric(length(thetagrid))
    mGrunwaldCosh = numeric(length(thetagrid))
    mNormRam = numeric(length(thetagrid))
    mNormCosh = numeric(length(thetagrid))
    wRamdasRam = numeric(length(thetagrid))
    wRamdasCosh = numeric(length(thetagrid))
    wGrunwaldRam = numeric(length(thetagrid))
    wGrunwaldCosh = numeric(length(thetagrid))
    wNormRam = numeric(length(thetagrid))
    wNormCosh = numeric(length(thetagrid))
    RamdasRamrec = numeric(B)
    RamdasCoshrec = numeric(B)
    GrunwaldRamrec = numeric(B)
    GrunwaldCoshrec = numeric(B)
    NormRamrec = numeric(B)
    NormCoshrec = numeric(B)
  }
  
  for(j in 1:length(thetagrid)){
    for(b in 1:B){
      a = rt(10000,3,thetagrid[j])
      c = rt(10000,3)
      z = a-c
      RamdasRamrec[b] = RamdasTest(z,Eval = 'Ram')[[3]]
      RamdasCoshrec[b] = RamdasTest(z,Eval = 'NewE')[[3]]
      GrunwaldRamrec[b] = EfronsGExp(z,type = 'Ram',twosided = T)[[1]]
      GrunwaldCoshrec[b] = EfronsGExp(z,type = 'NewE',twosided = T)[[1]]
      NormRamrec[b] = EfronsScaleInvNorm(z,1,type = 'Ram',twosided = TRUE)[[1]]
      NormCoshrec[b] = EfronsScaleInvNorm(z,1,type = 'NewE',twosided = TRUE)[[1]]
      print(b)
    }
    wRamdasRam[j] = quantile(RamdasRamrec,0.8)
    wRamdasCosh[j] = quantile(RamdasCoshrec,0.8)
    wGrunwaldRam[j] = quantile(GrunwaldRamrec,0.8)
    wGrunwaldCosh[j] = quantile(GrunwaldCoshrec,0.8)
    wNormRam[j] = quantile(NormRamrec,0.8)
    wNormCosh[j] = quantile(NormCoshrec,0.8)
    RamdasRamrec[RamdasRamrec>quantile(RamdasRamrec,0.8)] = quantile(RamdasRamrec,0.8)
    mRamdasRam[j] = mean(RamdasRamrec)
    RamdasCoshrec[RamdasCoshrec>quantile(RamdasCoshrec,0.8)] = quantile(RamdasCoshrec,0.8)
    mRamdasCosh[j] = mean(RamdasCoshrec)
    GrunwaldRamrec[GrunwaldRamrec>quantile(GrunwaldRamrec,0.8)] = quantile(GrunwaldRamrec,0.8)
    mGrunwaldRam[j] = mean(GrunwaldRamrec)
    GrunwaldCoshrec[GrunwaldCoshrec>quantile(GrunwaldCoshrec,0.8)] = quantile(GrunwaldCoshrec,0.8)
    mGrunwaldCosh[j] = mean(GrunwaldCoshrec)
    NormRamrec[NormRamrec>quantile(NormRamrec,0.8)] = quantile(NormRamrec,0.8)
    mNormRam[j] = mean(NormRamrec)
    NormCoshrec[NormCoshrec>quantile(NormCoshrec,0.8)] = quantile(NormCoshrec,0.8)
    mNormCosh[j] = mean(NormCoshrec)
    print(thetagrid[j])
  }
  {
    plot(thetagrid,mGrunwaldCosh/regmwurecrecT[3:length(regmwurecrecT)],xlim = c(0.5,3),ylim = c(0,6),col = 'blue',xlab = 'Delta',ylab = 'Comparative average stopping time',type = 'l',lwd = 3,main = 'X ~ Student`s t(ncp = Delta, df = 3)')
    lines(thetagrid,mGrunwaldRam/regmwurecrecT[3:length(regmwurecrecT)],col = 'green',lwd =2)
    lines(thetagrid,mRamdasRam/regmwurecrecT[3:length(regmwurecrecT)],col = 'black' ,lwd = 2)
    lines(thetagrid,mRamdasCosh/regmwurecrecT[3:length(regmwurecrecT)],col = 'red',lwd = 2)
    lines(thetagrid,mNormRam/regmwurecrecT[3:length(regmwurecrecT)],col = 'orange',lwd = 2)
    lines(thetagrid,mNormCosh/regmwurecrecT[3:length(regmwurecrecT)], col = 'purple',lwd = 2)
    legend(1.75,6,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 1)
  }
  TMSTOP = recordPlot()
  {
    plot(thetagrid,wGrunwaldCosh/regmwurecrecT[3:length(regmwurecrecT)],xlim = c(0.5,3),ylim = c(0,8.5),col = 'blue',xlab = 'Delta',ylab = 'Comparative worst case stopping time',type = 'l',lwd = 3,main = 'X ~ Student`s t(ncp = Delta, df = 3)')
    lines(thetagrid,wGrunwaldRam/regmwurecrecT[3:length(regmwurecrecT)],col = 'green',lwd =2)
    lines(thetagrid,wRamdasRam/regmwurecrecT[3:length(regmwurecrecT)],col = 'black' ,lwd = 2)
    lines(thetagrid,wRamdasCosh/regmwurecrecT[3:length(regmwurecrecT)],col = 'red',lwd = 2)
    lines(thetagrid,wNormRam/regmwurecrecT[3:length(regmwurecrecT)],col = 'orange',lwd = 2)
    lines(thetagrid,wNormCosh/regmwurecrecT[3:length(regmwurecrecT)], col = 'purple',lwd = 2)
    legend(1.75,8.5,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 1)
  }
  
  TWSTOP = recordPlot()
}


## HEAVY TAILED DISTRIBUTIONS ##
# FIGURE 3.11
{
{
  library(invgamma)
  library(rootSolve)
  B = 1000
  N = 50
  RamRam = matrix(0,nrow = N,ncol = B)
  RamCosh = matrix(0,nrow = N,ncol = B)
  GrunwaldRam = matrix(0,nrow = N,ncol = B)
  GrunwaldCosh = matrix(0,nrow = N,ncol = B)
  NormRam = matrix(0,nrow = N,ncol = B)
  NormCosh = matrix(0,nrow = N,ncol = B)
  Allard = matrix(0,nrow = N,ncol = B)
}
for(b in 1:B){
  z = rcauchy(N,1)
  RamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
  RamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
  GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
  GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
  NormRamt = EfronsScaleInvNormlog(z,1,type = 'Ram',twosided = TRUE)[[6]]
  NormCosht = EfronsScaleInvNormlog(z,1,type = 'NewE',twosided = TRUE)[[6]]
  Allardt = Ebernoullilog(z)[[3]]
  for(n in 1:N){
    RamRam[n,b] = prod(RamRamt[1:n])
    RamCosh[n,b] = prod(RamCosht[1:n])
    GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
    GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
    NormRam[n,b] = prod(NormRamt[1:n])
    NormCosh[n,b] = prod(NormCosht[1:n])
    Allard[n,b] = prod(Allardt[1:n])
  }
  print(b)
}
{
  plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X ~ Cauchy(location = 1, dispersion = 1)',ylim = c(-1,5),lwd = 2)
  lines(rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
  lines(rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
  lines(rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
  lines(rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
  lines(rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
  lines(rowSums(log(Allard)/B),col = 'brown',lwd = 2)
}
Cauchy1 = recordPlot()

legend(0,5,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh','E_Bern'),text.col = c('black','red','green','blue','orange','purple','brown'),text.width = 14)
}

# FIGURE 3.12
{
  {
    library(invgamma)
    library(rootSolve)
    B = 10000
    N = 50
    RamRam = matrix(0,nrow = N,ncol = B)
    RamCosh = matrix(0,nrow = N,ncol = B)
    GrunwaldRam = matrix(0,nrow = N,ncol = B)
    GrunwaldCosh = matrix(0,nrow = N,ncol = B)
    NormRam = matrix(0,nrow = N,ncol = B)
    NormCosh = matrix(0,nrow = N,ncol = B)
    Allard = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rcauchy(N,3)
    RamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
    NormRamt = EfronsScaleInvNormlog(z,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(z,1,type = 'NewE',twosided = TRUE)[[6]]
    Allardt = Ebernoullilog(z)[[3]]
    for(n in 1:N){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
      Allard[n,b] = prod(Allardt[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X ~ Cauchy(location = 3, dispersion = 1)',ylim = c(-1,12),lwd = 2)
    lines(rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
    lines(rowSums(log(Allard)/B),col = 'brown',lwd = 2)
  }
  Cauchy3 = recordPlot()
  
  legend(0,12,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh','E_Bern'),text.col = c('black','red','green','blue','orange','purple','brown'),text.width = 14)
}

# FIGURE 3.13
{
{
  thetagrid = seq(0.6,3,0.2)
  library(invgamma)
  library(rootSolve)
  B = 500
  mRamdasRam = numeric(length(thetagrid))
  mRamdasCosh = numeric(length(thetagrid))
  mGrunwaldRam = numeric(length(thetagrid))
  mGrunwaldCosh = numeric(length(thetagrid))
  mNormRam = numeric(length(thetagrid))
  mNormCosh = numeric(length(thetagrid))
  mAllard = numeric(length(thetagrid))
  wRamdasRam = numeric(length(thetagrid))
  wRamdasCosh = numeric(length(thetagrid))
  wGrunwaldRam = numeric(length(thetagrid))
  wGrunwaldCosh = numeric(length(thetagrid))
  wNormRam = numeric(length(thetagrid))
  wNormCosh = numeric(length(thetagrid))
  wAllard = numeric(length(thetagrid))
  RamdasRamrec = numeric(B)
  RamdasCoshrec = numeric(B)
  GrunwaldRamrec = numeric(B)
  GrunwaldCoshrec = numeric(B)
  NormRamrec = numeric(B)
  NormCoshrec = numeric(B)
  Allardrec = numeric(B)
}

for(j in 1:length(thetagrid)){
  for(b in 1:B){
    a = rcauchy(1000,thetagrid[j])
    c = rcauchy(1000)
    z = a-c
    #RamdasRamrec[b] = RamdasTest(z,Eval = 'Ram')[[3]]
    #RamdasCoshrec[b] = RamdasTest(z,Eval = 'NewE')[[3]]
    #GrunwaldRamrec[b] = EfronsGExp(z,type = 'Ram',twosided = T)[[1]]
    #GrunwaldCoshrec[b] = EfronsGExp(z,type = 'NewE',twosided = T)[[1]]
    #NormRamrec[b] = EfronsScaleInvNorm(z,1,type = 'Ram',twosided = TRUE)[[1]]
    #NormCoshrec[b] = EfronsScaleInvNorm(z,1,type = 'NewE',twosided = TRUE)[[1]]
    Allardrec[b] = Ebernoulli(z)[[1]]
    print(b)
  }
  wRamdasRam[j] = quantile(RamdasRamrec,0.8)
  wRamdasCosh[j] = quantile(RamdasCoshrec,0.8)
  wGrunwaldRam[j] = quantile(GrunwaldRamrec,0.8)
  wGrunwaldCosh[j] = quantile(GrunwaldCoshrec,0.8)
  wNormRam[j] = quantile(NormRamrec,0.8)
  wNormCosh[j] = quantile(NormCoshrec,0.8)
  wAllard[j] = quantile(Allardrec,0.8)
  RamdasRamrec[RamdasRamrec>quantile(RamdasRamrec,0.8)] = quantile(RamdasRamrec,0.8)
  mRamdasRam[j] = mean(RamdasRamrec)
  RamdasCoshrec[RamdasCoshrec>quantile(RamdasCoshrec,0.8)] = quantile(RamdasCoshrec,0.8)
  mRamdasCosh[j] = mean(RamdasCoshrec)
  GrunwaldRamrec[GrunwaldRamrec>quantile(GrunwaldRamrec,0.8)] = quantile(GrunwaldRamrec,0.8)
  mGrunwaldRam[j] = mean(GrunwaldRamrec)
  GrunwaldCoshrec[GrunwaldCoshrec>quantile(GrunwaldCoshrec,0.8)] = quantile(GrunwaldCoshrec,0.8)
  mGrunwaldCosh[j] = mean(GrunwaldCoshrec)
  NormRamrec[NormRamrec>quantile(NormRamrec,0.8)] = quantile(NormRamrec,0.8)
  mNormRam[j] = mean(NormRamrec)
  NormCoshrec[NormCoshrec>quantile(NormCoshrec,0.8)] = quantile(NormCoshrec,0.8)
  mNormCosh[j] = mean(NormCoshrec)
  Allardrec[Allardrec>quantile(Allardrec,0.8)] = quantile(Allardrec,0.8)
  mAllard[j] = mean(Allardrec)
  print(thetagrid[j])
}
  mGrunwaldCosh[wGrunwaldCosh = 3000] = 3000 #CHECKING IF LIMIT REACHED
  mGrunwaldRam[wGrunwaldRam = 3000] = 3000
  mRamdasRam[wRamdasRam = 3000] = 3000
  mRamdasCosh[wRamdasCosh = 3000] = 3000
  mNormRam[wNormRam = 3000] = 3000
  mNormCosh[wNormCosh = 3000] = 3000
  mAllard[wAllard = 3000] = 3000
  
{
  plot(thetagrid,mGrunwaldCosh[1:13]/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],xlim = c(0.5,3),ylim = c(1,8),col = 'blue',xlab = 'Delta' ,ylab = 'Average stopping time n',type = 'l',lwd = 2,main = 'Xa ~ Cauchy(location = 0), Xb ~ Cauchy(location = Delta)')
  lines(thetagrid,mGrunwaldRam[1:13]/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],col = 'green',lwd =2)
  lines(thetagrid,mRamdasRam[1:13]/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],col = 'black' ,lwd = 2)
  lines(thetagrid,mRamdasCosh[1:13]/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],col = 'red',lwd = 2)
  lines(thetagrid,mNormRam[1:13]/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],col = 'orange',lwd = 2)
  lines(thetagrid,mNormCosh[1:13]/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)], col = 'purple',lwd = 2)
  lines(thetagrid,mAllard[1:13]/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)], col = 'brown',lwd = 2)
  legend(2.2,8,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh','Allard'),text.col = c('black','red','green','blue','orange','purple','brown'),text.width = 0.65)
}

CauchyMStop = recordPlot()

{
  plot(thetagrid,wGrunwaldCosh/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],xlim = c(0.5,3),ylim = c(1,20),col = 'blue',xlab = 'Delta',ylab = 'Worst case stopping time n',type = 'l',lwd = 2,main ='Xa ~ Cauchy(location = 0), Xb ~ Cauchy(location = Delta)')
  lines(thetagrid,wGrunwaldRam/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],col = 'green',lwd =2)
  lines(thetagrid,wRamdasRam/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],col = 'black' ,lwd = 2)
  lines(thetagrid,wRamdasCosh/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],col = 'red',lwd = 2)
  lines(thetagrid,wNormRam/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)],col = 'orange',lwd = 2)
  lines(thetagrid,wNormCosh/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)], col = 'purple',lwd = 2)
  lines(thetagrid,wAllard/regmwurecrecCauch[seq(5,length(regmwurecrecCauch),2)], col = 'brown',lwd = 2)
  legend(0.5,20.5,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh','E_Bern'),text.col = c('black','red','green','blue','orange','purple','brown'),text.width = 0.65)
}

CauchyWStop = recordPlot()
}

## HEDGE METHOD ASSESSMENT ##

# FIGURE 3.14
{
  {
    library(invgamma)
    library(rootSolve)
    B = 2000
    N = 50
    hGrunCosh =  matrix(0,nrow = N,ncol = B)
    hGrunCosh01 =  matrix(0,nrow = N,ncol = B)
    hGrunCosh02 =  matrix(0,nrow = N,ncol = B)
    hGrunCosh04 =  matrix(0,nrow = N,ncol = B)
    Bern = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rcauchy(N,1)
    hGrunwaldCosht = EfronsGExpWlog(z,type = 'NewE',twosided = T,0)[[6]]
    hGrunwaldCosh01t = EfronsGExpWlog(z,type = 'NewE',twosided = T,0.1)[[6]]
    hGrunwaldCosh02t = EfronsGExpWlog(z,type = 'NewE',twosided = T,0.2)[[6]]
    hGrunwaldCosh04t = EfronsGExpWlog(z,type = 'NewE',twosided = T,0.4)[[6]]
    Bernt = Ebernoullilog(z)[[3]]
    for(n in 1:N){
      hGrunCosh[n,b] = prod(hGrunwaldCosht[1:n])
      hGrunCosh01[n,b] = prod(hGrunwaldCosh01t[1:n])
      hGrunCosh02[n,b] = prod(hGrunwaldCosh02t[1:n])
      hGrunCosh04[n,b] = prod(hGrunwaldCosh04t[1:n])
      Bern[n,b] = prod(Bernt[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(hGrunCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Cauchy(location = 1, dispersion = 1)',ylim = c(-1,5),lwd = 2)
    lines(rowSums(log(hGrunCosh01)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(hGrunCosh02)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(hGrunCosh04)/B),col = 'blue' ,lwd = 2)
    lines(rowSums(log(Bern)/B), col = 'brown', lwd = 2)
  }
  FinalCauchyCompare = recordPlot()
  legend(0,5,c('Gr.E_Cosh','Gr.E_Cosh(h=0.1)','Gr.E_Cosh(h=0.2)','Gr.E_Cosh(h=0.4)','E_Bern'),text.col = c('black','red','green','blue','brown'),text.width = 14)
}

#update likelihood without hedge?
# FIGURE 3.15
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    hGrunRam =  matrix(0,nrow = N,ncol = B)
    hGrunCosh =  matrix(0,nrow = N,ncol = B)
    hGrunRamNolik =  matrix(0,nrow = N,ncol = B)
    hGrunCoshNolik =  matrix(0,nrow = N,ncol = B)
#    hRamRam =  matrix(0,nrow = N,ncol = B)
#    hRamCosh =  matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rcauchy(N,1)
#    hRamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
#    hRamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
    hGrunwaldRamt = EfronsGExpWlog(z,type = 'Ram',twosided = T,0.2)[[6]]
    hGrunwaldCosht = EfronsGExpWlog(z,type = 'NewE',twosided = T,0.2)[[6]]
    hGrunwaldRamtNolik = EfronsGExpWlognolik(z,type = 'Ram',twosided = T,0.2)[[6]]
    hGrunwaldCoshtNolik = EfronsGExpWlognolik(z,type = 'NewE',twosided = T,0.2)[[6]]
    for(n in 1:N){
#      RamRam[n,b] = prod(RamRamt[1:n])
#      RamCosh[n,b] = prod(RamCosht[1:n])
      hGrunRam[n,b] = prod(hGrunwaldRamt[1:n])
      hGrunCosh[n,b] = prod(hGrunwaldCosht[1:n])
      hGrunRamNolik[n,b] = prod(hGrunwaldRamtNolik[1:n])
      hGrunCoshNolik[n,b] = prod(hGrunwaldCoshtNolik[1:n])
      
    }
    print(b)
  }
  {
    plot(rowSums(log(hGrunRam)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Cauchy(location = 1, dispersion = 1)',ylim = c(-1,5),lwd = 2)
    lines(rowSums(log(hGrunCosh)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(hGrunRamNolik)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(hGrunCoshNolik)/B),col = 'blue' ,lwd = 2)
  }
  Hedgenolik = recordPlot()
  
  legend(0,5,c('Gr.E_Ram.(h=0.2)','Gr.E_Cosh(h=0.2)','Gr.E_Ram.(No Lik)','Gr.E_Cosh(No Lik)'),text.col = c('black','red','green','blue'),text.width = 14)
}

#update E-variable without hedge?
# FIGURE 3.16
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    hGrunRam =  matrix(0,nrow = N,ncol = B)
    hGrunCosh =  matrix(0,nrow = N,ncol = B)
    hGrunRamNolik =  matrix(0,nrow = N,ncol = B)
    hGrunCoshNolik =  matrix(0,nrow = N,ncol = B)
    hGrunRam3 =  matrix(0,nrow = N,ncol = B)
    hGrunCosh3 =  matrix(0,nrow = N,ncol = B)
    #    hRamRam =  matrix(0,nrow = N,ncol = B)
    #    hRamCosh =  matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rcauchy(N,1)
    #    hRamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
    #    hRamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
    hGrunwaldRamt = EfronsGExpWlog(z,type = 'Ram',twosided = T,0.2)[[6]]
    hGrunwaldCosht = EfronsGExpWlog(z,type = 'NewE',twosided = T,0.2)[[6]]
    hGrunwaldRamtNolik = EfronsGExpWlognolik(z,type = 'Ram',twosided = T,0.2)[[6]]
    hGrunwaldCoshtNolik = EfronsGExpWlognolik(z,type = 'NewE',twosided = T,0.2)[[6]]
    hGrunwaldRamt3 = EfronsGExpWlog3(z,type = 'Ram',twosided = T,0.2)[[6]]
    hGrunwaldCosht3 = EfronsGExpWlog3(z,type = 'NewE',twosided = T,0.2)[[6]]
    for(n in 1:N){
      #      RamRam[n,b] = prod(RamRamt[1:n])
      #      RamCosh[n,b] = prod(RamCosht[1:n])
      hGrunRam[n,b] = prod(hGrunwaldRamt[1:n])
      hGrunCosh[n,b] = prod(hGrunwaldCosht[1:n])
      hGrunRamNolik[n,b] = prod(hGrunwaldRamtNolik[1:n])
      hGrunCoshNolik[n,b] = prod(hGrunwaldCoshtNolik[1:n])
      hGrunRam3[n,b] = prod(hGrunwaldRamt3[1:n])
      hGrunCosh3[n,b] = prod(hGrunwaldCosht3[1:n])
      
    }
    print(b)
  }
  {
    plot(rowSums(log(hGrunRam)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Cauchy(location = 1, dispersion = 1)',ylim = c(-1,5),lwd = 2)
    lines(rowSums(log(hGrunCosh)/B),col = 'red' ,lwd = 2)
#    lines(rowSums(log(hGrunRamNolik)/B),col = 'green' ,lwd = 2)
#    lines(rowSums(log(hGrunCoshNolik)/B),col = 'blue' ,lwd = 2)
    lines(rowSums(log(hGrunRam3)/B),col = 'orange' ,lwd = 2)
    lines(rowSums(log(hGrunCosh3)/B),col = 'purple' ,lwd = 2)
  }
  NoHedge = recordPlot()

  legend(0,5,c('Gr.E_Ram.(h=0.2)','Gr.E_Cosh(h=0.2)','Gr.E_Ram.(No Lik)','Gr.E_Cosh(No Lik)','Gr.E_Ram.(No hedge)','Gr.E_Cosh(No hedge)'),text.col = c('black','red','green','blue','orange','purple'),text.width = 14)
}

#estimate h from data?
# FIGURE 3.17
{
  {
    library(invgamma)
    library(rootSolve)
    B = 2000
    N = 50
    hGrunRam =  matrix(0,nrow = N,ncol = B)
    hGrunCosh =  matrix(0,nrow = N,ncol = B)
    hGrunRamlearn =  matrix(0,nrow = N,ncol = B)
    hGrunCoshlearn =  matrix(0,nrow = N,ncol = B)
    GrunRam =  matrix(0,nrow = N,ncol = B)
    GrunCosh =  matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rcauchy(N,1)
    GrunRamt = EfronsGExpWlog(z,type = 'Ram',twosided = T,0)[[6]]
    GrunCosht = EfronsGExpWlog(z,type = 'NewE',twosided = T,0)[[6]]
    hGrunwaldRamt = EfronsGExpWlog(z,type = 'Ram',twosided = T,0.2)[[6]]
    hGrunwaldCosht = EfronsGExpWlog(z,type = 'NewE',twosided = T,0.2)[[6]]
    hGrunRamlearnt = EfronsGExpVlog(z,type = 'Ram',twosided = T,0.2)[[6]]
    hgrunCoshlearnt = EfronsGExpVlog(z,type = 'NewE',twosided = T,0.2)[[6]]
    for(n in 1:N){
      GrunRam[n,b] = prod(GrunRamt[1:n])
      GrunCosh[n,b] = prod(GrunCosht[1:n])
      hGrunRam[n,b] = prod(hGrunwaldRamt[1:n])
      hGrunCosh[n,b] = prod(hGrunwaldCosht[1:n])
      hGrunRamlearn[n,b] = prod(hGrunRamlearnt[1:n])
      hGrunCoshlearn[n,b] = prod(hgrunCoshlearnt[1:n])
      
    }
    print(b)
  }
  {
    plot(rowSums(log(hGrunRam)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Cauchy(location = 1, dispersion = 1)',ylim = c(-1,5),lwd = 2)
    lines(rowSums(log(hGrunCosh)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(hGrunRamlearn)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(hGrunCoshlearn)/B),col = 'blue' ,lwd = 2)
    lines(rowSums(log(GrunRam)/B),col = 'orange' ,lwd = 2)
    lines(rowSums(log(GrunCosh)/B),col = 'purple' ,lwd = 2)
  }
  Estimateh = recordPlot()
  
  legend(0,5,c('Gr.E_Ram.(h=0.2)','Gr.E_Cosh(h=0.2)','Gr.E_Ram.(learn h)','Gr.E_Cosh(learn h)','Gr.E_Ram.(h=0)','Gr.E_Cosh(h=0)'),text.col = c('black','red','green','blue','orange','purple'),text.width = 14)
}

# FIGURE 3.18
{
  {
    library(invgamma)
    library(rootSolve)
    B = 2000
    N = 100
    hGrunRam =  matrix(0,nrow = N,ncol = B)
    hGrunCosh =  matrix(0,nrow = N,ncol = B)
    hGrunRam02 =  matrix(0,nrow = N,ncol = B)
    hGrunCosh02 =  matrix(0,nrow = N,ncol = B)
    hGrunRamlearn =  matrix(0,nrow = N,ncol = B)
    hgrunCoshlearn =  matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rcauchy(N,1)
    hGrunwaldRamt = EfronsGExpWlog(z,type = 'Ram',twosided = T,0)[[6]]
    hGrunwaldCosht = EfronsGExpWlog(z,type = 'NewE',twosided = T,0)[[6]]
    hGrunwaldRam02t = EfronsGExpWlog(z,type = 'Ram',twosided = T,0.2)[[6]]
    hGrunwaldCosh02t = EfronsGExpWlog(z,type = 'NewE',twosided = T,0.2)[[6]]
    hGrunRamlearnt = EfronsGExpVlog2point(z,type = 'Ram',twosided = T)[[6]]
    hgrunCoshlearnt = EfronsGExpVlog2point(z,type = 'NewE',twosided = T)[[6]]
    for(n in 1:N){
      hGrunRam02[n,b] = prod(hGrunwaldRam02t[1:n])
      hGrunCosh02[n,b] = prod(hGrunwaldCosh02t[1:n])
      hGrunRam[n,b] = prod(hGrunwaldRamt[1:n])
      hGrunCosh[n,b] = prod(hGrunwaldCosht[1:n])
      hGrunRamlearn[n,b] = prod(hGrunRamlearnt[1:n])
      hgrunCoshlearn[n,b] = prod(hgrunCoshlearnt[1:n])
      
    }
    print(b)
  }
  {
    plot(rowSums(log(hGrunRam)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Cauchy(location = 1, dispersion = 1)',ylim = c(-1,5.5),lwd = 2)
    lines(rowSums(log(hGrunCosh)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(hGrunRam02)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(hGrunCosh02)/B),col = 'blue' ,lwd = 2)
    lines(rowSums(log(hGrunRamlearn)/B),col = 'brown' ,lwd = 2)
    lines(rowSums(log(hgrunCoshlearn)/B),col = 'purple' ,lwd = 2)
  }
  Hedge2pointcompare = recordPlot()
  legend(0,5.5,c('Gr.E_Ram.','Gr.E_Cosh','Gr.E_Ram.(h=0.2)','Gr.E_Cosh(h=0.2)','Gr.E_Ram.(2 point h)','Gr.E_Cosh(2 point h)'),text.col = c('black','red','green','blue','brown','purple'),text.width = 28)
}

#Maximum likelihood with h?
# FIGURE 3.19
{
  {
    library(invgamma)
    library(rootSolve)
    B = 2000
    N = 50
    hGrunCosh =  matrix(0,nrow = N,ncol = B)
    hGrunCosh02 =  matrix(0,nrow = N,ncol = B)
    hMaxL =  matrix(0,nrow = N,ncol = B)
    hMaxL02 =  matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rcauchy(N,1)
    hGrunwaldCosht = EfronsGExpWlog(z,type = 'NewE',twosided = T,0)[[6]]
    hGrunwaldCosh02t = EfronsGExpWlog(z,type = 'NewE',twosided = T,0.2)[[6]]
    hMaxLt = Efronspowertestmaxllog(z, twosided = T,h=0)[[6]]
    hMaxL02t = Efronspowertestmaxllog(z, twosided = T,h=0.2)[[6]]
    for(n in 1:N){
      hGrunCosh02[n,b] = prod(hGrunwaldCosh02t[1:n])
      hGrunCosh[n,b] = prod(hGrunwaldCosht[1:n])
      hMaxL02[n,b] = prod(hMaxL02t[1:n])
      hMaxL[n,b] = prod(hMaxLt[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(hGrunCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Cauchy(location = 1, dispersion = 1)',ylim = c(-1,5),lwd = 2)
    lines(rowSums(log(hGrunCosh02)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(hMaxL)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(hMaxL02)/B),col = 'blue' ,lwd = 2)
  }
  Hedge2pointcompare = recordPlot()
  legend(0,5,c('Gr.E_Cosh','Gr.E_Cosh(h=0.2)','MaxL.E_Cosh','MaxL.E_Cosh(h=0.2)'),text.col = c('black','red','green','blue'),text.width = 14)
}

# FIGURE 3.20
{
  {
    library(invgamma)
    library(rootSolve)
    B = 2000
    N = 50
    hGrunCosh =  matrix(0,nrow = N,ncol = B)
    hGrunCosh02 =  matrix(0,nrow = N,ncol = B)
    hMaxL =  matrix(0,nrow = N,ncol = B)
    hMaxL02 =  matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rnorm(N,1)
    hGrunwaldCosht = EfronsGExpWlog(z,type = 'NewE',twosided = T,0)[[6]]
    hGrunwaldCosh02t = EfronsGExpWlog(z,type = 'NewE',twosided = T,0.2)[[6]]
    hMaxLt = Efronspowertestmaxllog(z, twosided = T,h=0)[[6]]
    hMaxL02t = Efronspowertestmaxllog(z, twosided = T,h=0.2)[[6]]
    for(n in 1:N){
      hGrunCosh02[n,b] = prod(hGrunwaldCosh02t[1:n])
      hGrunCosh[n,b] = prod(hGrunwaldCosht[1:n])
      hMaxL02[n,b] = prod(hMaxL02t[1:n])
      hMaxL[n,b] = prod(hMaxLt[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(hGrunCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'black',main = 'X ~ Normal(mean = 1, variance = 1)',ylim = c(-1,15),lwd = 2)
    lines(rowSums(log(hGrunCosh02)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(hMaxL)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(hMaxL02)/B),col = 'blue' ,lwd = 2)
  }
  HedgeMaxLNorm = recordPlot()
  legend(0,15,c('Gr.E_Cosh','Gr.E_Cosh(h=0.2)','MaxL.E_Cosh','MaxL.E_Cosh(h=0.2)'),text.col = c('black','red','green','blue'),text.width = 14)
}

#FIGURE 3.21
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    RamRam = matrix(0,nrow = N,ncol = B)
    RamCosh = matrix(0,nrow = N,ncol = B)
    GrunwaldRam = matrix(0,nrow = N,ncol = B)
    GrunwaldCosh = matrix(0,nrow = N,ncol = B)
    NormRam = matrix(0,nrow = N,ncol = B)
    NormCosh = matrix(0,nrow = N,ncol = B)
    
    WRamRam = matrix(0,nrow = N,ncol = B)
    WRamCosh = matrix(0,nrow = N,ncol = B)
    WGrunwaldRam = matrix(0,nrow = N,ncol = B)
    WGrunwaldCosh = matrix(0,nrow = N,ncol = B)
    WNormRam = matrix(0,nrow = N,ncol = B)
    WNormCosh = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rnorm(N,1)
    RamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
    NormRamt = EfronsScaleInvNormlog(z,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(z,1,type = 'NewE',twosided = TRUE)[[6]]
    
    WRamRamt = RamdasTestWlog(z, Eval = 'Ram',W = 0.1)[[1]]
    WRamCosht = RamdasTestWlog(z, Eval = 'NewE',W = 0.1)[[1]]
    WGrunwaldRamt = EfronsGExpWlog(z,type = 'Ram',twosided = T,W=0.1)[[6]]
    WGrunwaldCosht = EfronsGExpWlog(z,type = 'NewE',twosided = T,W=0.1)[[6]]
    WNormRamt = EfronsScaleInvNormWlog(z,1,type = 'Ram',twosided = TRUE,W=0.1)[[6]]
    WNormCosht = EfronsScaleInvNormWlog(z,1,type = 'NewE',twosided = TRUE,W=0.1)[[6]]
    for(n in 1:N){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
      
      WRamRam[n,b] = prod(WRamRamt[1:n])
      WRamCosh[n,b] = prod(WRamCosht[1:n])
      WGrunwaldRam[n,b] = prod(WGrunwaldRamt[1:n])
      WGrunwaldCosh[n,b] = prod(WGrunwaldCosht[1:n])
      WNormRam[n,b] = prod(WNormRamt[1:n])
      WNormCosh[n,b] = prod(WNormCosht[1:n])
    }
    print(b)
  }
  {
    cols = rainbow(12)
    plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X ~ Normal(mean = 1, variance = 1)',ylim = c(0,20),lwd = 2)
    lines(rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
    
    lines(rowSums(log(WGrunwaldCosh)/B),col = 'blue' ,lwd = 2,type = 'b')
    lines(rowSums(log(WGrunwaldRam)/B),col = 'green' ,lwd = 2,type = 'b')
    lines(rowSums(log(WRamRam)/B),col = 'black' ,lwd = 2,type = 'b')
    lines(rowSums(log(WRamCosh)/B),col = 'red' ,lwd = 2,type = 'b')
    lines(rowSums(log(WNormRam)/B),col = 'orange' ,lwd = 2,type = 'b')
    lines(rowSums(log(WNormCosh)/B),col = 'purple' ,lwd = 2,type = 'b')
    AllhedgeNorm = recordPlot()
    legend(0,20,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 10)
  }
}


## LOGISTIC ALTERNATIVE ##

#LOGALT PLOTS SECTION 5

# FIGURE 3.22
{
{
  library(invgamma)
  library(rootSolve)
  B = 2000
  N = 50
  RamRam = matrix(0,nrow = N,ncol = B)
  RamCosh = matrix(0,nrow = N,ncol = B)
  GrunwaldRam = matrix(0,nrow = N,ncol = B)
  GrunwaldCosh = matrix(0,nrow = N,ncol = B)
  NormRam = matrix(0,nrow = N,ncol = B)
  NormCosh = matrix(0,nrow = N,ncol = B)
}
for(b in 1:B){
  z = rlogalt(N,1)
  RamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
  RamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
  GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
  GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
  NormRamt = EfronsScaleInvNormlog(z,1,type = 'Ram',twosided = TRUE)[[6]]
  NormCosht = EfronsScaleInvNormlog(z,1,type = 'NewE',twosided = TRUE)[[6]]
  for(n in 1:N){
    RamRam[n,b] = prod(RamRamt[1:n])
    RamCosh[n,b] = prod(RamCosht[1:n])
    GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
    GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
    NormRam[n,b] = prod(NormRamt[1:n])
    NormCosh[n,b] = prod(NormCosht[1:n])
  }
  print(b)
}
{
  plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X ~ Logistic alternative (Delta = 1)',ylim = c(0,5),lwd = 2)
  lines(rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
  lines(rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
  lines(rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
  lines(rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
  lines(rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
}
LogAlt1 = recordPlot()

legend(0,5,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 14)
}

# FIGURE 3.23
{
  {
    library(invgamma)
    library(rootSolve)
    B = 2000
    N = 50
    RamRam = matrix(0,nrow = N,ncol = B)
    RamCosh = matrix(0,nrow = N,ncol = B)
    GrunwaldRam = matrix(0,nrow = N,ncol = B)
    GrunwaldCosh = matrix(0,nrow = N,ncol = B)
    NormRam = matrix(0,nrow = N,ncol = B)
    NormCosh = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    z = rlogalt(N,0.5)
    RamRamt = RamdasTestlog(z, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(z, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(z,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(z,type = 'NewE',twosided = T)[[6]]
    NormRamt = EfronsScaleInvNormlog(z,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(z,1,type = 'NewE',twosided = TRUE)[[6]]
    for(n in 1:N){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(GrunwaldCosh)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X ~ Logistic alternative (Delta = 0.5)',ylim = c(-2,2),lwd = 2)
    lines(rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
  }
  LogAlt05 = recordPlot()
  
  legend(0,2,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 14)
}

# FIGURE 3.24
{
{
  thetagrid = seq(0.4,3,0.1)
  library(invgamma)
  library(rootSolve)
  B = 500
  mRamdasRam = numeric(length(thetagrid))
  mRamdasCosh = numeric(length(thetagrid))
  mGrunwaldRam = numeric(length(thetagrid))
  mGrunwaldCosh = numeric(length(thetagrid))
  mNormRam = numeric(length(thetagrid))
  mNormCosh = numeric(length(thetagrid))
  wRamdasRam = numeric(length(thetagrid))
  wRamdasCosh = numeric(length(thetagrid))
  wGrunwaldRam = numeric(length(thetagrid))
  wGrunwaldCosh = numeric(length(thetagrid))
  wNormRam = numeric(length(thetagrid))
  wNormCosh = numeric(length(thetagrid))
  RamdasRamrec = numeric(B)
  RamdasCoshrec = numeric(B)
  GrunwaldRamrec = numeric(B)
  GrunwaldCoshrec = numeric(B)
  NormRamrec = numeric(B)
  NormCoshrec = numeric(B)
}

for(j in 1:length(thetagrid)){
  for(b in 1:B){
    a = rlogalt(100000,thetagrid[j])
    c = rlogalt(100000)
    z = a-c
    RamdasRamrec[b] = RamdasTest(z,Eval = 'Ram')[[3]]
    RamdasCoshrec[b] = RamdasTest(z,Eval = 'NewE')[[3]]
    GrunwaldRamrec[b] = EfronsGExp(z,type = 'Ram',twosided = T)[[1]]
    GrunwaldCoshrec[b] = EfronsGExp(z,type = 'NewE',twosided = T)[[1]]
    NormRamrec[b] = EfronsScaleInvNorm(z,1,type = 'Ram',twosided = TRUE)[[1]]
    NormCoshrec[b] = EfronsScaleInvNorm(z,1,type = 'NewE',twosided = TRUE)[[1]]
    print(b)
  }
  wRamdasRam[j] = quantile(RamdasRamrec,0.8)
  wRamdasCosh[j] = quantile(RamdasCoshrec,0.8)
  wGrunwaldRam[j] = quantile(GrunwaldRamrec,0.8)
  wGrunwaldCosh[j] = quantile(GrunwaldCoshrec,0.8)
  wNormRam[j] = quantile(NormRamrec,0.8)
  wNormCosh[j] = quantile(NormCoshrec,0.8)
  RamdasRamrec[RamdasRamrec>quantile(RamdasRamrec,0.8)] = quantile(RamdasRamrec,0.8)
  mRamdasRam[j] = mean(RamdasRamrec)
  RamdasCoshrec[RamdasCoshrec>quantile(RamdasCoshrec,0.8)] = quantile(RamdasCoshrec,0.8)
  mRamdasCosh[j] = mean(RamdasCoshrec)
  GrunwaldRamrec[GrunwaldRamrec>quantile(GrunwaldRamrec,0.8)] = quantile(GrunwaldRamrec,0.8)
  mGrunwaldRam[j] = mean(GrunwaldRamrec)
  GrunwaldCoshrec[GrunwaldCoshrec>quantile(GrunwaldCoshrec,0.8)] = quantile(GrunwaldCoshrec,0.8)
  mGrunwaldCosh[j] = mean(GrunwaldCoshrec)
  NormRamrec[NormRamrec>quantile(NormRamrec,0.8)] = quantile(NormRamrec,0.8)
  mNormRam[j] = mean(NormRamrec)
  NormCoshrec[NormCoshrec>quantile(NormCoshrec,0.8)] = quantile(NormCoshrec,0.8)
  mNormCosh[j] = mean(NormCoshrec)
}
{
  plot(thetagrid,mGrunwaldCosh/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],xlim = c(0.5,3),ylim = c(0,8),col = 'blue',xlab = 'Location parameter' ,ylab = 'Comparative average stopping time',type = 'l',lwd = 2,main = 'X ~ Logistic alternative')
  lines(thetagrid,mGrunwaldRam/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],col = 'green',lwd =2)
  lines(thetagrid,mRamdasRam/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],col = 'black' ,lwd = 2)
  lines(thetagrid,mRamdasCosh/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],col = 'red',lwd = 2)
  lines(thetagrid,mNormRam/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],col = 'orange',lwd = 2)
  lines(thetagrid,mNormCosh/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)], col = 'purple',lwd = 2)
  legend(2.2,8,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 0.65)
}

LogAltMStop = recordPlot()

{
  plot(thetagrid,wGrunwaldCosh/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],xlim = c(0.5,3),ylim = c(0,11),col = 'blue',xlab = 'Location parameter',ylab = 'Comparative worst case stopping time',type = 'l',lwd = 2,main = 'X ~ Logistic alternative')
  lines(thetagrid,wGrunwaldRam/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],col = 'green',lwd =2)
  lines(thetagrid,wRamdasRam/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],col = 'black' ,lwd = 2)
  lines(thetagrid,wRamdasCosh/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],col = 'red',lwd = 2)
  lines(thetagrid,wNormRam/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)],col = 'orange',lwd = 2)
  lines(thetagrid,wNormCosh/regmwurecrecLogAlt[3:length(regmwurecrecLogAlt)], col = 'purple',lwd = 2)
  legend(2.2,11,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh'),text.col = c('black','red','green','blue','orange','purple'),text.width = 0.65)
}

LogAltWStop = recordPlot()
}

############# SECTION RANK BASED METHODS


### SPLIT MWU WHICH PRIOR?

# FIGURE 4.1
{
{
  library(invgamma)
  library(rootSolve)
  B = 1000
  N = 50
  Split1 = matrix(0,nrow = N,ncol = B)
  Split2 = matrix(0,nrow = N,ncol = B)
  Split3 = matrix(0,nrow = N,ncol = B)
}
for(b in 1:B){
  a = rnorm(N,0)
  c = rnorm(N,1)
  Split1t = SplitMWUDelta(a,c)[[1]]
  Split2t = SplitMWUDelta(a,c,prior = 'normal')[[1]]
  Split3t = SplitMWUDelta(a,c,prior = 'invnormal')[[1]]
  for(n in 1:N){
    Split1[n,b] = prod(Split1t[1:n])
    Split2[n,b] = prod(Split2t[1:n])
    Split3[n,b] = prod(Split3t[1:n])
  }
  print(b)
}
{
  cols = rainbow(4)
  plot(rowSums(log(Split1)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X_a ~ Normal(mean = 0,sd = 1), X_b ~ Normal(mean = delta, sd = 0)',ylim = c(0,6),lwd = 4)
  lines(rowSums(log(Split2)/B),col = 'black' ,lwd = 2)
  lines(rowSums(log(Split3)/B),col = 'red', lwd = 2)
  SPlitMWUGROWTH = recordPlot()
  legend(0,6,c('Split-MWU(prior = uniform)','Split-MWU(prior = normal)','Split-MWU(prior = 1/normal)'),text.col = c('blue','black','red'),text.width = 20)
}
}

# FIGURE 4.2
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    Split1 = matrix(0,nrow = N,ncol = B)
    Split2 = matrix(0,nrow = N,ncol = B)
    Split3 = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    a = rcauchy(N,0)
    c = rcauchy(N,2)
    Split1t = SplitMWUDelta(a,c)[[1]]
    Split2t = SplitMWUDelta(a,c,prior = 'normal')[[1]]
    Split3t = SplitMWUDelta(a,c,prior = 'invnormal')[[1]]
    for(n in 1:N){
      Split1[n,b] = prod(Split1t[1:n])
      Split2[n,b] = prod(Split2t[1:n])
      Split3[n,b] = prod(Split3t[1:n])
    }
    print(b)
  }
  {
    cols = rainbow(4)
    plot(rowSums(log(Split1)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X_a ~ Cauchy(location = 0), X_b ~ Cauchy(location = 2)',ylim = c(-0.5,6),lwd = 4)
    lines(rowSums(log(Split2)/B),col = 'black' ,lwd = 2)
    lines(rowSums(log(Split3)/B),col = 'red', lwd = 2)
    SPLITGROWTHCAUCH = recordPlot()
    legend(0,6,c('Split-MWU(prior = uniform)','Split-MWU(prior = normal)','Split-MWU(prior = 1/normal)'),text.col = c('blue','black','red'),text.width = 20)
  }
}

# FIGURE 4.3
{
  {
    library(invgamma)
    library(rootSolve)
    B = 1000
    N = 50
    Split1 = matrix(0,nrow = N,ncol = B)
    Split2 = matrix(0,nrow = N,ncol = B)
    Split3 = matrix(0,nrow = N,ncol = B)
  }
  for(b in 1:B){
    a = rlogalt(N,0)
    c = rlogalt(N,2)
    Split1t = SplitMWUDelta(a,c)[[1]]
    Split2t = SplitMWUDelta(a,c,prior = 'normal')[[1]]
    Split3t = SplitMWUDelta(a,c,prior = 'invnormal')[[1]]
    for(n in 1:N){
      Split1[n,b] = prod(Split1t[1:n])
      Split2[n,b] = prod(Split2t[1:n])
      Split3[n,b] = prod(Split3t[1:n])
    }
    print(b)
  }
  {
    cols = rainbow(4)
    plot(rowSums(log(Split1)/B),type = 'l',xlab = 'Sample size (N)',ylab = 'Expectation(log(E))', col = 'blue',main = 'X_a ~ LogAlt(delta = 0), X_b ~ LogAlt(delta = 1, sd = 0)',ylim = c(-0.5,8),lwd = 4)
    lines(rowSums(log(Split2)/B),col = 'black' ,lwd = 2)
    lines(rowSums(log(Split3)/B),col = 'red', lwd = 2)
    SPLITGROWTHLOGALT = recordPlot()
    legend(0,8,c('Split-MWU(prior = uniform)','Split-MWU(prior = normal)','Split-MWU(prior = 1/normal)'),text.col = c('blue','black','red'),text.width = 20)
  }
}

### SRT WHICH PRIOR?

# FIGURE 4.4
{
  B = 1000
  N = 100
  SRT1 = matrix(0,nrow = N,ncol = B)
  SRT2 = matrix(0,nrow = N,ncol = B)
  SRT3 = matrix(0,nrow = N,ncol = B)
  
  for(b in 1:B){
    a = rnorm(100,0)
    c = rnorm(100,1)
    data = cbind(c(rbind(a,c)),c(rbind(rep(1,100),rep(2,100))))
    SRT1t = MurielsPrior(data,prior = 'divided')[[1]]
    SRT2t = MurielsPrior(data,prior = 'normal')[[1]]
    SRT3t = MurielsPrior(data, prior = 'uniform')[[1]]
    for(n in 1:N){
      SRT1[n,b] = prod(SRT1t[1:n])
      SRT2[n,b] = prod(SRT2t[1:n])
      SRT3[n,b] = prod(SRT3t[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(SRT1)/B),col = 'black',lwd = 2,type = 'l',main = 'X_a ~ Normal(mean = 0, sd = 1), X_b ~ Normal(mean = 1, sd = 1)')
    lines(rowSums(log(SRT2)/B), col = 'blue',lwd = 2)
    lines(rowSums(log(SRT3)/B),col = 'green' ,lwd = 2)
    legend(0,6.4,c('SRT(prior = divided)','SRT(prior = normal)','SRT(prior = uniform)'),text.col = c('black','blue','green'),text.width = 25)
  }
  SRTCompare1 = recordPlot()
}

# FIGURE 4.5
{
  B = 1000
  N = 100
  SRT1 = matrix(0,nrow = N,ncol = B)
  SRT2 = matrix(0,nrow = N,ncol = B)
  SRT3 = matrix(0,nrow = N,ncol = B)
  
  for(b in 1:B){
    a = rcauchy(100,0)
    c = rcauchy(100,2)
    data = cbind(c(rbind(a,c)),c(rbind(rep(1,100),rep(2,100))))
    SRT1t = MurielsPrior(data,prior = 'divided')[[1]]
    SRT2t = MurielsPrior(data,prior = 'normal')[[1]]
    SRT3t = MurielsPrior(data, prior = 'uniform')[[1]]
    for(n in 1:N){
      SRT1[n,b] = prod(SRT1t[1:n])
      SRT2[n,b] = prod(SRT2t[1:n])
      SRT3[n,b] = prod(SRT3t[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(SRT1)/B),col = 'black',lwd = 2,type = 'l',main = 'X_a ~ Cauchy(location = 0), X_b ~ Cauchy(location = 2)')
    lines(rowSums(log(SRT2)/B), col = 'blue',lwd = 2)
    lines(rowSums(log(SRT3)/B),col = 'green' ,lwd = 2)
    legend(0,3,c('SRT(prior = divided)','SRT(prior = normal)','SRT(prior = uniform)'),text.col = c('black','blue','green'),text.width = 25)
  }
  SRTCompareCauchy = recordPlot()
}

# FIGURE 4.6
{
  B = 1000
  N = 100
  SRT1 = matrix(0,nrow = N,ncol = B)
  SRT2 = matrix(0,nrow = N,ncol = B)
  SRT3 = matrix(0,nrow = N,ncol = B)
  
  for(b in 1:B){
    a = rlogalt(100,0)
    c = rlogalt(100,2)
    data = cbind(c(rbind(a,c)),c(rbind(rep(1,100),rep(2,100))))
    SRT1t = MurielsPrior(data,prior = 'divided')[[1]]
    SRT2t = MurielsPrior(data,prior = 'normal')[[1]]
    SRT3t = MurielsPrior(data, prior = 'uniform')[[1]]
    for(n in 1:N){
      SRT1[n,b] = prod(SRT1t[1:n])
      SRT2[n,b] = prod(SRT2t[1:n])
      SRT3[n,b] = prod(SRT3t[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(SRT1)/B),col = 'black',lwd = 2,type = 'l',main = 'X_a ~ LogAlt(delta = 0), X_b ~ LogAlt(delta = 2)')
    lines(rowSums(log(SRT2)/B), col = 'blue',lwd = 2)
    lines(rowSums(log(SRT3)/B),col = 'green' ,lwd = 2)
    legend(0,8.5,c('SRT(prior = divided)','SRT(prior = normal)','SRT(prior = uniform)'),text.col = c('black','blue','green'),text.width = 25)
  }
  SRTCompareLogAlt = recordPlot()
}

######### FINAL COMPARISON: MURIELS VS SPLIT-MWU VS EFRON VS MWU VS SMWU
#Comparison for power
# FIGURE 5.1 & 5.5 & 5.9
{
  {
    thetagrid = seq(0.4,3,0.1)
    library(invgamma)
    library(rootSolve)
    B = 500
    mMuriel1 = numeric(length(thetagrid))
    wMuriel1 = numeric(length(thetagrid))
    mSplit = numeric(length(thetagrid))
    wSplit = numeric(length(thetagrid))
    mGrunCosh = numeric(length(thetagrid))
    wGrunCosh = numeric(length(thetagrid))
    mGrunCoshw = numeric(length(thetagrid))
    wGrunCoshw = numeric(length(thetagrid))
    mNormRam = numeric(length(thetagrid))
    wNormRam = numeric(length(thetagrid))
    Muriel1rec = numeric(B)
    Splitrec = numeric(B)
    GrunCoshrec = numeric(B)
    GrunCoshrecw = numeric(B)
    NormRamrec = numeric(B)
  }
  
  for(j in 1:length(thetagrid)){
    for(b in 1:B){
      a = rnorm(150,0)
      c = rnorm(150,thetagrid[j])
      data = cbind(c(rbind(a,c)),c(rbind(rep(1,150),rep(2,150))))
      Muriel1rec[b] = MurielsPriorStop(data,prior = 'normal')[[3]]
      Splitrec[b] = SplitMWUDeltaStop(a,c,prior = 'normal')[[2]]
      GrunCoshrec[b] = EfronsGExp(a-c,type = 'NewE',twosided = T)[[1]]
      GrunCoshrecw[b] = EfronsGExpW(a-c,type = 'NewE',twosided = T,W=0.1)[[1]]
      NormRamrec[b] = EfronsScaleInvNorm(a-c,type = 'Ram',twosided = T)[[1]]
      print(b)
    }
    wMuriel1[j] = quantile(Muriel1rec,0.8)
    Muriel1rec[Muriel1rec>quantile(Muriel1rec,0.8)] = quantile(Muriel1rec,0.8)
    mMuriel1[j] = mean(Muriel1rec)
    wSplit[j] = quantile(Splitrec,0.8)
    Splitrec[Splitrec>quantile(Splitrec,0.8)] = quantile(Splitrec,0.8)
    mSplit[j] = mean(Splitrec)
    wGrunCosh[j] = quantile(GrunCoshrec,0.8)
    GrunCoshrec[GrunCoshrec>quantile(GrunCoshrec,0.8)] = quantile(GrunCoshrec,0.8)
    mGrunCosh[j] = mean(GrunCoshrec)
    wGrunCoshw[j] = quantile(GrunCoshrecw,0.8)
    GrunCoshrecw[GrunCoshrecw>quantile(GrunCoshrecw,0.8)] = quantile(GrunCoshrecw,0.8)
    mGrunCoshw[j] = mean(GrunCoshrecw)
    wNormRam[j] = quantile(NormRamrec,0.8)
    NormRamrec[NormRamrec>quantile(NormRamrec,0.8)] = quantile(NormRamrec,0.8)
    mNormRam[j] = mean(NormRamrec)
    print(thetagrid[j])
  }
  
  {
    plot(thetagrid,mMuriel1,xlim = c(0.25,3),ylim = c(5,100),col = 'orange',xlab = 'Delta',ylab = 'Sample size needed for power 0.8 / Average stopping time n',type = 'l',lwd = 3,main = 'X_a ~ Normal(mean = 0), X_b ~ Normal(mean = delta)')
    lines(thetagrid,mSplit*2,col = 'blue',lwd = 2)
    lines(thetagrid,mGrunCosh*2,col = 'red',lwd = 2)
    lines(thetagrid,mGrunCoshw*2,col = 'brown',lwd = 2)
    lines(thetagrid,mNormRam*2,col = 'purple',lwd = 2)
    lines(thetagrid, regmwurecrec[seq(3,length(regmwurecrec))]*2,col = 'black',lwd = 2)
#    lines(thetagrid,safmwurecrec*2, col = 'purple',lwd = 2)
    FinalNormAvg = recordPlot()
    legend(1.8,100,c('Muriels(Prior = 1/theta)','Split-MWU','Gr.E_Cosh','Gr.E_Cosh(h = 0.1)','Norm.E_Ram','MWU'),text.col = c('orange','blue','red','brown','purple','black'),text.width = 1)
  }
  {
    plot(thetagrid,wMuriel1,xlim = c(0.25,3),ylim = c(5,100),col = 'orange',xlab = 'Delta',ylab = 'Sample size needed for power 0.8 / Worst case stopping time n',type = 'l',lwd = 3,main = 'X_a ~ Normal(mean = 0), X_b ~ Normal(mean = delta)')
    lines(thetagrid,wSplit*2,col = 'blue',lwd = 2)
    lines(thetagrid,wGrunCosh*2,col = 'red',lwd = 2)
    lines(thetagrid,wGrunCoshw*2,col = 'brown',lwd = 2)
    lines(thetagrid,wNormRam*2,col = 'purple',lwd = 2)
    lines(thetagrid, regmwurecrec[seq(3,length(regmwurecrec))]*2,col = 'black',lwd = 2)
    #    lines(thetagrid,safmwurecrec*2, col = 'purple',lwd = 2)
    FinalNormWorst = recordPlot()
    legend(1.8,100,c('Muriels(Prior = 1/theta)','Split-MWU','Gr.E_Cosh','Gr.E_Cosh(h = 0.1)','Norm.E_Ram','MWU'),text.col = c('orange','blue','red','brown','purple','black'),text.width = 1)
  }
}

# FIGURE 5.2 & 5.6 & 5.10
{
  {
    thetagrid = seq(0.4,3,0.1)
    library(invgamma)
    library(rootSolve)
    B = 500
    mMuriel1 = numeric(length(thetagrid))
    wMuriel1 = numeric(length(thetagrid))
    mSplit = numeric(length(thetagrid))
    wSplit = numeric(length(thetagrid))
    mGrunCosh = numeric(length(thetagrid))
    wGrunCosh = numeric(length(thetagrid))
    mGrunCoshw = numeric(length(thetagrid))
    wGrunCoshw = numeric(length(thetagrid))
    mNormRam = numeric(length(thetagrid))
    wNormRam = numeric(length(thetagrid))
    Muriel1rec = numeric(B)
    Splitrec = numeric(B)
    GrunCoshrec = numeric(B)
    NormRamrec = numeric(B)
  }
  
  for(j in 1:length(thetagrid)){
    for(b in 1:B){
      a = rt(150,3,0)
      c = rt(150,3,thetagrid[j])
      data = cbind(c(rbind(a,c)),c(rbind(rep(1,150),rep(2,150))))
      Muriel1rec[b] = MurielsPriorStop(data,prior = 'normal')[[3]]
      Splitrec[b] = SplitMWUDeltaStop(a,c,prior = 'normal')[[2]]
      GrunCoshrec[b] = EfronsGExp(a-c,type = 'NewE',twosided = T)[[1]]
      GrunCoshrecw[b] = EfronsGExpW(a-c,type = 'NewE',twosided = T,W=0.1)[[1]]
      NormRamrec[b] = EfronsScaleInvNorm(a-c,type = 'Ram',twosided = T)[[1]]
      print(b)
    }
    wMuriel1[j] = quantile(Muriel1rec,0.8)
    Muriel1rec[Muriel1rec>quantile(Muriel1rec,0.8)] = quantile(Muriel1rec,0.8)
    mMuriel1[j] = mean(Muriel1rec)
    wSplit[j] = quantile(Splitrec,0.8)
    Splitrec[Splitrec>quantile(Splitrec,0.8)] = quantile(Splitrec,0.8)
    mSplit[j] = mean(Splitrec)
    wGrunCosh[j] = quantile(GrunCoshrec,0.8)
    GrunCoshrec[GrunCoshrec>quantile(GrunCoshrec,0.8)] = quantile(GrunCoshrec,0.8)
    mGrunCosh[j] = mean(GrunCoshrec)
    wGrunCoshw[j] = quantile(GrunCoshrecw,0.8)
    GrunCoshrecw[GrunCoshrecw>quantile(GrunCoshrecw,0.8)] = quantile(GrunCoshrecw,0.8)
    mGrunCoshw[j] = mean(GrunCoshrecw)
    wNormRam[j] = quantile(NormRamrec,0.8)
    NormRamrec[NormRamrec>quantile(NormRamrec,0.8)] = quantile(NormRamrec,0.8)
    mNormRam[j] = mean(NormRamrec)
    print(thetagrid[j])
  }

  {
    plot(thetagrid,mMuriel1,xlim = c(0.25,3),ylim = c(5,150),col = 'orange',xlab = 'Delta',ylab = 'Sample size needed for power 0.8 / Average stopping time n',type = 'l',lwd = 3,main = 'X_a ~ Student`s t (ncp = 0), X_b ~ Student`s t (ncp = Delta)')
    lines(thetagrid,mSplit*2,col = 'blue',lwd = 2)
    lines(thetagrid,mGrunCosh*2,col = 'red',lwd = 2)
    lines(thetagrid,mGrunCoshw*2,col = 'brown',lwd = 2)
    lines(thetagrid,mNormRam*2,col = 'purple',lwd = 2)
    lines(thetagrid, regmwurecrecT[seq(3,length(regmwurecrecT))]*2,col = 'black',lwd = 2)
    #    lines(thetagrid,safmwurecrec*2, col = 'purple',lwd = 2)
    FinalTAvg = recordPlot()
    legend(1.8,150,c('SRT','Split-MWU','Gr.E_Cosh','Gr.E_Cosh(h = 0.1)','Norm.E_Ram','MWU'),text.col = c('orange','blue','red','brown','purple','black'),text.width = 1)
  }
  {
    plot(thetagrid,wMuriel1,xlim = c(0.25,3),ylim = c(5,150),col = 'orange',xlab = 'Delta',ylab = 'Sample size needed for power 0.8 / Worst case stopping time n',type = 'l',lwd = 3,main = 'X_a ~ Student`s t (ncp = 0), X_b ~ Student`s t (ncp = Delta)')
    lines(thetagrid,wSplit*2,col = 'blue',lwd = 2)
    lines(thetagrid,wGrunCosh*2,col = 'red',lwd = 2)
    lines(thetagrid,wGrunCoshw*2,col = 'brown',lwd = 2)
    lines(thetagrid,wNormRam*2,col = 'purple',lwd = 2)
    lines(thetagrid, regmwurecrecT[seq(3,length(regmwurecrecT))]*2,col = 'black',lwd = 2)
    #    lines(thetagrid,safmwurecrec*2, col = 'purple',lwd = 2)
    FinalTWorst = recordPlot()
    legend(1.8,150,c('SRT','Split-MWU','Gr.E_Cosh','Gr.E_Cosh(h = 0.1)','Norm.E_Ram','MWU'),text.col = c('orange','blue','red','brown','purple','black'),text.width = 1)
  }
  
}

# FIGURE 5.3 & 5.7
{
  {
    thetagrid = seq(0.4,3,0.1)
    library(invgamma)
    library(rootSolve)
    B = 500
    mMuriel1 = numeric(length(thetagrid))
    wMuriel1 = numeric(length(thetagrid))
    mSplit = numeric(length(thetagrid))
    wSplit = numeric(length(thetagrid))
    mGrunCosh = numeric(length(thetagrid))
    wGrunCosh = numeric(length(thetagrid))
    mGrunCoshw = numeric(length(thetagrid))
    wGrunCoshw = numeric(length(thetagrid))
    mNormRam = numeric(length(thetagrid))
    wNormRam = numeric(length(thetagrid))
    Muriel1rec = numeric(B)
    Splitrec = numeric(B)
    GrunCoshrec = numeric(B)
    NormRamrec = numeric(B)
  }
  
  for(j in 1:length(thetagrid)){
    for(b in 1:B){
      a = rcauchy(150,0)
      c = rcauchy(150,thetagrid[j])
      data = cbind(c(rbind(a,c)),c(rbind(rep(1,150),rep(2,150))))
      Muriel1rec[b] = MurielsPriorStop(data,prior = 'normal')[[3]]
      Splitrec[b] = SplitMWUDeltaStop(a,c,prior = 'normal')[[2]]
      GrunCoshrec[b] = EfronsGExp(a-c,type = 'NewE',twosided = T)[[1]]
      GrunCoshrecw[b] = EfronsGExpW(a-c,type = 'NewE',twosided = T,W=0.1)[[1]]
      NormRamrec[b] = EfronsScaleInvNorm(a-c,type = 'Ram',twosided = T)[[1]]
      print(b)
    }
    wMuriel1[j] = quantile(Muriel1rec,0.8)
    Muriel1rec[Muriel1rec>quantile(Muriel1rec,0.8)] = quantile(Muriel1rec,0.8)
    mMuriel1[j] = mean(Muriel1rec)
    wSplit[j] = quantile(Splitrec,0.8)
    Splitrec[Splitrec>quantile(Splitrec,0.8)] = quantile(Splitrec,0.8)
    mSplit[j] = mean(Splitrec)
    wGrunCosh[j] = quantile(GrunCoshrec,0.8)
    GrunCoshrec[GrunCoshrec>quantile(GrunCoshrec,0.8)] = quantile(GrunCoshrec,0.8)
    mGrunCosh[j] = mean(GrunCoshrec)
    wGrunCoshw[j] = quantile(GrunCoshrecw,0.8)
    GrunCoshrecw[GrunCoshrecw>quantile(GrunCoshrecw,0.8)] = quantile(GrunCoshrecw,0.8)
    mGrunCoshw[j] = mean(GrunCoshrecw)
    wNormRam[j] = quantile(NormRamrec,0.8)
    NormRamrec[NormRamrec>quantile(NormRamrec,0.8)] = quantile(NormRamrec,0.8)
    mNormRam[j] = mean(NormRamrec)
    print(thetagrid[j])
  }
  {
    plot(thetagrid,mMuriel1,xlim = c(0.0,3),ylim = c(5,160),col = 'orange',xlab = 'Delta',ylab = 'Sample size needed for power 0.8 / Average stopping time n',type = 'l',lwd = 3,main = 'X_a ~ Cauchy (location = 0), X_b ~ Cauchy (location = Delta)')
    lines(thetagrid,mSplit*2,col = 'blue',lwd = 2)
    lines(thetagrid,mGrunCosh*2,col = 'red',lwd = 2)
    lines(thetagrid,mGrunCoshw*2,col = 'brown',lwd = 2)
    lines(thetagrid,mNormRam*2,col = 'purple',lwd = 2)
    lines(thetagrid, regmwurecrecCauch[seq(3,length(regmwurecrecCauch))]*2,col = 'black',lwd = 2)
    #    lines(thetagrid,safmwurecrec*2, col = 'purple',lwd = 2)
    FinalCauchAvg = recordPlot()
    legend(0,130,c('Muriels(Prior = 1/theta)','Split-MWU','Gr.E_Cosh','Gr.E_Cosh(h = 0.1)','Norm.E_Ram','MWU'),text.col = c('orange','blue','red','brown','purple','black'),text.width = 1)
  }
  {
    plot(thetagrid,wMuriel1,xlim = c(0,3),ylim = c(5,150),col = 'orange',xlab = 'Delta',ylab = 'Sample size needed for power 0.8 / Worst case stopping time n',type = 'l',lwd = 3,main = 'X_a ~ Cauchy (location = 0), X_b ~ Cauchy (location = Delta)')
    lines(thetagrid,wSplit*2,col = 'blue',lwd = 2)
    lines(thetagrid,wGrunCosh*2,col = 'red',lwd = 2)
    lines(thetagrid,wGrunCoshw*2,col = 'brown',lwd = 2)
    lines(thetagrid,wNormRam*2,col = 'purple',lwd = 2)
    lines(thetagrid, regmwurecrecCauch[seq(3,length(regmwurecrecCauch))]*2,col = 'black',lwd = 2)
    #    lines(thetagrid,safmwurecrec*2, col = 'purple',lwd = 2)
    FinalCauchWorst = recordPlot()
    legend(0,130,c('Muriels(Prior = 1/theta)','Split-MWU','Gr.E_Cosh','Gr.E_Cosh(h = 0.1)','Norm.E_Ram','MWU'),text.col = c('orange','blue','red','brown','purple','black'),text.width = 1)
  }
}

# FIGURE 5.4 & 5.8
{
  {
    thetagrid = seq(0.4,3,0.1)
    library(invgamma)
    library(rootSolve)
    B = 500
    mMuriel1 = numeric(length(thetagrid))
    wMuriel1 = numeric(length(thetagrid))
    mSplit = numeric(length(thetagrid))
    wSplit = numeric(length(thetagrid))
    mGrunCosh = numeric(length(thetagrid))
    wGrunCosh = numeric(length(thetagrid))
    mGrunCoshw = numeric(length(thetagrid))
    wGrunCoshw = numeric(length(thetagrid))
    mNormRam = numeric(length(thetagrid))
    wNormRam = numeric(length(thetagrid))
    Muriel1rec = numeric(B)
    Splitrec = numeric(B)
    GrunCoshrec = numeric(B)
    NormRamrec = numeric(B)
  }
  
  for(j in 1:length(thetagrid)){
    for(b in 1:B){
      a = rlogalt(150,0)
      c = rlogalt(150,thetagrid[j])
      data = cbind(c(rbind(a,c)),c(rbind(rep(1,150),rep(2,150))))
      Muriel1rec[b] = MurielsPriorStop(data,prior = 'normal')[[3]]
      Splitrec[b] = SplitMWUDeltaStop(a,c,prior = 'normal')[[2]]
      GrunCoshrec[b] = EfronsGExp(a-c,type = 'NewE',twosided = T)[[1]]
      GrunCoshrecw[b] = EfronsGExpW(a-c,type = 'NewE',twosided = T,W=0.1)[[1]]
      NormRamrec[b] = EfronsScaleInvNorm(a-c,type = 'Ram',twosided = T)[[1]]
      print(b)
    }
    wMuriel1[j] = quantile(Muriel1rec,0.8)
    Muriel1rec[Muriel1rec>quantile(Muriel1rec,0.8)] = quantile(Muriel1rec,0.8)
    mMuriel1[j] = mean(Muriel1rec)
    wSplit[j] = quantile(Splitrec,0.8)
    Splitrec[Splitrec>quantile(Splitrec,0.8)] = quantile(Splitrec,0.8)
    mSplit[j] = mean(Splitrec)
    wGrunCosh[j] = quantile(GrunCoshrec,0.8)
    GrunCoshrec[GrunCoshrec>quantile(GrunCoshrec,0.8)] = quantile(GrunCoshrec,0.8)
    mGrunCosh[j] = mean(GrunCoshrec)
    wGrunCoshw[j] = quantile(GrunCoshrecw,0.8)
    GrunCoshrecw[GrunCoshrecw>quantile(GrunCoshrecw,0.8)] = quantile(GrunCoshrecw,0.8)
    mGrunCoshw[j] = mean(GrunCoshrecw)
    wNormRam[j] = quantile(NormRamrec,0.8)
    NormRamrec[NormRamrec>quantile(NormRamrec,0.8)] = quantile(NormRamrec,0.8)
    mNormRam[j] = mean(NormRamrec)
    print(thetagrid[j])
  }
  {
    plot(thetagrid,mMuriel1,xlim = c(0.0,3),ylim = c(5,160),col = 'orange',xlab = 'Delta',ylab = 'Sample size needed for power 0.8 / Average stopping time n',type = 'l',lwd = 3,main = 'X_a ~ LogAlt (location = 0), X_b ~ LogAlt (location = Delta)')
    lines(thetagrid,mSplit*2,col = 'blue',lwd = 2)
    lines(thetagrid,mGrunCosh*2,col = 'red',lwd = 2)
    lines(thetagrid,mGrunCoshw*2,col = 'brown',lwd = 2)
    lines(thetagrid,mNormRam*2,col = 'purple',lwd = 2)
    lines(thetagrid, regmwurecrecLogAlt[seq(3,length(regmwurecrecLogAlt))]*2,col = 'black',lwd = 2)
    #    lines(thetagrid,safmwurecrec*2, col = 'purple',lwd = 2)
    FinalLogAltAvg = recordPlot()
    legend(1.8,150,c('SRT','Split-MWU','Gr.E_Cosh','Gr.E_Cosh(h = 0.1)','Norm.E_Ram','MWU'),text.col = c('orange','blue','red','brown','purple','black'),text.width = 1)
  }
  {
    plot(thetagrid,wMuriel1,xlim = c(0,3),ylim = c(5,150),col = 'orange',xlab = 'Delta',ylab = 'Sample size needed for power 0.8 / Worst case stopping time n',type = 'l',lwd = 3,main = 'X_a ~ LogAlt (location = 0), X_b ~ LogAlt (location = Delta)')
    lines(thetagrid,wSplit*2,col = 'blue',lwd = 2)
    lines(thetagrid,wGrunCosh*2,col = 'red',lwd = 2)
    lines(thetagrid,wGrunCoshw*2,col = 'brown',lwd = 2)
    lines(thetagrid,wNormRam*2,col = 'purple',lwd = 2)
    lines(thetagrid, regmwurecrecLogAlt[seq(3,length(regmwurecrecLogAlt))]*2,col = 'black',lwd = 2)
    #    lines(thetagrid,safmwurecrec*2, col = 'purple',lwd = 2)
    FinalLogAltWorst = recordPlot()
    legend(1.8,150,c('SRT','Split-MWU','Gr.E_Cosh','Gr.E_Cosh(h = 0.1)','Norm.E_Ram','MWU'),text.col = c('orange','blue','red','brown','purple','black'),text.width = 1)
  }
}

#Comparison for growth

# FIGURE 5.11
{
  B = 1000
  N = 200
  MurielPrior = matrix(0,nrow = N,ncol = B)
  Split = matrix(0,nrow = N/2,ncol = B)
  RamRam = matrix(0,nrow = N/2,ncol = B)
  RamCosh = matrix(0,nrow = N/2,ncol = B)
  GrunwaldRam = matrix(0,nrow = N/2,ncol = B)
  GrunwaldCosh = matrix(0,nrow = N/2,ncol = B)
  wGrunwaldCosh = matrix(0,nrow = N/2,ncol = B)
  NormRam = matrix(0,nrow = N/2,ncol = B)
  NormCosh = matrix(0,nrow = N/2,ncol = B)
  for(b in 1:B){
    a = rnorm(100,0)
    c = rnorm(100,1)
    RamRamt = RamdasTestlog(a-c, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(a-c, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(a-c,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(a-c,type = 'NewE',twosided = T)[[6]]
    wGrunwaldCosht = EfronsGExpWlog(a-c,type = 'NewE',twosided = T,W = 0.1)[[6]]
    NormRamt = EfronsScaleInvNormlog(a-c,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(a-c,1,type = 'NewE',twosided = TRUE)[[6]]
    data = cbind(c(rbind(a,c)),c(rbind(rep(1,100),rep(2,100))))
    MurielPriort = MurielsPrior(data,prior = 'normal')[[1]]
    Split1t = SplitMWUDelta(a,c,prior = 'normal')[[1]]
    for(n in 1:N){
      MurielPrior[n,b] = prod(MurielPriort[1:n])
    }
    for(n in 1:N/2){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      wGrunwaldCosh[n,b] = prod(wGrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
      Split[n,b] = prod(Split1t[1:n])
    }
    print(b)
  }
  
  {
    plot(rowSums(log(MurielPrior)/B),col = 'brown',lwd = 2,type = 'l',main = 'X ~ Normal(mean = 1, sd = 1)',xlab = 'Total sample size', ylab = 'Expectation[log(E)]' )
    lines(seq(2,200,2),rowSums(log(Split)/B),lwd = 2,col = 'gray')
    lines(seq(2,200,2),rowSums(log(GrunwaldCosh)/B), col = 'blue',lwd = 2)
    lines(seq(2,200,2),rowSums(log(wGrunwaldCosh)/B), col = 'blue',lwd = 2,type ='o')
    lines(seq(2,200,2),rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
    legend(0,15,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh','SRT','Split-MWU'),text.col = c('black','red','green','blue','orange','purple','brown','gray'),text.width = 40)
  }
  FinalGrowthNorm = recordPlot()
}

# FIGURE 5.12
{
  B = 1000
  N = 200
  MurielPrior = matrix(0,nrow = N,ncol = B)
  Split = matrix(0,nrow = N/2,ncol = B)
  RamRam = matrix(0,nrow = N/2,ncol = B)
  RamCosh = matrix(0,nrow = N/2,ncol = B)
  GrunwaldRam = matrix(0,nrow = N/2,ncol = B)
  GrunwaldCosh = matrix(0,nrow = N/2,ncol = B)
  wGrunwaldCosh = matrix(0,nrow = N/2,ncol = B)
  NormRam = matrix(0,nrow = N/2,ncol = B)
  NormCosh = matrix(0,nrow = N/2,ncol = B)
  for(b in 1:B){
    a = rt(100,3,0)
    c = rt(100,3,1)
    RamRamt = RamdasTestlog(a-c, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(a-c, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(a-c,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(a-c,type = 'NewE',twosided = T)[[6]]
    wGrunwaldCosht = EfronsGExpWlog(a-c,type = 'NewE',twosided = T,W = 0.1)[[6]]
    NormRamt = EfronsScaleInvNormlog(a-c,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(a-c,1,type = 'NewE',twosided = TRUE)[[6]]
    data = cbind(c(rbind(a,c)),c(rbind(rep(1,100),rep(2,100))))
    MurielPriort = MurielsPrior(data,prior = 'normal')[[1]]
    Split1t = SplitMWUDelta(a,c,prior = 'normal')[[1]]
    for(n in 1:N){
      MurielPrior[n,b] = prod(MurielPriort[1:n])
    }
    for(n in 1:N/2){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      wGrunwaldCosh[n,b] = prod(wGrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
      Split[n,b] = prod(Split1t[1:n])
    }
    print(b)
  }
  {
    plot(rowSums(log(MurielPrior)/B),col = 'brown',ylim = c(-2,14),lwd = 2,type = 'l',xlab = 'Total sample size',ylab = 'Expectation[log(E)]',main = 'X ~ Student`s t (ncp = 1, df = 3)')
    lines(seq(2,200,2),rowSums(log(Split)/B),lwd = 2,col = 'gray')
    lines(seq(2,200,2),rowSums(log(GrunwaldCosh)/B), col = 'blue',lwd = 2)
    lines(seq(2,200,2),rowSums(log(wGrunwaldCosh)/B), col = 'blue',lwd = 2,type ='b')
    lines(seq(2,200,2),rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
    legend(0,14,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh','SRT','Split-MWU'),text.col = c('black','red','green','blue','orange','purple','brown','gray'),text.width = 40)
  }
  FinalGrowthT = recordPlot()
}

# FIGURE 5.13
{
  B = 1000
  N = 200
  MurielPrior = matrix(0,nrow = N,ncol = B)
  Split = matrix(0,nrow = N/2,ncol = B)
  RamRam = matrix(0,nrow = N/2,ncol = B)
  RamCosh = matrix(0,nrow = N/2,ncol = B)
  GrunwaldRam = matrix(0,nrow = N/2,ncol = B)
  GrunwaldCosh = matrix(0,nrow = N/2,ncol = B)
  wGrunwaldCosh = matrix(0,nrow = N/2,ncol = B)
  NormRam = matrix(0,nrow = N/2,ncol = B)
  NormCosh = matrix(0,nrow = N/2,ncol = B)
  for(b in 1:B){
    a = rcauchy(100,0)
    c = rcauchy(100,2)
    RamRamt = RamdasTestlog(a-c, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(a-c, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(a-c,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(a-c,type = 'NewE',twosided = T)[[6]]
    wGrunwaldCosht = EfronsGExpWlog(a-c,type = 'NewE',twosided = T,W = 0.1)[[6]]
    NormRamt = EfronsScaleInvNormlog(a-c,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(a-c,1,type = 'NewE',twosided = TRUE)[[6]]
    data = cbind(c(rbind(a,c)),c(rbind(rep(1,100),rep(2,100))))
    MurielPriort = MurielsPrior(data,prior = 'normal')[[1]]
    Split1t = SplitMWUDelta(a,c,prior = 'normal')[[1]]
    for(n in 1:N){
      MurielPrior[n,b] = prod(MurielPriort[1:n])
    }
    for(n in 1:N/2){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      wGrunwaldCosh[n,b] = prod(wGrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
      Split[n,b] = prod(Split1t[1:n])
    }
    print(b)
  }
  
  {
    plot(rowSums(log(MurielPrior)/B),col = 'brown',ylim = c(-2,14),lwd = 2,type = 'l',xlab = 'Total sample size',ylab = 'Expectation[log(E)]',main = 'X ~ Cauchy (location = 2, dispersion = 1)')
    lines(seq(2,200,2),rowSums(log(Split)/B),lwd = 2,col = 'gray')
    lines(seq(2,200,2),rowSums(log(GrunwaldCosh)/B), col = 'blue',lwd = 2)
    lines(seq(2,200,2),rowSums(log(wGrunwaldCosh)/B), col = 'blue',lwd = 2,type ='b')
    lines(seq(2,200,2),rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
    legend(0,13,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh','SRT','Split-MWU'),text.col = c('black','red','green','blue','orange','purple','brown','gray'),text.width = 40)
  }
  FinalGrowthCauchy = recordPlot()
}

# FIGURE 5.14
{
  B = 1000
  N = 200
  MurielPrior = matrix(0,nrow = N,ncol = B)
  Split = matrix(0,nrow = N/2,ncol = B)
  RamRam = matrix(0,nrow = N/2,ncol = B)
  RamCosh = matrix(0,nrow = N/2,ncol = B)
  GrunwaldRam = matrix(0,nrow = N/2,ncol = B)
  GrunwaldCosh = matrix(0,nrow = N/2,ncol = B)
  wGrunwaldCosh = matrix(0,nrow = N/2,ncol = B)
  NormRam = matrix(0,nrow = N/2,ncol = B)
  NormCosh = matrix(0,nrow = N/2,ncol = B)
  for(b in 1:B){
    a = rlogalt(100,0)
    c = rlogalt(100,2)
    RamRamt = RamdasTestlog(a-c, Eval = 'Ram')[[1]]
    RamCosht = RamdasTestlog(a-c, Eval = 'NewE')[[1]]
    GrunwaldRamt = EfronsGExplog(a-c,type = 'Ram',twosided = T)[[6]]
    GrunwaldCosht = EfronsGExplog(a-c,type = 'NewE',twosided = T)[[6]]
    wGrunwaldCosht = EfronsGExpWlog(a-c,type = 'NewE',twosided = T,W = 0.1)[[6]]
    NormRamt = EfronsScaleInvNormlog(a-c,1,type = 'Ram',twosided = TRUE)[[6]]
    NormCosht = EfronsScaleInvNormlog(a-c,1,type = 'NewE',twosided = TRUE)[[6]]
    data = cbind(c(rbind(a,c)),c(rbind(rep(1,100),rep(2,100))))
    MurielPriort = MurielsPrior(data,prior = 'normal')[[1]]
    Split1t = SplitMWUDelta(a,c,prior = 'normal')[[1]]
    for(n in 1:N){
      MurielPrior[n,b] = prod(MurielPriort[1:n])
    }
    for(n in 1:N/2){
      RamRam[n,b] = prod(RamRamt[1:n])
      RamCosh[n,b] = prod(RamCosht[1:n])
      GrunwaldRam[n,b] = prod(GrunwaldRamt[1:n])
      GrunwaldCosh[n,b] = prod(GrunwaldCosht[1:n])
      wGrunwaldCosh[n,b] = prod(wGrunwaldCosht[1:n])
      NormRam[n,b] = prod(NormRamt[1:n])
      NormCosh[n,b] = prod(NormCosht[1:n])
      Split[n,b] = prod(Split1t[1:n])
    }
    print(b)
  }
  
  {
    plot(rowSums(log(MurielPrior)/B),col = 'brown',ylim = c(-2,20),lwd = 2,type = 'l',xlab = 'Sample size n',ylab = 'Expectation[log(E)]',main = 'X ~ Logistic Alternative (delta = 2)')
    lines(seq(2,200,2),rowSums(log(Split)/B),lwd = 2,col = 'gray')
    lines(seq(2,200,2),rowSums(log(GrunwaldCosh)/B), col = 'blue',lwd = 2)
    lines(seq(2,200,2),rowSums(log(wGrunwaldCosh)/B), col = 'blue',lwd = 2,type ='b')
    lines(seq(2,200,2),rowSums(log(GrunwaldRam)/B),col = 'green' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(RamRam)/B),col = 'black' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(RamCosh)/B),col = 'red' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(NormRam)/B),col = 'orange' ,lwd = 2)
    lines(seq(2,200,2),rowSums(log(NormCosh)/B),col = 'purple' ,lwd = 2)
    legend(0,20,c('Ram.E_Ram','Ram.E_Cosh','Gr.E_Ram.','Gr.E_Cosh','Norm.E_Ram','Norm.E_Cosh','SRT','Split-MWU'),text.col = c('black','red','green','blue','orange','purple','brown','gray'),text.width = 40)
  }
  FinalGrowthLogAlt = recordPlot()
}

# [1] from https://stackoverflow.com/questions/57785222/avoiding-overflow-in-logcoshx
