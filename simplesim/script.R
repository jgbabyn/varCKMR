library(TMB)
  compile("./simplesim/models/CMP.cpp")
  dyn.load(dynlib("./simplesim/models/CMP"))
  compile("./simplesim/models/CMPSim.cpp")
  dyn.load(dynlib("./simplesim/models/CMPSim"))
compile("./simplesim/models/mprobsim.cpp")
dyn.load(dynlib("./simplesim/models/mprobsim"))
compile("./simplesim/models/mprobsimAdj.cpp")
dyn.load(dynlib("./simplesim/models/mprobsimAdj"))


library(extraDistr)

  data = list()
  data$n = 5000

  parm = list()
  parm$mean = 5
  parm$nu = 1

cmp = MakeADFun(data,parm,type="Fun",DLL="CMPSim")
cmpNS = MakeADFun(data,parm,type="Fun",DLL="CMP")

  ball_sim_diffCMP <- function(nn,n1,n2,meanV,nu){
    onesim <- replicate(nn,{
      balls1 = cmp$simulate(c(meanV,nu))$sims
      balls2 = cmp$simulate(c(meanV,nu))$sims
      draw1 = rmvhyper(1,balls1,n1)
      draw2 = rmvhyper(1,balls2,n2)
      sum(draw1*draw2)
    })
    onesim
  }

  ball_sim_sameCMP <- function(nn,n,meanV,nu){
    res = replicate(nn,{
      balls = cmp$simulate(c(meanV,nu))$sims
      draws = rmvhyper(1,balls,n)
      pairs = rowSums(choose(draws,2))
      pairs
    })
    res
  }

  ##62*125 is same as choose(125,2)

  exp_diff <- function(pairs,num_par){
    pairs*(1/num_par)
  }

  exp_same_given_CMP <- function(n,N,meanV,nu){
    densCMP = N*cmpNS$report(c(meanV,nu))$compPDF
    totballs = N*meanV
    exppairs =  lapply(0:4999,function(x){
      pairs = choose(0:x,2)
      sum(pairs*dhyper(0:x,x,max(0,totballs-x),n))
    })
    exppairsV = do.call(rbind,exppairs)
    sum(densCMP*exppairsV)
  }





  sim_nus = c(exp(seq(-5,5,by=0.1)))

  exp_CMPs = sapply(sim_nus,function(x){
    exp_same_given_CMP(125,5000,5,x)
  })

  same_case = list()
  diff_case = list()
  for(i in 1:length(sim_nus)){
    print(i)
    same_case[[i]] = mean(ball_sim_sameCMP(1000,125,5,sim_nus[i]))
    diff_case[[i]] = mean(ball_sim_diffCMP(1000,62,125,5,sim_nus[i]))
  }

  scc =  do.call(rbind,same_case)
  dcc = do.call(rbind,diff_case)

  df = data.frame(nus=log(sim_nus),same_case_theo=exp_CMPs,diff_case_theo=exp_diff(62*125,5000),same_case=scc,diff_case=dcc)
  saveRDS(df,file="./simplesim/sameVsDiff.rds")

  exp_diff(62*125,5000)

  dum_fix <- function(x){
    ifelse(x == 1,0,x)
  }

library(extraDistr)

  od_rpois <- function(n,lambda,theta){
    E <- rgamma(n,theta)/theta
    Y <- rpois(n,lambda*E)
  }

  od_ball_sim_same <- function(nn,n,N,lambda,theta){
    res = replicate(nn,{
      balls = od_rpois(N,lambda,theta)
      draws = rmvhyper(1,balls,n)
      pairs =colSums(apply(draws,1,dum_fix))
      pairs
    })
    res
  }


  od_ball_sim_samePairs <- function(nn,n,N,lambda,theta){
    res = replicate(nn,{
      balls = od_rpois(N,lambda,theta)
      draws = rmvhyper(1,balls,n)
      pairs = colSums(apply(draws,1,function(x){choose(x,2)}))
      pairs
    })
    res
  }


  thetas = c(0.1,0.25,0.5,0.75,1,1.25,2,5)
  set.seed(42)
  superDats = replicate(1000,{dats = lapply(thetas,function(x){
    n = sample(100:200,10,TRUE)
    dat = data.frame(n=n)
    dat$N = 5000
    dat$theta = x
    dat$lambda = 10
    dat$obs = apply(dat,1,function(x){
      od_ball_sim_same(1,x[1],x[2],x[4],x[3])})
    dat
  })})

 superDatsPairs = replicate(1000,{dats = lapply(thetas,function(x){
    n = sample(100:200,10,TRUE)
    dat = data.frame(n=n)
    dat$N = 5000
    dat$theta = x
    dat$lambda = 10
    dat$obs = apply(dat,1,function(x){
      od_ball_sim_samePairs(1,x[1],x[2],x[4],x[3])})
    dat
  })})



  compile("./simplesim/models/EXPOG.cpp")
  dyn.load(dynlib("./simplesim/models/EXPOG"))


  resultsTMB <- lapply(superDats,function(x){
    parmI = list()
    parmI$log_theta = log(1)
    parmI$log_N = log(5000)
    parmI$log_lambda = log(10)
    data = list()
    data$datan = as.matrix(x)
    mapp = list(log_N = as.factor(c(NA)),log_lambda= as.factor(c(NA)))
    obj = TMB::MakeADFun(data=data,parameters = parmI,map=mapp,DLL="EXPOG")
    opt = optim(obj$par,obj$fn,obj$gr,method="Brent",lower=-10,upper=10)
    ret = c(est=exp(opt$par),theta=unique(x$theta))
    ret
  })

resultsTMBMP <- lapply(superDatsPairs,function(x){
  parmMP = list()
  parmMP$log_V_A = log(1)
  parmMP$log_N = log(5000)
  parmMP$log_E_A = log(10)
  data = list()
  data$datan = as.matrix(x)
  mapp = list(log_N = as.factor(c(NA)),log_E_A= as.factor(c(NA)))
  obj = TMB::MakeADFun(data=data,parameters = parmMP,method="Brent",lower=-10,upper=10,map=mapp,DLL="mprobsim")
  opt = optim(obj$par,obj$fn,obj$gr)
  lam = unique(x$lambda)
  tvar = lam+lam^2/unique(x$theta)
  ret = c(V_A=exp(opt$par),tV_A=tvar)
})

resultsTMBMPFreeN <- lapply(superDatsPairs,function(x){
  parmMP = list()
  parmMP$log_V_A = log(1)
  parmMP$log_N = log(5000)
  parmMP$log_E_A = log(10)
  data = list()
  data$datan = as.matrix(x)
  mapp = list(log_N = as.factor(c(NA)),log_E_A= as.factor(c(NA)))
  obj = TMB::MakeADFun(data=data,parameters = parmMP,map=mapp,DLL="mprobsim")
  opt = optim(obj$par,obj$fn,obj$gr)
  rep = obj$report()
  lam = unique(x$lambda)
  tvar = lam+lam^2/unique(x$theta)
  ret = c(V_A=rep$V_A,N=rep$N,E_A=rep$E_A,tV_A=tvar)
})


dfsMP = as.data.frame(do.call(rbind,resultsTMBMP)) |>
  group_by(tV_A) |>
  summarise(quantile = scales::percent(c(0.05, 0.5, 0.95)),
	    thetahat = quantile(V_A.log_V_A, c(0.05, 0.5, 0.95)))


  dfs3MP =  spread(dfsMP,key="quantile",value="thetahat")


  dfs =  do.call(rbind,resultsTMB)
  library(tidyverse)
  dfs2 = as.data.frame(dfs) |>
    group_by(theta) |>
    summarise(quantile = scales::percent(c(0.05, 0.5, 0.95)),
	      thetahat = quantile(est, c(0.05, 0.5, 0.95)))
  saveRDS(dfs,file="./simplesim/thetaTest.rds")

  dfs3 =  spread(dfs2,key="quantile",value="thetahat")

  resultsTMBFreeN <- lapply(superDats,function(x){
    parmI = list()
    parmI$log_theta = log(1)
    parmI$log_N = log(2000)
    parmI$log_lambda = log(10)
    ##parmI$log_theta_sum = log(1)
    data = list()
    data$datan = as.matrix(x)
    mapp = list(log_N = as.factor(c(1)),log_lambda= as.factor(c(NA)))
    obj = TMB::MakeADFun(data=data,parameters = parmI,map=mapp,DLL="EXPOG")
    opt = optim(obj$par,obj$fn,obj$gr)
    ret = c(est=exp(opt$par),theta=unique(x$theta))
    ret
  })

  dfsFreeN =  do.call(rbind,resultsTMBFreeN)
  saveRDS(dfsFreeN,file="./simplesim/thetaNtest.rds")



  dfsFreeN2 = as.data.frame(dfsFreeN) |>
    group_by(theta) |>
    summarise(quantile = scales::percent(c(0.05, 0.5, 0.95)),
	      thetahat = quantile(est.log_theta, c(0.05, 0.5, 0.95)),
	      Nhat = quantile(est.log_N,c(0.05,0.5,0.95)))

library(ggplot2)

dfsFN3 = as.data.frame(dfsFreeN)

dfsFN3sp = split(dfsFN3,dfsFN3$theta)

ggplot(dfsFN3sp[[3]],aes(x=log(dfsFN3sp[[3]]$est.log_theta))) + geom_density()

densTheta3 = density(log(dfsFN3sp[[3]]$est.log_theta))


densMAX = function(dat){
  dens = density(dat)
  wmax = which.max(dens$y)
  xmax =dens$x[wmax]
  ret = c(xmax,wmax)
  ret
}

maxThetas = lapply(dfsFN3sp,function(x){
  densMAX(log(x$est.log_theta))
  })
