library(readr)

  folders = list.dirs("../sims2",full.names = TRUE)
  ffolders = folders[grep("/fitted",folders)]

  rmse = function(x,xtrue){
    sqrt(sum((x-xtrue)^2)/length(x))
  }

  popbit =  regexec("pop[0-9][0-9][0-9][0-9]",ffolders)
  mlength = sapply(popbit,function(x){attributes(x)$match.length})
  popID = unique(do.call(substring,list(text=ffolders,first=popbit,last=unlist(popbit)+mlength-1)))

  get_dat_from_fitted <- function(file,sim){
    dat = readRDS(file)
    dat$Mest$sim = sim
    dat$Nest$sim = sim
    dat$Nest$year = sort(rep(1:11,3))
    dat$FecEst$sim = sim
    dat$FecEst$age = 1:3
    dat$thetaEst$sim = sim
    dat$VkEst$sim = sim
    dat$growth_rateEst$sim = sim

    dat$MestSC$sim = sim
    dat$NestSC$sim = sim
    dat$NestSC$year = sort(rep(1:11,3))
    dat$FecEstSC$sim = sim
    dat$FecEstSC$age = 1:3
    dat$growth_rateEstSC$sim = sim
    dat$N_eEst$sim = sim
    dat$N_eEst$year = 1:11

    ##RMSE
    dat$Mest$RMSE = apply(dat$Mest,1,function(x){
      rmse(x[7],x[8])
      })
    dat$Nest$RMSE = apply(dat$Nest,1,function(x){
      rmse(x[7],x[8])
      })
    dat$FecEst$RMSE = apply(dat$FecEst,1,function(x){
      rmse(x[7],x[8])
    })
    dat$thetaEst$RMSE = apply(dat$thetaEst,1,function(x){
      rmse(x[7],x[8])
    })
    dat$VkEst$RMSE = apply(dat$VkEst,1,function(x){
      rmse(x[1],x[5])
    })
    dat$growth_rateEst$RMSE = apply(dat$growth_rateEst,1,function(x){
      rmse(x[1],x[5])
    })
    Nrmse = rmse(dat$Nest$Nest,dat$Nest$Ntrue)

    dat$MestSC$RMSE = apply(dat$MestSC,1,function(x){
      rmse(x[7],x[8])
    })
    dat$NestSC$RMSE = apply(dat$NestSC,1,function(x){
      rmse(x[7],x[8])
    })
    dat$FecEstSC$RMSE = apply(dat$FecEstSC,1,function(x){
      rmse(x[7],x[8])
    })
    dat$growth_rateEstSC$RMSE = apply(dat$growth_rateEstSC,1,function(x){
      rmse(x[1],x[5])
    })
    NrmseSC = rmse(dat$NestSC$NestSC,dat$NestSC$Ntrue)
    ret = with(dat,list(Mest=Mest,Nest=Nest,FecEst=FecEst,
			thetaEst=thetaEst,VkEst=VkEst,
			growth_rateEst=growth_rateEst,
			MestSC=MestSC,
			NestSC=NestSC,
			FecEstSC=FecEstSC,
			growth_rateEstSC=growth_rateEstSC,
			Nrmse=Nrmse,NrmseSC=NrmseSC,
			N_eEst = N_eEst,
			converged=dat$converged,convergedSC=dat$convergedSC))
    ret
  }



  get_dat_from_fittedNEW <- function(file,sim){
    dat = readRDS(file)
    dat$Mest$sim = sim
    dat$Nest$sim = sim
    dat$Nest$year = sort(rep(1:11,3))
    dat$FecEst$sim = sim
    dat$FecEst$age = 1:3
    dat$thetaEst$sim = sim
    dat$VkEst$sim = sim
    dat$growth_rateEst$sim = sim
    dat$N_eEst$sim = sim
    dat$N_eEst$year = 1:11



    dat$MestSC$sim = sim
    dat$NestSC$sim = sim
    dat$NestSC$year = sort(rep(1:11,3))
    dat$FecEstSC$sim = sim
    dat$FecEstSC$age = 1:3
    dat$growth_rateEstSC$sim = sim

    dat$MestMB$sim = sim
    dat$NestMB$sim = sim
    dat$NestMB$year = sort(rep(1:11,3))
    dat$FecEstMB$sim = sim
    dat$FecEstMB$age = 1:3
    dat$thetaEstMB$sim = sim
    dat$VkEstMB$sim = sim
    dat$growth_rateEstMB$sim = sim
    dat$N_eEstMB$sim = sim
    dat$N_eEstMB$year = 1:11



    ##RMSE
    dat$Mest$RMSE = apply(dat$Mest,1,function(x){
      rmse(x[7],x[8])
      })
    dat$Nest$RMSE = apply(dat$Nest,1,function(x){
      rmse(x[7],x[8])
      })
    dat$FecEst$RMSE = apply(dat$FecEst,1,function(x){
      rmse(x[7],x[8])
    })
    dat$thetaEst$RMSE = apply(dat$thetaEst,1,function(x){
      rmse(x[7],x[8])
    })
    dat$VkEst$RMSE = apply(dat$VkEst,1,function(x){
      rmse(x[1],x[5])
    })
    dat$growth_rateEst$RMSE = apply(dat$growth_rateEst,1,function(x){
      rmse(x[1],x[5])
    })
    Nrmse = rmse(dat$Nest$Nest,dat$Nest$Ntrue)

    dat$MestSC$RMSE = apply(dat$MestSC,1,function(x){
      rmse(x[7],x[8])
    })
    dat$NestSC$RMSE = apply(dat$NestSC,1,function(x){
      rmse(x[7],x[8])
    })
    dat$FecEstSC$RMSE = apply(dat$FecEstSC,1,function(x){
      rmse(x[7],x[8])
    })
    dat$growth_rateEstSC$RMSE = apply(dat$growth_rateEstSC,1,function(x){
      rmse(x[1],x[5])
    })
    NrmseSC = rmse(dat$NestSC$NestSC,dat$NestSC$Ntrue)

    dat$MestMB$RMSE = apply(dat$MestMB,1,function(x){
      rmse(x[7],x[8])
      })
    dat$NestMB$RMSE = apply(dat$NestMB,1,function(x){
      rmse(x[7],x[8])
      })
    dat$FecEstMB$RMSE = apply(dat$FecEstMB,1,function(x){
      rmse(x[7],x[8])
    })
    dat$thetaEstMB$RMSE = apply(dat$thetaEstMB,1,function(x){
      rmse(x[7],x[8])
    })
    dat$VkEstMB$RMSE = apply(dat$VkEstMB,1,function(x){
      rmse(x[1],x[5])
    })
    dat$growth_rateEstMB$RMSE = apply(dat$growth_rateEstMB,1,function(x){
      rmse(x[1],x[5])
    })
    NrmseMB = rmse(dat$NestMB$Nest,dat$NestMB$Ntrue)



    ret = with(dat,list(Mest=Mest,Nest=Nest,FecEst=FecEst,
			thetaEst=thetaEst,VkEst=VkEst,
			growth_rateEst=growth_rateEst,
			MestSC=MestSC,
			NestSC=NestSC,
			FecEstSC=FecEstSC,
			growth_rateEstSC=growth_rateEstSC,
			Nrmse=Nrmse,NrmseSC=NrmseSC,
			N_eEst = N_eEst,
			converged=dat$converged,convergedSC=dat$convergedSC,
			MestMB=MestMB,NestMB=NestMB,
			FecEstMB=FecEstMB,thetaEstMB=thetaEstMB,
			VkEstMB=VkEstMB,growth_rateEstMB=growth_rateEstMB,
			NrmseMB=NrmseMB))

    ret
  }

  get_data_from_files <- function(folder){
    files = list.files(folder,"*.rds")
    filesFN = list.files(folder,"*.rds",full.names = TRUE)
    sims = parse_number(files)

    ##Do the first one
    dat =  get_dat_from_fitted(filesFN[1],sims[1])

    for(i in 2:length(filesFN)){
      datTem = get_dat_from_fitted(filesFN[i],sims[i])
      dat = mapply(function(d,t){
	rbind(d,t)},d=dat,t=datTem)
    }
    dat
  }

  get_data_from_filesNEW <- function(folder){
    files = list.files(folder,"*.rds")
    filesFN = list.files(folder,"*.rds",full.names = TRUE)
    sims = parse_number(files)

    ##Do the first one
    dat =  get_dat_from_fittedNEW(filesFN[1],sims[1])

    for(i in 2:length(filesFN)){
      datTem = get_dat_from_fittedNEW(filesFN[i],sims[i])
      dat = mapply(function(d,t){
	rbind(d,t)},d=dat,t=datTem)
    }
    dat
  }


  library(ggplot2)

  quantileNE = readRDS("~/varCKMR3/quantileNE.rds")

  make_plots_from_data <- function(data,scenario){
      grdensitySC = ggplot(data$growth_rateEstSC,aes(x=estSCgrowth_rate))+
      geom_density() + geom_vline(aes(xintercept=unique(data$growth_rateEst$trueGR)),color="blue",linetype="dashed") + xlab("Growth Rate") + ggtitle(paste("Growth Rate - No same cohort  - ",scenario))
    MestDensitySC = ggplot(data$MestSC,aes(x=MestSC))+ geom_density() + geom_vline(aes(xintercept=data$Mest$Mtrue),color="blue",linetype="dashed") + xlab("M") + ggtitle(paste("M - No same cohort - ",scenario)) + facet_wrap(~age)
    FecEstDensitySC = ggplot(data$FecEstSC,aes(x=FecEstSC))+ geom_density() + geom_vline(aes(xintercept=data$FecEst$FecTrue),color="blue",linetype="dashed") + xlab("Fecundity") + ggtitle(paste("Fecundity - No Same Cohort - ",scenario)) + facet_wrap(~age,scales="free")
    NestDensitySC = ggplot(data$NestSC[data$NestSC$year == 1,],aes(x=NestSC)) + geom_density() + geom_vline(aes(xintercept=data$Nest$Ntrue[data$Nest$year == 1]),
												   color="blue",linetype="dashed") + xlab("Numbers at age")+
      ggtitle(paste("N - No same Cohort - ",scenario)) + facet_wrap(~age,scales="free")

      get_xy_lims <- function(plot){
	orb = ggplot_build(plot)$layout
	xranges = lapply(orb$panel_params,function(x){x$x.range})
	yranges = lapply(orb$panel_params,function(x){x$y.range})
	ret = list(xranges=xranges,yranges=yranges)
	ret
	}





    grdensity = ggplot(data$growth_rateEst,aes(x=estgrowth_rate))+
      geom_density() + geom_vline(aes(xintercept=unique(data$growth_rateEst$trueGR)),color="blue",linetype="dashed") + xlab("Growth Rate") + ggtitle(paste("Growth Rate - ",scenario))
      thetadensity = ggplot(data$thetaEst[data$thetaEst$thetaEst < 20,],aes(x=thetaEst))+
	geom_density() + geom_vline(aes(xintercept=unique(data$thetaEst$thetaTrue)),color="blue",linetype="dashed") + xlab("Theta") + ggtitle(paste("Theta - ",scenario))
    MestDensity = ggplot(data$Mest,aes(x=Mest))+ geom_density() + geom_vline(aes(xintercept=data$Mest$Mtrue),color="blue",linetype="dashed") + xlab("M") + ggtitle(paste("M - ",scenario)) + facet_wrap(~age)
    FecEstDensity = ggplot(data$FecEst,aes(x=FecEst))+ geom_density() + geom_vline(aes(xintercept=data$FecEst$FecTrue),color="blue",linetype="dashed") + xlab("Fecundity") + ggtitle(paste("Fecundity - ",scenario)) + facet_wrap(~age,scales="free")
    VkEstDensity = ggplot(data$VkEst,aes(x=estVk))+ geom_density() + geom_vline(aes(xintercept=data$VkEst$fromSim),color="blue",linetype="dashed") + xlab("Variance of Lifetime Reproductive Success") + ggtitle(paste("Vk - ",scenario))
    NestDensity = ggplot(data$Nest[data$Nest$year == 1,],aes(x=Nest)) + geom_density() + geom_vline(aes(xintercept=data$Nest$Ntrue[data$Nest$year == 1]),
												   color="blue",linetype="dashed") + xlab("Numbers at age")+
      ggtitle(paste("N - ",scenario)) + facet_wrap(~age,scales="free")


    ret = list(grdensity=grdensity,thetadensity=thetadensity,MestDensity=MestDensity,
	       FecEstDensity=FecEstDensity,VkEstDensity=VkEstDensity,NestDensity=NestDensity,
	       grdensitySC=grdensitySC,MestDensitySC=MestDensitySC,FecEstDensitySC=FecEstDensitySC,
	       NestDensitySC=NestDensitySC)
    ret
  }

  check_rmse_diffs <- function(data){
    grrmse = sum(data$growth_rateEst$RMSE < data$growth_rateEstMB$RMSE)/length(data$growth_rateEst$RMSE) 
    Mrmse = sum(data$Mest$RMSE < data$MestMB$RMSE)/length(data$Mest$RMSE)
    FecRmse = sum(data$FecEst$RMSE < data$FecEstMB$RMSE)/length(data$FecEst$RMSE)
    Nrmse = sum(data$Nest$RMSE < data$NestMB$RMSE)/length(data$Nest$RMSE)
    Noverall = sum(data$Nrmse < data$NrmseMB)/length(data$Nrmse)
    theta = sum(data$thetaEst$RMSE < data$thetaEstMB$RMSE)

    ret = data.frame(growth=grrmse,M=Mrmse,Fec=FecRmse,N=Nrmse,Noverall=Noverall)
    ret
  }


  quantile_bits <- function(data){
    probs = c(0.05,0.5,0.95)
    fecA = split(data$FecEst,data$FecEst$age)
    fecAQ = t(sapply(fecA,function(x){quantile(x$FecEst,probs)}))
    fecAQ = as.data.frame(cbind(fecAQ,true=unique(data$FecEst$FecTrue)))

    fecASC = split(data$FecEstSC,data$FecEstSC$age)
    fecAQSC = as.data.frame(t(sapply(fecASC,function(x){quantile(x$FecEst,probs)})))

    fecC = data.frame(parm=c("fec1","fec2","fec3"))
    fecC = cbind(fecC,fecAQ,fecAQSC)

    MA = split(data$Mest,data$Mest$age)
    MAQ = t(sapply(MA,function(x){quantile(x$Mest,probs)}))
    MAQ = as.data.frame(cbind(MAQ,true=unique(data$Mest$Mtrue)))

    MASC = split(data$MestSC,data$MestSC$age)
    MAQSC = as.data.frame(t(sapply(MASC,function(x){quantile(x$MestSC,probs)})))

    MC = data.frame(parm=c("M1","M2"))
    MC = cbind(MC,MAQ,MAQSC)

    NestA = split(data$Nest[data$Nest$year == 5,],data$Nest$age[data$Nest$year == 5])
    NestAQ = t(sapply(NestA,function(x){quantile(x$Nest,probs)}))
    NestAQ = as.data.frame(cbind(NestAQ,true=unique(data$Nest$Ntrue[data$Nest$year == 5])))

  NestASC = split(data$NestSC[data$NestSC$year == 5,],data$NestSC$age[data$NestSC$year == 5])
    NestAQSC = as.data.frame(t(sapply(NestASC,function(x){quantile(x$NestSC,probs)})))

    NestC = data.frame(parm=c("N1Y5","N2Y5","N3Y5"))
    NestC = cbind(NestC,NestAQ,NestAQSC)


    thetaA = c("theta",quantile(data$thetaEst$thetaEst,probs),unique(data$thetaEst$thetaTrue),NA,NA,NA)

    combo = rbind(fecC,MC,NestC,thetaA)
    combo[,2:8] = apply(combo[,2:8],2,as.numeric)
    combo[,2:8] = round(combo[,2:8],3)
    combo
    }



  quantile_bits2 <- function(data){
    probs = c(0.05,0.5,0.95)
    fecA = split(data$FecEst,data$FecEst$age)
    fecAQ = t(sapply(fecA,function(x){quantile(x$FecEst,probs)}))
    fecAQ = as.data.frame(cbind(fecAQ,true=unique(data$FecEst$FecTrue)))

    fecAMB = split(data$FecEstMB,data$FecEstMB$age)
    fecAQMB = as.data.frame(t(sapply(fecAMB,function(x){quantile(x$FecEst,probs)})))

    fecC = data.frame(parm=c("fec1","fec2","fec3"))
    fecC = cbind(fecC,fecAQ,fecAQMB)

    MA = split(data$Mest,data$Mest$age)
    MAQ = t(sapply(MA,function(x){quantile(x$Mest,probs)}))
    MAQ = as.data.frame(cbind(MAQ,true=unique(data$Mest$Mtrue)))

    MAMB = split(data$MestMB,data$MestMB$age)
    MAQMB = as.data.frame(t(sapply(MAMB,function(x){quantile(x$MestMB,probs)})))

    MC = data.frame(parm=c("M1","M2"))
    MC = cbind(MC,MAQ,MAQMB)

    NestA = split(data$Nest[data$Nest$year == 5,],data$Nest$age[data$Nest$year == 5])
    NestAQ = t(sapply(NestA,function(x){quantile(x$Nest,probs)}))
    NestAQ = as.data.frame(cbind(NestAQ,true=unique(data$Nest$Ntrue[data$Nest$year == 5])))

  NestAMB = split(data$NestMB[data$NestMB$year == 5,],data$NestMB$age[data$NestMB$year == 5])
    NestAQMB = as.data.frame(t(sapply(NestAMB,function(x){quantile(x$NestMB,probs)})))

    NestC = data.frame(parm=c("N1Y5","N2Y5","N3Y5"))
    NestC = cbind(NestC,NestAQ,NestAQMB)

    Vk = quantile(data$VkEst$estVk,probs) 
    VkMB = quantile(data$VkEstMB$estVk,probs)
    VkA = c("Vk",Vk,unique(data$VkEst$fromSimPar),VkMB)
  
    thetaA = c("theta",quantile(data$thetaEst$thetaEst,probs),unique(data$thetaEst$thetaTrue),quantile(data$thetaEstMB$thetaEstMB,probs))

    combo = rbind(fecC,MC,NestC,thetaA,VkA)
    combo[,2:8] = apply(combo[,2:8],2,as.numeric)
    combo[,2:8] = round(combo[,2:8],3)
    combo
    }



  quantile_ldne <- function(data){
    ldne = data$N_eEst
    ldney = split(ldne,ldne$year)
    qldne = sapply(ldney,function(x){
      quantile(x$LDNE,c(0.05,0.5,0.95))
    })

    qldneE = sapply(ldney,function(x){
      quantile(x$estN_e,c(0.05,0.5,0.95))
    })
    tqldne = as.data.frame(t(qldne))
    tqldne$type = "AdjNe"
    tqldneE = as.data.frame(t(qldneE))
    tqldneE$type = "ModelEstimate"
    cbind(tqldne,tqldneE)
  }

  quantile_neMB <- function(data){
    ldne = data$N_eEst
    ldney = split(ldne,ldne$year)
    ## qldne = sapply(ldney,function(x){
    ##   quantile(x$LDNE,c(0.05,0.5,0.95))
    ## })

    qldneE = sapply(ldney,function(x){
      quantile(x$estN_e,c(0.05,0.5,0.95))
    })
  #  tqldne = as.data.frame(t(qldne))
   # tqldne$type = "AdjNe"
    tqldneE = as.data.frame(t(qldneE))
    tqldneE$type = "MBModelEstimate"
    cbind(tqldne,tqldneE)
  }
M




  plots = list()
  rmses = list()
  converged = list()
  convergedSC = list()
  quantiles = list()
  quantileNE = list()
  popPars = list()
  for(i in 1:length(ffolders)){
    data = get_data_from_filesNEW(ffolders[i])
    quantiles[[i]] = quantile_bits2(data)
   # quantileNE[[i]] = quantile_ldne(data)
  }

  #saveRDS(quantileNE,file="~/varCKMR3/quantileNE.rds")


  for(i in 1:length(ffolders)){
    indiv = readRDS(paste0("../populations/",popID[i],".rds"))
    popPars[[i]] = indiv$popPars
    sim_scenario = paste("Growth Rate: ",indiv$popPars$growth_rate,"Theta: ",indiv$popPars$theta)
    data = get_data_from_filesNEW(ffolders[i])
    converged[[i]] = data$converged
    convergedSC[[i]] = data$convergedSC
    plots[[i]] = make_plots_from_data(data,sim_scenario)
    rmses[[i]] = check_rmse_diffs(data)
    rmses[[i]]$scenario = sim_scenario
  }


  pPars = do.call(rbind,popPars)
  library(tidyverse)

  gen_length2 <- function(fec,surv,lambda){
    ells = cumprod(c(1,surv))
    nn = length(fec)
    T = 0.0
    for(i in 1:nn){
      T = T + i*lambda^(-i)*ells[i]*fec[i]
    }
    T
  }

  dprobs = function(surv){
    dprob = numeric(length(surv)+1)
    for(i in 1:length(dprob)){
      if(i == 1){
	dprob[i] = 1-surv[i]
      }else if(i == length(surv)+1){
	dprob[i] = tail(cumprod(surv[1:(i-1)]),1)
      }else{
	dprob[i] = tail(cumprod(surv[1:(i-1)]),1)*(1-surv[i])
      }

    }
    dprob
  }


  crowNe <- function(N,k,L,Vk){
    top = (2*N-1)*k*L
    bot = 2*(1+Vk/k)
    top/bot
  }


  Neplots = list()
  Neplots2 = list()
  for(i in 1:length(ffolders)){
    data = get_data_from_files(ffolders[i])
    NTru = data$Nest |>
      group_by(year) |>
      distinct(Ntrue) |>
      summarise(Ntrue=sum(Ntrue))

    NTru1 = data$Nest |>
      filter(age == 1) |>
      group_by(year) |>
      distinct(Ntrue)

    Vk = unique(data$VkEst$fromSimPar)
    fff = unique(data$FecEst$FecTrue)
    sss = exp(-unique(data$Mest$Mtrue))
    dp = dprobs(sss)
    kbar = sum(cumsum(fff*2)*dp)
  

    sim = 500
    zorp = data$N_eEst[data$N_eEst$sim == sim,]
    gr = pPars[i,2]
    theta = pPars[i,3]
    gl = gen_length2(fff,sss,gr)

  crowNes = numeric(nrow(NTru))
    for(j in 1:length(crowNes)){
      crowNes[j] = crowNe(NTru1[j,1],kbar,gl,Vk)
    }

    crowNes = as.vector(do.call(rbind,crowNes))
  
    qNe = quantileNE[[i]]
    N = c(crowNes,qNe[,2],qNe[,6])
    Ntype = c(rep("FromSimPar",11),qNe[,4],qNe[,8])
    NLo = c(rep(NA,11),qNe[,1],qNe[,5])
    NHi = c(rep(NA,11),qNe[,3],qNe[,7])
    df = data.frame(year=rep(40:50,3),N=N,Ntype=Ntype,NLo=NLo,NHi=NHi)
    Neplots[[i]] = ggplot(df) + geom_line(aes(x=year,y=N,group=Ntype,color=Ntype)) + geom_ribbon(aes(x=year,ymin=NLo,ymax=NHi,group=Ntype,fill=Ntype),alpha=0.2) +ggtitle(paste("Vk=",round(Vk,2),"Growth Rate=",round(gr,2),"Theta=",round(theta,2)))

    df2 = data.frame(year=rep(40:50,4),N=c(zorp[,1],zorp[,6],NTru$Ntrue,qNe[,2]),Ntype=c(rep("estNe",11),rep("LDNE1",11),rep("N",11),rep("AdjNe",11)),
		     NLo = c(zorp[,3],rep(NA,11),rep(NA,11),qNe[,1]),NHi = c(zorp[,4],rep(NA,11),rep(NA,11),qNe[,3]))
      Neplots2[[i]] = ggplot(df2) + geom_line(aes(x=year,y=N,group=Ntype,color=Ntype)) + geom_ribbon(aes(x=year,ymin=NLo,ymax=NHi,group=Ntype,fill=Ntype),alpha=0.2) +ggtitle(paste("Vk=",round(Vk,2),"Growth Rate=",round(gr,2),"Theta=",round(theta,2)))

  }



  pdf(file="Neplots.pdf")
  for(i in 1:length(Neplots)){
    print(Neplots[[i]])
  }
  dev.off()

  coup = do.call(rbind,rmses)
  sim_scenario = coup$scenario

  conv = data.frame(scenario=sim_scenario,conv=sapply(converged,sum))


  RMSEtab = coup[,c(6,5,1,2,3)]
  RMSEtab[,2:5] = round(RMSEtab[,2:5],2)




  ##Make 
  data = get_data_from_files(ffolders[4])

  scenario = paste("Growth Rate: ",0.95," Theta: ", 0.75)
    NestDensitySC = ggplot(data$NestSC[data$NestSC$year == 1,],aes(x=NestSC)) + geom_density() + geom_vline(aes(xintercept=data$Nest$Ntrue[data$Nest$year == 1]),
												   color="blue",linetype="dashed") + xlab("Numbers at age")+
      ggtitle(paste("N - within-cohort excluded- ",scenario)) + facet_wrap(~age,scales="free")

  library(cowplot)
