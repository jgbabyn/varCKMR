library(TMB)
  library(tidyverse)
  library(stringr)
  library(varCKMR)

  if(!file.exists("../models/ssmodelAdj.o")){
    compile("../models/ssmodelAdj.cpp")
  }
  if(!file.exists("../models/ssmodelNoSC.o")){
    compile("../models/ssmodelNoSC.cpp")
  }
  if(!file.exists("../models/ssmodelMVP.o")){
    compile("../models/ssmodelMVP.cpp")
  }



  dyn.load(dynlib("../models/ssmodelAdj"))
  dyn.load(dynlib("../models/ssmodelNoSC"))
  dyn.load(dynlib("../models/ssmodelMVP"))


  rep_df = readRDS(repdf)


      for(k in 1:1){

	files = list.files(paste0(rep_df[[k]]$foldname[1]),pattern = "*.rds",full.names = TRUE)
	popbit =  regexec("pop[0-9][0-9][0-9][0-9]",files)
	mlength = sapply(popbit,function(x){attributes(x)$match.length})
	popID = unique(do.call(substring,list(text=files,first=popbit,last=unlist(popbit)+mlength-1)))

	for(rr in 1:length(files)){
	  print(rr)
	  stuff = readRDS(files[rr])

	  dat = stuff$data
	  dat$skipsame = 0
	  parm = stuff$parm
	  parm$log_init_tot_N = 10

	mapp2 = list(log_M=as.factor(c(1,2,NA)),log_init_tot_N=as.factor(c(1)),
		     log_fecundity=as.factor(c(1,2,3)),log_theta=as.factor(c(1)))


	##No SC cohort
	datSC = stuff$data
	datSC$SIBsSC = NULL
	parmSC = parm
	parmSC$log_theta = NULL
	mapp2SC = mapp2
	mapp2SC$log_theta = NULL
	objSC = MakeADFun(datSC,parmSC,map = mapp2SC, DLL="ssmodelNoSC")
	optSC = nlminb(objSC$par,objSC$fn,objSC$gr,control=list(iter.max=5000))
	reportSC = objSC$report()
	sdrSC = sdreport(objSC)
	  ssdrSC = summary(sdrSC)

	  ##MB version
	  datMB = stuff$data
	  datMB$skipsame = 0
	  parmMB = stuff$parm

	  objMB = MakeADFun(datMB,parmMB,map=mapp2,DLL="ssmodelAdj")
	  optMB = nlminb(objMB$par,objMB$fn,objMB$gr,control=list(iter.max=5000))
	  reportMB = objMB$report()
	  sdrMB = sdreport(objMB)
	  ssdrMB = summary(sdrMB)


	  ##MVP sensitivity analysis
	  ##Theta here is "c", var_power is gamma
	  datMVP2 = stuff$data
	  datMVP2$skipsame = 0
	  parmMVP2 = stuff$parm
	  parmMVP2$var_power = 2
	  mapp3 = mapp2
	  mapp3$log_theta = as.factor(1)
	  mapp3$var_power = as.factor(NA)

	  objMVP2 = MakeADFun(datMVP2,parmMVP2,map=mapp3,DLL="ssmodelMVP")
	  optMVP2 = nlminb(objMVP2$par,objMVP2$fn,objMVP2$gr,control=list(iter.max=5000))
	  reportMVP2 = objMVP2$report()
	  sdrMVP2 = sdreport(objMVP2)
	  ssdrMVP2 = summary(sdrMVP2)


	  datMVP1 = stuff$data
	  datMVP1$skipsame = 0
	  parmMVP1 = stuff$parm
	  parmMVP1$var_power = 2
	  mapp3 = mapp2
	  mapp3$log_theta = as.factor(1)
	  mapp3$var_power = as.factor(NA)

	  objMVP1 = MakeADFun(datMVP1,parmMVP1,map=mapp3,DLL="ssmodelMVP")
	  optMVP1 = nlminb(objMVP1$par,objMVP1$fn,objMVP1$gr,control=list(iter.max=5000))
	  reportMVP1 = objMVP1$report()
	  sdrMVP1 = sdreport(objMVP1)
	  ssdrMVP1 = summary(sdrMVP1)



	##Find variance 
	indiv = readRDS(paste0("../populations/",popID,".rds"))
	ped = make_ped(indiv)


	  pp = ped$pedigree

	  momkids = pp |>
	    group_by(mother) |>
	    summarise(n=n()) |>
	    dplyr::rename(ID=mother)

	  dadkids = pp |>
	    group_by(father) |>
	    summarise(n=n()) |>
	    dplyr::rename(ID=father)

	  zorp = data.frame(ID=pp$ID,birth_year=pp$birth_year)
	  zorp = left_join(zorp,momkids)
	  zorp = left_join(zorp,dadkids,by="ID")
	  zorp[is.na(zorp)] = 0
	  zorp$kids = zorp$n.x+zorp$n.y

	  zack = zorp |>
	    group_by(birth_year) |>
	    summarise(var=var(kids),mean=mean(kids))

	    estvar = mean(zack[35:46,2,TRUE])
	    estmean = mean(zack[35:46,3,TRUE])

	  tN = n_at_age(indiv)[,]

	  ##This just takes a very long time...
	  ##N_e sampled
	  ## library(abind)
	  ## library(adegenet)
	  ## library(strataG)

	  ## sampsForLDNe = varCKMR::rough_sample_size(t(tN),function(y){0.15*y})
	  ## sampForLDNe = varCKMR::sample_age_by_year(indiv$population,sampsForLDNe,40:50,1)
	  ## samp2ForLDNe = lapply(sampForLDNe,function(x){
	  ##   do.call(rbind,x)
	  ## })

	  ## getLDNE = function(ped){
	  ##   pp = lapply(ped,function(by50){
	  ##     Agroup = select(by50,ends_with("_A"))[,1:100]
	  ##     Bgroup = select(by50,ends_with("_B"))[,1:100]

	  ##     ABpaste = mapply(function(x,y){
	  ##       paste0(x,y)},Agroup,Bgroup)

	  ##     genindy = df2genind(ABpaste,ploidy=2,ncode=1)
	  ##     gi.g <- genind2gtypes(genindy)
	  ##     ldNeey = ldNe(gi.g,maf.threshold = 0.05,num.cores = 8)
	  ##     ldNeey
	  ##   })
	  ##   do.call(rbind,pp)
	  ## }

	  ## getLDNE2 = function(N){
	  ##   ldNe = list()
	  ##   for(i in 1:N){
	  ##     sampsForLDNe = varCKMR::rough_sample_size(t(tN),function(y){0.15*y})
	  ##     sampForLDNe = varCKMR::sample_age_by_year(indiv$population,sampsForLDNe,40:50,1)
	  ##     samp2ForLDNe = lapply(sampForLDNe,function(x){
	  ##       do.call(rbind,x)
	  ##     })
	  ##     ldNeT = getLDNE(samp2ForLDNe)
	  ##     ldNeT$rawNb = 1/(3*(ldNeT[,4]-ldNeT[,5]))
	  ##     CVf = sd(fec)/mean(fec)
	  ##     AL = 3-1+1
	  ##     ldNeT$Nbadj = ldNeT$rawNb/(0.991-0.206*log10(AL)+0.256*log10(3)+0.137*CVf)
	  ##     ldNeT$Neadj = ldNeT$Nbadj/(0.833+0.637*log10(AL)-0.793*log10(3)-0.423*CVf)
	  ##     ldNe[[i]] = ldNeT
	  ##   }
	  ##   ldNe
	  ## }
	  ## dorp = getLDNE2(10)
	  ## zap = lapply(1:10,function(x){
	  ##   y = dorp[[x]][,"Neadj"]
	  ##   y
	  ##   })
	  ## zap2 = do.call(rbind,zap)
	  ## qzap2 = apply(zap2,2,function(x){
	  ##   quantile(x,c(0.05,0.5,0.95))})


	  ## dfdum = data.frame(ModelN=colSums(exp(report$log_N)),TrueN=rowSums(tN[,]),ModelNe=report$N_e,CrowNe=cNe,LDNE50=qzap2[2,],LDNE5=qzap2[1,],LDNE95=qzap2[3,],year=40:50)
	  ## dfdum2 = dfdum |>
	  ##   gather(key=Ntype,value=N,-LDNE5,-LDNE95,-year)

	  ## library(ggplot2)
	  ## dp1 = ggplot(dfdum2) + geom_line(aes(x=year,y=N,group=Ntype,color=Ntype)) + geom_ribbon(aes(x=year,ymin=LDNE5,ymax=LDNE95),alpha=0.3)

	  popPar = indiv$popPars  

	  dprob = c(1-popPar$surv1,popPar$surv1*(1-popPar$surv2),popPar$surv2*popPar$surv1)
	  fec = c(popPar$fec1,popPar$fec2,popPar$fec3)



	  ## ldNeEst = getLDNE(samp2ForLDNe)
	  ## ldNeEst$rawNb = 1/(3*(ldNeEst[,4]-ldNeEst[,5]))
	  ## CVf = sd(fec)/mean(fec)
	  ## AL = 3-1+1
	  ## ldNeEst$Nbadj = ldNeEst$rawNb/(0.991-0.206*log10(AL)+0.256*log10(3)+0.137*CVf)
	  ## ldNeEst$Neadj = ldNeEst$Nbadj/(0.833+0.637*log10(AL)-0.793*log10(3)-0.423*CVf)



	condvar <- function(dprob,fec,theta){
	  firstsum = 0.0
	  varsy = fec+fec^2/theta
	  for(i in 1:3){
	    firstsum = firstsum + dprob[i]*sum(varsy[1:i])
	  }
	  secsum = 0.0
	  for(i in 1:3){
	    secsum = secsum + (sum(fec[1:i]))^2*(1-dprob[i])*dprob[i]
	  }
	  lastsum = 0.0
	  for(i in 2:3){
	    for(j in 1:(i-1)){
	      lastsum = lastsum + sum(fec[1:i])*dprob[i]*sum(fec[1:j])*dprob[j]
	    }
	  }
	  tot = firstsum+secsum-2*lastsum
	  tot
	}

	popPar = indiv$popPars  

	dprob = c(1-popPar$surv1,popPar$surv1*(1-popPar$surv2),popPar$surv2*popPar$surv1)
	fec = c(popPar$fec1,popPar$fec2,popPar$fec3)
	surv = c(popPar$surv1,popPar$surv2)
	les = make_les(fec,surv)

	  doub_fec = 2*fec
	  thet = popPar$theta
	  varsy = doub_fec+doub_fec^2/thet


	  growth_rate = max(Re(eigen(les,only.values = TRUE)$values))
	  surv2 = cumprod(c(1,surv,0))
	  surv3 = c(surv,0)
	  fec2 = c(fec,0)

	  repro_value <- function(a){
	    lx = surv2[a]
	    vx = 0.0
	    for(i in a:4){
	      vx = vx+growth_rate^(-i-a+1)*surv2[i]/lx*fec2[i]
	    }
	    vx
	  }

	  rvs = c(repro_value(1),repro_value(2),repro_value(3))
	  denom = 0.0
	  for(i in 1:3){
	    out = surv2[i]*growth_rate^(-i-1)
	    iner = varsy[i]+surv3[i]*(1-surv3[i])*rvs[i]^2
	    denom = denom + out*iner
	    }


gen_length2 <- function(fec,surv,lambda){
    ells = cumprod(c(1,surv))
    nn = length(fec)
    T = 0.0
    for(i in 1:nn){
	T = T + i*lambda^(-i)*ells[i]*fec[i]
    }
    T
}

	  gl2 = gen_length2(fec,surv,popPar$growth_rate)
	  crowNe <- function(N,kbar,Vk,L){
	    nume = (2*N-1)*kbar
	    denom = 2*(1+Vk/kbar)
	    nume*L/denom
	  }

	  tN = tN[41:51,]
NeEstSim = (4*tN[,1]*(gl2))/(estvar+estmean)
cNe = crowNe(tN[,1],estmean,estvar,gl2)
cNeNOL = crowNe(tN[,1],estmean,estvar,1)
	realVk = condvar(dprob,fec*2,popPar$theta)
	realKbar = sum(cumsum(2*fec)*dprob)

	    library(tidyverse)



	  NestSC = data.frame(log_N=ssdrSC[rownames(ssdrSC) == "log_N",1],stdErr =ssdrSC[rownames(ssdrSC) == "log_N",2])
	NestSC$Llower = NestSC$log_N - qnorm(0.05/2,lower.tail = FALSE)*NestSC$stdErr
	NestSC$Lupper = NestSC$log_N + qnorm(0.05/2,lower.tail = FALSE)*NestSC$stdErr
	NestSC$lower = exp(NestSC$Llower)
	NestSC$upper = exp(NestSC$Lupper)
	NestSC$NestSC = exp(NestSC$log_N)
	NestSC$Ntrue = as.vector(t(tN))
	NestSC$age = rep(c(1,2,3),11)

	MestSC = data.frame(log_M=ssdrSC[rownames(ssdrSC) == "log_M",1],stdErr =ssdrSC[rownames(ssdrSC) == "log_M",2])
	MestSC$Llower = MestSC$log_M - qnorm(0.05/2,lower.tail = FALSE)*MestSC$stdErr
	MestSC$Lupper = MestSC$log_M + qnorm(0.05/2,lower.tail = FALSE)*MestSC$stdErr
	MestSC$lower = exp(MestSC$Llower)
	MestSC$upper = exp(MestSC$Lupper)
	MestSC$MestSC = exp(MestSC$log_M)
	MestSC$Mtrue = -log(c(popPar$surv1,popPar$surv2))
	MestSC$age = c(1,2)

	FecEstSC = data.frame(log_fec=ssdrSC[rownames(ssdrSC) == "log_fecundity",1],stdErr =ssdrSC[rownames(ssdrSC) == "log_fecundity",2])
	FecEstSC$Llower = FecEstSC$log_fec - qnorm(0.05/2,lower.tail = FALSE)*FecEstSC$stdErr
	FecEstSC$Lupper = FecEstSC$log_fec + qnorm(0.05/2,lower.tail = FALSE)*FecEstSC$stdErr
	FecEstSC$lower = exp(FecEstSC$Llower)
	FecEstSC$upper = exp(FecEstSC$Lupper)
	FecEstSC$FecEstSC = exp(FecEstSC$log_fec)
	FecEstSC$FecTrue = c(popPar$fec1,popPar$fec2,popPar$fec3)

	growth_rateEstSC = data.frame(estSCgrowth_rate = ssdrSC[rownames(ssdrSC) == "growth_rate",1],stdErr = ssdrSC[rownames(ssdrSC) == "growth_rate",2])
	growth_rateEstSC$lower = growth_rateEstSC$estSCgrowth_rate - qnorm(0.05/2,lower.tail = FALSE)*growth_rateEstSC$stdErr
	growth_rateEstSC$upper = growth_rateEstSC$estSCgrowth_rate + qnorm(0.05/2,lower.tail = FALSE)*growth_rateEstSC$stdErr
	growth_rateEstSC$trueGR = popPar$growth_rate


	NestMB = data.frame(log_N=ssdrMB[rownames(ssdrMB) == "log_N",1],stdErr =ssdrMB[rownames(ssdrMB) == "log_N",2])
	NestMB$Llower = NestMB$log_N - qnorm(0.05/2,lower.tail = FALSE)*NestMB$stdErr
	NestMB$Lupper = NestMB$log_N + qnorm(0.05/2,lower.tail = FALSE)*NestMB$stdErr
	NestMB$lower = exp(NestMB$Llower)
	NestMB$upper = exp(NestMB$Lupper)
	NestMB$NestMB = exp(NestMB$log_N)
	NestMB$Ntrue = as.vector(t(tN))
	NestMB$age = rep(c(1,2,3),11)

	MestMB = data.frame(log_M=ssdrMB[rownames(ssdrMB) == "log_M",1],stdErr =ssdrMB[rownames(ssdrMB) == "log_M",2])
	MestMB$Llower = MestMB$log_M - qnorm(0.05/2,lower.tail = FALSE)*MestMB$stdErr
	MestMB$Lupper = MestMB$log_M + qnorm(0.05/2,lower.tail = FALSE)*MestMB$stdErr
	MestMB$lower = exp(MestMB$Llower)
	MestMB$upper = exp(MestMB$Lupper)
	MestMB$MestMB = exp(MestMB$log_M)
	MestMB$Mtrue = -log(c(popPar$surv1,popPar$surv2))
	MestMB$age = c(1,2)

	FecEstMB = data.frame(log_fec=ssdrMB[rownames(ssdrMB) == "log_fecundity",1],stdErr =ssdrMB[rownames(ssdrMB) == "log_fecundity",2])
	FecEstMB$Llower = FecEstMB$log_fec - qnorm(0.05/2,lower.tail = FALSE)*FecEstMB$stdErr
	FecEstMB$Lupper = FecEstMB$log_fec + qnorm(0.05/2,lower.tail = FALSE)*FecEstMB$stdErr
	FecEstMB$lower = exp(FecEstMB$Llower)
	FecEstMB$upper = exp(FecEstMB$Lupper)
	FecEstMB$FecEstMB = exp(FecEstMB$log_fec)
	FecEstMB$FecTrue = c(popPar$fec1,popPar$fec2,popPar$fec3)

	thetaEstMB = data.frame(log_theta=ssdrMB[rownames(ssdrMB) == "log_theta",1],stdErr =ssdrMB[rownames(ssdrMB) == "log_theta",2])
	thetaEstMB$Llower = thetaEstMB$log_theta - qnorm(0.05/2,lower.tail = FALSE)*thetaEstMB$stdErr
	thetaEstMB$Lupper = thetaEstMB$log_theta + qnorm(0.05/2,lower.tail = FALSE)*thetaEstMB$stdErr
	thetaEstMB$lower = exp(thetaEstMB$Llower)
	thetaEstMB$upper = exp(thetaEstMB$Lupper)
	thetaEstMB$thetaEstMB = exp(thetaEstMB$log_theta)
	thetaEstMB$thetaTrue = c(popPar$theta)

	VkEstMB = data.frame(estVk = ssdrMB[rownames(ssdrMB) == "Vk",1],stdErr = ssdrMB[rownames(ssdrMB) == "Vk",2])
	VkEstMB$lower = VkEstMB$estVk - qnorm(0.05/2,lower.tail = FALSE)*VkEstMB$stdErr
	VkEstMB$upper = VkEstMB$estVk + qnorm(0.05/2,lower.tail = FALSE)*VkEstMB$stdErr
	VkEstMB$fromSim = estvar
	VkEstMB$fromSimPar = realVk

	growth_rateEstMB = data.frame(estgrowth_rate = ssdrMB[rownames(ssdrMB) == "growth_rate",1],stdErr = ssdrMB[rownames(ssdrMB) == "growth_rate",2])
	growth_rateEstMB$lower = growth_rateEstMB$estgrowth_rate - qnorm(0.05/2,lower.tail = FALSE)*growth_rateEstMB$stdErr
	growth_rateEstMB$upper = growth_rateEstMB$estgrowth_rate + qnorm(0.05/2,lower.tail = FALSE)*growth_rateEstMB$stdErr
	growth_rateEstMB$trueGR = popPar$growth_rate

	N_eEstMB = data.frame(estN_e = ssdrMB[rownames(ssdrMB) == "N_e",1],stdErr = ssdrMB[rownames(ssdrMB) == "N_e",2])
	N_eEstMB$lower = N_eEstMB$estN_e - qnorm(0.05/2,lower.tail = FALSE)*N_eEstMB$stdErr
	N_eEstMB$upper = N_eEstMB$estN_e + qnorm(0.05/2,lower.tail = FALSE)*N_eEstMB$stdErr
	  N_eEstMB$fromSim = cNe
	  #N_eEstMB$LDNE = ldNeEst$Neadj


	  ##Sensitivity Analysis Additions
	  ##For 2
		NestMVP2 = data.frame(log_N=ssdrMVP2[rownames(ssdrMVP2) == "log_N",1],stdErr =ssdrMVP2[rownames(ssdrMVP2) == "log_N",2])
	NestMVP2$Llower = NestMVP2$log_N - qnorm(0.05/2,lower.tail = FALSE)*NestMVP2$stdErr
	NestMVP2$Lupper = NestMVP2$log_N + qnorm(0.05/2,lower.tail = FALSE)*NestMVP2$stdErr
	NestMVP2$lower = exp(NestMVP2$Llower)
	NestMVP2$upper = exp(NestMVP2$Lupper)
	NestMVP2$NestMVP2 = exp(NestMVP2$log_N)
	NestMVP2$Ntrue = as.vector(t(tN))
	NestMVP2$age = rep(c(1,2,3),11)

	MestMVP2 = data.frame(log_M=ssdrMVP2[rownames(ssdrMVP2) == "log_M",1],stdErr =ssdrMVP2[rownames(ssdrMVP2) == "log_M",2])
	MestMVP2$Llower = MestMVP2$log_M - qnorm(0.05/2,lower.tail = FALSE)*MestMVP2$stdErr
	MestMVP2$Lupper = MestMVP2$log_M + qnorm(0.05/2,lower.tail = FALSE)*MestMVP2$stdErr
	MestMVP2$lower = exp(MestMVP2$Llower)
	MestMVP2$upper = exp(MestMVP2$Lupper)
	MestMVP2$MestMVP2 = exp(MestMVP2$log_M)
	MestMVP2$Mtrue = -log(c(popPar$surv1,popPar$surv2))
	MestMVP2$age = c(1,2)

	FecEstMVP2 = data.frame(log_fec=ssdrMVP2[rownames(ssdrMVP2) == "log_fecundity",1],stdErr =ssdrMVP2[rownames(ssdrMVP2) == "log_fecundity",2])
	FecEstMVP2$Llower = FecEstMVP2$log_fec - qnorm(0.05/2,lower.tail = FALSE)*FecEstMVP2$stdErr
	FecEstMVP2$Lupper = FecEstMVP2$log_fec + qnorm(0.05/2,lower.tail = FALSE)*FecEstMVP2$stdErr
	FecEstMVP2$lower = exp(FecEstMVP2$Llower)
	FecEstMVP2$upper = exp(FecEstMVP2$Lupper)
	FecEstMVP2$FecEstMVP2 = exp(FecEstMVP2$log_fec)
	FecEstMVP2$FecTrue = c(popPar$fec1,popPar$fec2,popPar$fec3)

	thetaEstMVP2 = data.frame(log_theta=ssdrMVP2[rownames(ssdrMVP2) == "log_theta",1],stdErr =ssdrMVP2[rownames(ssdrMVP2) == "log_theta",2])
	thetaEstMVP2$Llower = thetaEstMVP2$log_theta - qnorm(0.05/2,lower.tail = FALSE)*thetaEstMVP2$stdErr
	thetaEstMVP2$Lupper = thetaEstMVP2$log_theta + qnorm(0.05/2,lower.tail = FALSE)*thetaEstMVP2$stdErr
	thetaEstMVP2$lower = exp(thetaEstMVP2$Llower)
	thetaEstMVP2$upper = exp(thetaEstMVP2$Lupper)
	thetaEstMVP2$thetaEstMVP2 = exp(thetaEstMVP2$log_theta)
	thetaEstMVP2$thetaTrue = c(popPar$theta)

	VkEstMVP2 = data.frame(estVk = ssdrMVP2[rownames(ssdrMVP2) == "Vk",1],stdErr = ssdrMVP2[rownames(ssdrMVP2) == "Vk",2])
	VkEstMVP2$lower = VkEstMVP2$estVk - qnorm(0.05/2,lower.tail = FALSE)*VkEstMVP2$stdErr
	VkEstMVP2$upper = VkEstMVP2$estVk + qnorm(0.05/2,lower.tail = FALSE)*VkEstMVP2$stdErr
	VkEstMVP2$fromSim = estvar
	VkEstMVP2$fromSimPar = realVk

	growth_rateEstMVP2 = data.frame(estgrowth_rate = ssdrMVP2[rownames(ssdrMVP2) == "growth_rate",1],stdErr = ssdrMVP2[rownames(ssdrMVP2) == "growth_rate",2])
	growth_rateEstMVP2$lower = growth_rateEstMVP2$estgrowth_rate - qnorm(0.05/2,lower.tail = FALSE)*growth_rateEstMVP2$stdErr
	growth_rateEstMVP2$upper = growth_rateEstMVP2$estgrowth_rate + qnorm(0.05/2,lower.tail = FALSE)*growth_rateEstMVP2$stdErr
	growth_rateEstMVP2$trueGR = popPar$growth_rate

	N_eEstMVP2 = data.frame(estN_e = ssdrMVP2[rownames(ssdrMVP2) == "N_e",1],stdErr = ssdrMVP2[rownames(ssdrMVP2) == "N_e",2])
	N_eEstMVP2$lower = N_eEstMVP2$estN_e - qnorm(0.05/2,lower.tail = FALSE)*N_eEstMVP2$stdErr
	N_eEstMVP2$upper = N_eEstMVP2$estN_e + qnorm(0.05/2,lower.tail = FALSE)*N_eEstMVP2$stdErr
	  N_eEstMVP2$fromSim = cNe

	  ##For 1
		NestMVP1 = data.frame(log_N=ssdrMVP1[rownames(ssdrMVP1) == "log_N",1],stdErr =ssdrMVP1[rownames(ssdrMVP1) == "log_N",2])
	NestMVP1$Llower = NestMVP1$log_N - qnorm(0.05/2,lower.tail = FALSE)*NestMVP1$stdErr
	NestMVP1$Lupper = NestMVP1$log_N + qnorm(0.05/2,lower.tail = FALSE)*NestMVP1$stdErr
	NestMVP1$lower = exp(NestMVP1$Llower)
	NestMVP1$upper = exp(NestMVP1$Lupper)
	NestMVP1$NestMVP1 = exp(NestMVP1$log_N)
	NestMVP1$Ntrue = as.vector(t(tN))
	NestMVP1$age = rep(c(1,2,3),11)

	MestMVP1 = data.frame(log_M=ssdrMVP1[rownames(ssdrMVP1) == "log_M",1],stdErr =ssdrMVP1[rownames(ssdrMVP1) == "log_M",2])
	MestMVP1$Llower = MestMVP1$log_M - qnorm(0.05/2,lower.tail = FALSE)*MestMVP1$stdErr
	MestMVP1$Lupper = MestMVP1$log_M + qnorm(0.05/2,lower.tail = FALSE)*MestMVP1$stdErr
	MestMVP1$lower = exp(MestMVP1$Llower)
	MestMVP1$upper = exp(MestMVP1$Lupper)
	MestMVP1$MestMVP1 = exp(MestMVP1$log_M)
	MestMVP1$Mtrue = -log(c(popPar$surv1,popPar$surv2))
	MestMVP1$age = c(1,2)

	FecEstMVP1 = data.frame(log_fec=ssdrMVP1[rownames(ssdrMVP1) == "log_fecundity",1],stdErr =ssdrMVP1[rownames(ssdrMVP1) == "log_fecundity",2])
	FecEstMVP1$Llower = FecEstMVP1$log_fec - qnorm(0.05/2,lower.tail = FALSE)*FecEstMVP1$stdErr
	FecEstMVP1$Lupper = FecEstMVP1$log_fec + qnorm(0.05/2,lower.tail = FALSE)*FecEstMVP1$stdErr
	FecEstMVP1$lower = exp(FecEstMVP1$Llower)
	FecEstMVP1$upper = exp(FecEstMVP1$Lupper)
	FecEstMVP1$FecEstMVP1 = exp(FecEstMVP1$log_fec)
	FecEstMVP1$FecTrue = c(popPar$fec1,popPar$fec2,popPar$fec3)

	thetaEstMVP1 = data.frame(log_theta=ssdrMVP1[rownames(ssdrMVP1) == "log_theta",1],stdErr =ssdrMVP1[rownames(ssdrMVP1) == "log_theta",2])
	thetaEstMVP1$Llower = thetaEstMVP1$log_theta - qnorm(0.05/2,lower.tail = FALSE)*thetaEstMVP1$stdErr
	thetaEstMVP1$Lupper = thetaEstMVP1$log_theta + qnorm(0.05/2,lower.tail = FALSE)*thetaEstMVP1$stdErr
	thetaEstMVP1$lower = exp(thetaEstMVP1$Llower)
	thetaEstMVP1$upper = exp(thetaEstMVP1$Lupper)
	thetaEstMVP1$thetaEstMVP1 = exp(thetaEstMVP1$log_theta)
	thetaEstMVP1$thetaTrue = c(popPar$theta)

	VkEstMVP1 = data.frame(estVk = ssdrMVP1[rownames(ssdrMVP1) == "Vk",1],stdErr = ssdrMVP1[rownames(ssdrMVP1) == "Vk",2])
	VkEstMVP1$lower = VkEstMVP1$estVk - qnorm(0.05/2,lower.tail = FALSE)*VkEstMVP1$stdErr
	VkEstMVP1$upper = VkEstMVP1$estVk + qnorm(0.05/2,lower.tail = FALSE)*VkEstMVP1$stdErr
	VkEstMVP1$fromSim = estvar
	VkEstMVP1$fromSimPar = realVk

	growth_rateEstMVP1 = data.frame(estgrowth_rate = ssdrMVP1[rownames(ssdrMVP1) == "growth_rate",1],stdErr = ssdrMVP1[rownames(ssdrMVP1) == "growth_rate",2])
	growth_rateEstMVP1$lower = growth_rateEstMVP1$estgrowth_rate - qnorm(0.05/2,lower.tail = FALSE)*growth_rateEstMVP1$stdErr
	growth_rateEstMVP1$upper = growth_rateEstMVP1$estgrowth_rate + qnorm(0.05/2,lower.tail = FALSE)*growth_rateEstMVP1$stdErr
	growth_rateEstMVP1$trueGR = popPar$growth_rate

	N_eEstMVP1 = data.frame(estN_e = ssdrMVP1[rownames(ssdrMVP1) == "N_e",1],stdErr = ssdrMVP1[rownames(ssdrMVP1) == "N_e",2])
	N_eEstMVP1$lower = N_eEstMVP1$estN_e - qnorm(0.05/2,lower.tail = FALSE)*N_eEstMVP1$stdErr
	N_eEstMVP1$upper = N_eEstMVP1$estN_e + qnorm(0.05/2,lower.tail = FALSE)*N_eEstMVP1$stdErr
	  N_eEstMVP1$fromSim = cNe



	ret = list(optSC=optSC,convergedSC=optSC$convergence,reportSC=reportSC,NestSC=NestSC,MestSC=MestSC,FecEstSC=FecEstSC,
		   growth_rateEstSC=growth_rateEstSC,
		   optMB=optMB,convergedMB=optMB$convergence,reportMB=reportMB,NestMB=NestMB,MestMB=MestMB,FecEstMB=FecEstMB,thetaEstMB=thetaEstMB,VkEstMB=VkEstMB,
		   growth_rateEstMB=growth_rateEstMB,N_eEstMB=N_eEstMB,
		    optMVP2=optMVP2,convergedMVP2=optMVP2$convergence,reportMVP2=reportMVP2,NestMVP2=NestMVP2,MestMVP2=MestMVP2,FecEstMVP2=FecEstMVP2,thetaEstMVP2=thetaEstMVP2,VkEstMVP2=VkEstMVP2,
		   growth_rateEstMVP2=growth_rateEstMVP2,N_eEstMVP2=N_eEstMVP2,
		    optMVP1=optMVP1,convergedMVP1=optMVP1$convergence,reportMVP1=reportMVP1,NestMVP1=NestMVP1,MestMVP1=MestMVP1,FecEstMVP1=FecEstMVP1,thetaEstMVP1=thetaEstMVP1,VkEstMVP1=VkEstMVP1,
		   growth_rateEstMVP1=growth_rateEstMVP1,N_eEstMVP1=N_eEstMVP1)


	  fitteddir = paste0(rep_df[[k]]$foldname[1],"/fitted/")
	  if(!dir.exists(fitteddir)){
	    dir.create(fitteddir)
	  }

	  saveRDS(ret,file=paste0(fitteddir,"fittedR",str_pad(rr,4,pad="0"),".rds"))
	}
      }
