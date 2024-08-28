##This just reruns model with N_e to new form +estimate.

repdf <- "../sims/rep_df1.rds"
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



dyn.load(dynlib("../models/ssmodelAdj"))
dyn.load(dynlib("../models/ssmodelNoSC"))


rep_df = readRDS(repdf)
LDNEQC = readRDS("../LDNEQuantilesC.rds")
NESC = list()

for(k in 1:length(rep_df)){

    files = list.files(paste0(rep_df[[k]]$foldname[1]),pattern = "*.rds",full.names = TRUE)
    popbit =  regexec("pop[0-9][0-9][0-9][0-9]",files)
    mlength = sapply(popbit,function(x){attributes(x)$match.length})
    popID = unique(do.call(substring,list(text=files,first=popbit,last=unlist(popbit)+mlength-1)))

    indiv = readRDS(paste0("../populations/",popID,".rds"))
    ped = make_ped(indiv)

    popPar = indiv$popPars
    tN = n_at_age(indiv)[,]
    tN = tN[41:51,]
    NESCIn = list()

    popLDNE50 = LDNEQC |>
        filter(population == popID) |>
        filter(q == 0.5)

    popLDNE5 = LDNEQC |>
        filter(population == popID) |>
        filter(q == 0.05)

    popLDNE95 = LDNEQC |>
        filter(population == popID) |>
        filter(q == 0.95)

    
    for(rr in 1:length(files)){
    ##for(rr in 1:2){
        print(rr)
        stuff = readRDS(files[rr])

        dat = stuff$data
        dat$skipsame = 0
        parm = stuff$parm
        parm$log_init_tot_N = 10

	mapp2 = list(log_M=as.factor(c(1,2,NA)),log_init_tot_N=as.factor(c(1)),
		     log_fecundity=as.factor(c(1,2,3)),log_theta=as.factor(c(1)))


        datMB = stuff$data
        datMB$skipsame = 0
        parmMB = stuff$parm

        objMB = MakeADFun(datMB,parmMB,map=mapp2,DLL="ssmodelAdj")
        optMB = nlminb(objMB$par,objMB$fn,objMB$gr,control=list(iter.max=5000))
        reportMB = objMB$report()
        sdrMB = sdreport(objMB)
        ssdrMB = summary(sdrMB)


        dprob = c(1-popPar$surv1,popPar$surv1*(1-popPar$surv2),popPar$surv2*popPar$surv1)
        fec = c(popPar$fec1,popPar$fec2,popPar$fec3)

                  engenNe <- function(N,surv,ells,vars,repro_value,growth_rate,glen){
              Ne_engen = numeric(length(N))
              NeE_denom = 0
              for(a in 1:length(ells)){
                  front = ells[a]*growth_rate^(-a)
                  back = vars[a]+surv[a]*(1-surv[a])*repro_value[a+1]^2
                  NeE_denom = NeE_denom + front*back
              }
              for(y in 1:length(N)){
                  top = glen*N[y]
                  Ne_engen[y] = top/NeE_denom
              }
              Ne_engen
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

	fec = c(popPar$fec1,popPar$fec2,popPar$fec3)
        surv = c(popPar$surv1,popPar$surv2)
        ells = c(1,popPar$surv1,popPar$surv2*popPar$surv1)  
	les = make_les(fec,surv)

        
                  
        repro_value = Re(eigen(t(les))$vectors[,1]/eigen(t(les))$vectors[1,1])
        	dprob = c(1-popPar$surv1,popPar$surv1*(1-popPar$surv2),popPar$surv2*popPar$surv1)
	  doub_fec = 2*fec
	  thet = popPar$theta
	  varsy = doub_fec+doub_fec^2/thet


	  growth_rate = max(Re(eigen(les,only.values = TRUE)$values))
	  surv2 = cumprod(c(1,surv,0))
	  surv3 = c(surv,0)
	  fec2 = c(fec,0)

        

        gl2 = gen_length2(fec,surv,popPar$growth_rate)

      eNe  =  engenNe(tN[,1],surv3,ells,varsy,c(repro_value,0),growth_rate,gl2)        
        modelNe = ssdrMB[rownames(ssdrMB) == "Ne_engen",]
        NeDF = data.frame(sim_Ne=eNe,modelNe=reportMB$Ne_engen)
        NeDF$upper = NeDF$modelNe+1.96*modelNe[,2]
        NeDF$lower = NeDF$modelNe-1.96*modelNe[,2]
        NeDF$year = 1:11
        NeDF$sim = rr
        NeDF$pop = popID
        NeDF$convergence = optMB$convergence
        NeDF$LD5 = popLDNE5$Ne
        NeDF$LD50 = popLDNE50$Ne
        NeDF$LD95 = popLDNE95$Ne
        NESCIn[[rr]] = NeDF
    }
    NESC[[k]] = do.call(rbind,NESCIn)
}
NESCcom = do.call(rbind,NESC)
saveRDS(NESCcom,"../NESCcom1B.rds")
