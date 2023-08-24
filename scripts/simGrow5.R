repdf <- "../sims/rep_df5.rds"
 library(tidyverse)
  library(stringr)
  library(varCKMR)

  sim_df = read.csv("../sims2/sim_df.csv")
  sim_loc = "../"

  rep_df = readRDS(repdf)

  nyears = 50
  sampyears = 10
  sim1y = nyears-sampyears


  for(k in 1:length(rep_df)){

    indiv = readRDS(paste0(sim_loc,"populations/","pop",str_pad(rep_df[[k]]$X[1],4,pad="0"),".rds"))
ped = make_ped(indiv)
pp = ped$pedigree

    for(rr in 1:nrow(rep_df[[k]])){
      print(rr)
      set.seed(rep_df[[k]]$seed[rr])
      library(tidyverse)

      nyears = 50
      sim_years = 41:50

      tN = t(n_at_age(indiv)[sim_years,])

      ##Get a sample some multiple of true size
      fun_fun = function(multiplier){
	fun = function(y){
	  multiplier*sqrt(y)
	}
	fun
      }

      fun = fun_fun(2.5)

      samp = rough_sample_size(tN,fun)
      samps = sample_age_by_year(indiv$population,samp,sim_years,3)

      library(tidyverse)
      library(reshape)

      ##Try and find POPs and HSPs
      sampC = do.call(rbind,unlist(samps,FALSE))
      sampC1 = sampC[,1:8]
      sampC2 = sampC[,1:8]

      sampCDF1 = sampC |>
	group_by(birth_year,cur_year) |>
	summarise(n_samp=n())

      sampCDF2 = sampCDF1
      names(sampCDF1) = paste0(names(sampCDF1),".1")
      names(sampCDF2) = paste0(names(sampCDF2),".2")
      sampCDFF = expand.grid.df(sampCDF1,sampCDF2)

      sampCDFF = sampCDFF |>
	filter(birth_year.1 > birth_year.2) |>
	filter(birth_year.1 - birth_year.2 <= 3) |>
	mutate(n_comp=n_samp.1*n_samp.2) 

      ##Find POPs quicker
      fath = data.frame(ID=sampC$father,row=1:nrow(sampC))
      moth = data.frame(ID=sampC$mother,row=1:nrow(sampC))
      ID = data.frame(ID=sampC$ID,col=1:nrow(sampC))
      FOPsQ = inner_join(fath,ID)
      MOPsQ = inner_join(moth,ID)
      FOPs = cbind(sampC[FOPsQ[,"row"],c("birth_year","cur_year")],sampC[FOPsQ[,"col"],c("birth_year","cur_year")])
      MOPs = cbind(sampC[MOPsQ[,"row"],c("birth_year","cur_year")],sampC[MOPsQ[,"col"],c("birth_year","cur_year")])
      POPsC = rbind(FOPs,MOPs) 
      names(POPsC) = c("birth_year.1","cur_year.1","birth_year.2","cur_year.2")
      POPsC = POPsC |>
	group_by(birth_year.1,cur_year.1,birth_year.2,cur_year.2) |>
	summarize(POPs=n())
      POPsC2 = POPsC

      POPs = left_join(sampCDFF,POPsC)
      POPs$POPs[is.na(POPs$POPs)] = 0
      POPs = filter(POPs,birth_year.1 > (sim1y -1))
      POPs = select(POPs,birth_year.1,cur_year.1,birth_year.2,cur_year.2,n_comp,POPs)
      POPs2 = POPs

      ##Find Sibs quicker
      sampCDFF2A = expand.grid.df(sampCDF1,sampCDF2)
      sampCDFF2 = sampCDFF2A |>
	filter(birth_year.1 >= birth_year.2) |>
	filter(birth_year.1 - birth_year.2 <= 2) |>
	##Doesn't work for SC
	mutate(n_comp=ifelse(birth_year.1 == birth_year.2 & cur_year.1 == cur_year.2,choose(n_samp.1,2),n_samp.1*n_samp.2)) |>
	select(birth_year.1,cur_year.1,birth_year.2,cur_year.2,n_comp)

      bbs1 = data.frame(birth_year=sampC$birth_year,row=1:nrow(sampC))
      bbs2 = data.frame(birth_year=sampC$birth_year,col=1:nrow(sampC))
      Births = inner_join(bbs1,bbs2) |>
	filter(row !=col)

      BirthsSC = cbind(sampC[Births[,"row"],c("birth_year","cur_year","ID","mother","father")],sampC[Births[,"col"],c("birth_year","cur_year","ID","mother","father")])
      names(BirthsSC) = c("birth_year.1","cur_year.1","ID.1","mother.1","father.1","birth_year.2","cur_year.2","ID.2","mother.2","father.2")
      BirthsSC = BirthsSC |>
	filter(ID.1 != ID.2) |>
	filter(birth_year.1 >= birth_year.2) |>
	distinct()

      sampCDFF3 = BirthsSC |>
	group_by(birth_year.1,cur_year.1,birth_year.2,cur_year.2) |>
	summarize(n_comp=n())

      fath1 = data.frame(ID=sampC$father,row=1:nrow(sampC))
      moth1 = data.frame(ID=sampC$mother,row=1:nrow(sampC))
      fath2 = data.frame(ID=sampC$father,col=1:nrow(sampC))
      moth2 = data.frame(ID=sampC$mother,col=1:nrow(sampC))


      FSPsArr = inner_join(fath1,fath2) |>
	filter(row != col)
      MSPsArr = inner_join(moth1,moth2) |>
	filter(row != col)
      FSPs = cbind(sampC[FSPsArr[,"row"],c("birth_year","cur_year","ID","mother","father")],sampC[FSPsArr[,"col"],c("birth_year","cur_year","ID","mother","father")])
      MSPs = cbind(sampC[MSPsArr[,"row"],c("birth_year","cur_year","ID","mother","father")],sampC[MSPsArr[,"col"],c("birth_year","cur_year","ID","mother","father")])
      SPsA = rbind(FSPs,MSPs)
      names(SPsA) = c("birth_year.1","cur_year.1","ID.1","mother.1","father.1","birth_year.2","cur_year.2","ID.2","mother.2","father.2")
      SPs = SPsA |>
	filter(ID.1 != ID.2) |>
	filter(birth_year.1 >= birth_year.2) |>
	distinct()

      SPsC = SPs |>
	group_by(birth_year.1,cur_year.1,birth_year.2,cur_year.2) |>
	summarize(Sibs=n())


      SPsBP = SPsA |>
	filter(ID.1 != ID.2) |>
	filter(birth_year.1 == birth_year.2) |>
	distinct()

      SPsBP$tups = apply(SPsBP, 1, function(x) {
	paste(sort(c(x[3], x[8])), collapse = "|")
      })

      SPsBP = SPsBP |>
	distinct(tups,.keep_all=TRUE) |>
	group_by(birth_year.1) |>
	summarise(nsib = n()) |>
	select(birth_year=birth_year.1,nsib)


      SibsC = left_join(sampCDFF2,SPsC)
      SibsC$Sibs[is.na(SibsC$Sibs)] = 0
      SibsC = filter(SibsC,birth_year.1 > (sim1y -1))
      SibsC = filter(SibsC,birth_year.2 > (sim1y -1))

      SibsXC = SibsC |>
	filter(birth_year.1 != birth_year.2)


      SibsXC2 = SibsXC
##      SibsSCP2 = SibsSCP

      ##Same Cohort, using clusters
      SibtyComps = sampC |>
	distinct(ID,.keep_all=TRUE) |>
	group_by(birth_year) |>
	summarize(n_draws=n()*2)

      SibtySCM = sampC |>
	distinct(ID,.keep_all=TRUE) |>
	group_by(birth_year,mother) |>
	summarize(n_MSB=n()) |>
	filter(n_MSB > 1) |>
	group_by(birth_year) |>
	summarise(n_MSB=sum(n_MSB))

      SibtySCF = sampC |>
	distinct(ID,.keep_all=TRUE) |> 
	group_by(birth_year,father) |>
	summarize(n_FSB=n()) |>
	filter(n_FSB > 1) |>
	group_by(birth_year) |>
	summarize(n_FSB=sum(n_FSB))

      SibtySC = left_join(SibtySCM,SibtySCF)
      SibtySC = left_join(SibtySC,SibtyComps)
      SibsSC = SibtySC |>
	mutate(n_SB = n_MSB+n_FSB) |>
	select(birth_year,n_SB,n_draws) |>
	filter(birth_year > (sim1y -1))


sNcomp = sampC1 |>
    group_by(birth_year) |>
    summarize(ncomp=n()) |>
    mutate(ncomp=choose(ncomp,2))

  SCPgorp2 = sNcomp |>
    left_join(SPsBP) |>
    mutate(aprob = nsib/ncomp)

  SCP = SCPgorp2 |>
    filter(birth_year >= sim1y) |>
    mutate(birth_year = birth_year-sim1y) |>
    select(birth_year,nsib,ncomp)



      ##Correct bounds
      SibsSC[,1] = SibsSC[,1]-sim1y
      SibsXC[,1:4] = SibsXC[,1:4]-sim1y
      POPs[,1:4] = POPs[,1:4]-sim1y

      data = list()
      data$Y = 10+1
      data$A = 3
      data$POPs = as.matrix(POPs)
      data$SIBsXC = as.matrix(SibsXC)
      data$SIBsSC = as.matrix(SibsSC)
      data$SIBsSCP = as.matrix(SCP)


      parm = list()
      parm$log_init_tot_N = log(3000)
      parm$log_M = log(c(0.2,0.2,0))
      parm$log_fecundity = log(c(1,1,1))
      parm$log_theta = log(1)




      ret = list(data=data,parm=parm)
      saveRDS(ret,file=paste0(rep_df[[k]]$foldname[1],"/","clean",str_pad(rr,4,pad="0"),".rds"))
    }
  }
