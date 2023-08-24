library(varCKMR)
  library(stringr)

  set.seed(1234)
  thetas = c(0.1,0.25,0.5,0.75,1,1.25,2,5)
  growth_rate = c(0.95,1,1.01)
  pops = expand.grid(theta=thetas,growth_rate=growth_rate)

  sim_seeds = sample(1:100000,nrow(pops))

  sim_df = data.frame(seed=sim_seeds)
  sim_df$growth_rate = pops$growth_rate
  sim_df$theta = pops$theta
  sim_df$fec1 = 0
  sim_df$fec2 = 0
  sim_df$fec3 = 0
  sim_df$surv1 = 0
  sim_df$surv2 = 0


  target = numeric(length(sim_seeds))
  for(i in 1:nrow(sim_df)){
    set.seed(sim_df$seed[i])
    target[i] = sim_df$growth_rate[i]
    gr = -1
  while(!isTRUE(all.equal(gr,target[i]))){
      print(paste0("start",gr))
      parms = get_pars_with_target(target[i])
      gr = parms$growth_rate
      print(gr)
  }
    sim_df$growth_rate[i] = parms$growth_rate
    sim_df$fec1[i] = parms$fec[1]
    sim_df$fec2[i] = parms$fec[2]
    sim_df$fec3[i] = parms$fec[3]
    sim_df$surv1[i] = parms$surv[1]
    sim_df$surv2[i] = parms$surv[2]
  }

if(!dir.exists("../populations/")){
  dir.create("../populations/")
}

if(!dir.exists("../sims/")){
  dir.create("../sims/")
}



  write.csv(sim_df,file="../sims/sim_df.csv")

  sim_loc = "../populations/"

  ## Now to actually generate the simulations
  library(varCKMR)
  library(tidyverse)

  nyears = 50
  sampyears = 10
  growth_rate = 1
  init_N = c(7000,2500,500)
  genoM=100

  sim1y = nyears-sampyears

  gts = rep(2,genoM)
  leper = 1:26
  genos = unlist(lapply(gts,function(x){leper[seq_len(x)]}))
  locus = numeric(0)
  freqs = numeric(0)
  for(i in 1:genoM){
    locus = c(locus,rep(i,gts[i]))
    freqs = c(freqs,c(0.50,0.50))
  }
  marks = data.frame(geno=genos,locus=locus,freqs = freqs)

  ##Functions to make functions
  surv_fun_fun_sex <- function(probs_m,probs_f){
    surv_fun <- function(age,sex){
      age = pmin(age,length(probs_m))
      ret = ifelse(sex == "male",probs_m[age],probs_f[age])
      ret
    }    
    surv_fun
  }

  mature_fun_fun <- function(probs){
    mature_fun <- function(age){
      age = pmin(age+1,length(probs))
      probs[age]
    }
    mature_fun
  }

  fecun_fun_fun_sex <- function(probs_m,probs_f){
    fecun_fun <- function(age,sex){
      age = pmin(age,length(probs_m))
      ret = ifelse(sex == "male",probs_m[age],probs_f[age])
      ret
    }
    fecun_fun
  }

  sex_fraction = 0.5


  N = data.frame(Age = rep(1:(length(init_N)),2),Sex=c(rep("male",length(init_N)),rep("female",length(init_N))),N=c(floor(init_N*sex_fraction),(init_N-init_N*sex_fraction)))

  for(i in 1:nrow(sim_df)){
    set.seed(sim_df$seed[i])
    targ = get_pars_with_target(1)
    theta = sim_df$theta[i]

    mature_fun <- mature_fun_fun(c(1,1,1,1))

    ageMort = surv_fun_fun_sex(c(targ$surv,0),c(targ$surv,0))
    ageFec = fecun_fun_fun_sex(targ$fec,targ$fec)

    indiv <- create_founding_pop(N,marks,mature_fun,ageMort,ageFec)

    for(y in 1:(sim1y-10)){
      indiv <- breed_one_year_od(indiv,y-1,c(theta,theta),sex_fraction)
    }

    fec = c(sim_df$fec1[i],sim_df$fec2[i],sim_df$fec3[i])
    surv = c(sim_df$surv1[i],sim_df$surv2[i],0)
    theta = sim_df$theta[i]

    ageMort2 = surv_fun_fun_sex(surv,surv)
    ageFec2 = fecun_fun_fun_sex(fec,fec)
    indiv$fecund_fun = ageFec2
    indiv$surv_fun = ageMort2


    for(y in (sim1y+1-10):nyears){
      indiv <- breed_one_year_od(indiv,y-1,c(theta,theta),sex_fraction)
    }

    indiv$popPars = sim_df[i,]

    saveRDS(indiv,file=paste0(sim_loc,"pop",str_pad(i,4,pad="0"),".rds"))
  }
