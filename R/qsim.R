#' Create a founding population
#'
#' @param init_N dataframe with the initial population details, must be Age, Sex then N
#' @param markers dataframe with information on markers must be geno, locus then freqs
#' @param maturity_fun The function that specifies at what age they become mature
#' @param surv_fun the survival function
#' @param fecund_fun the mean fecundity at age function
#'
#' @import data.table
#' @export
create_founding_pop <- function(init_N,markers,maturity_fun,surv_fun,fecund_fun){
    dt = data.table(init_N)
    smt = simulate_markers(markers,1)
    dt = dt[rep(seq_len(nrow(dt)),N)]
    dt[, `:=`(birth_year=0-Age,
              mother=NA,
              father=NA,
              cur_year = 0,
              alive = TRUE,
              u = runif(.N),
              ID=seq_len(nrow(dt)))]
    dt[,mature := (u <= maturity_fun(Age))]
    dt[,N:=NULL]
    dt[,Age:=NULL]
    dt[,u:=NULL]

    dt[,names(smt) := simulate_markers(markers,nrow(dt)),]

    smtlist = split(names(smt),ceiling(seq_along(names(smt))/2))

    ret = list(population= dt,markers=markers,init_N=init_N,
               maturity_fun=maturity_fun,surv_fun= surv_fun,fecund_fun= fecund_fun,
               smtlist = smtlist,namsmt=names(smt))
    
    ret
}

#' Death function
#'
#' Given a population and survival function choose who dies
#'
#' @param population the population to have individuals die in
#' @param surv_fun the survival function giving the probability of an individual making it to next year
#' @export
yearly_deaths <- function(population,year_now){

    setkey(population$population,cur_year,alive)
    
    new_year = data.table::copy(population$population[.(year_now-1,TRUE)])

    new_year[,cur_year := year_now]
    
    setkey(new_year,alive,cur_year)
    new_year[,age := cur_year - birth_year][.(TRUE,year_now),u := runif(.N)][.(TRUE,year_now),alive := (u <= population$surv_fun(age))]
    setkey(new_year,alive,mature,cur_year)
    new_year[.(TRUE,FALSE,year_now),w := runif(.N)][.(TRUE,FALSE,year_now),mature := (w <= population$maturity_fun(age))] 
    new_year[,age := NULL][,u := NULL][,w := NULL]

    population$population = rbind(population$population,new_year)
    population
    

}

#' Breed 1 year for females only
#'
#' Select mature inviduals breed them together to create new kids!
#'
#' @param population the population to breed
#' @export
breed_one_year_fonly <- function(population,year_now,n_babies){

      setkey(population$population,alive,mature,cur_year,Sex)
    population$population[,age := cur_year-birth_year]

    mums = population$population[.(TRUE,TRUE,year_now,"female"),]
    mums = mums[,nkids := n_babies(age+1,.N,population)]

    faths = population$population[.(TRUE,TRUE,year_now,"male"),]
    IDmax = max(population$population$ID)+1

    glue = c("nkids","ID",population$namsmt)
    
    new_kids = mums[,..glue]
    setnames(new_kids,"ID","mother")
    new_kids = new_kids[rep(seq_len(nrow(new_kids)),nkids)]
    new_kids[, `:=`(birth_year=year_now+1,
                    father=NA,
                    cur_year=year_now+1,
                    alive = TRUE,
                    mature = FALSE,
                    Sex = sample(c("female"),.N,TRUE),
                    ID = seq(IDmax,length.out = nrow(new_kids)))]
    mums[,nkids := NULL]
    population$population[,age := NULL]
    new_kids[,nkids := NULL]
    setcolorder(new_kids,names(population$population))

    population$population = rbind(population$population,new_kids)
    population

}

#' Deaths and then randomly breed 2 sexes
#' @param population  the population to die and breed in
#' @param year_now the current year of the simulation
#' @param n_babies function describing the distribution how offspring are generated
#' @export
deaths_breed_one_year <- function(population,year_now,n_babies){
    setkey(population$population,cur_year,alive)
    new_deaths = data.table::copy(population$population[.(year_now-1,TRUE)])

    new_deaths[,cur_year := year_now]


    setkey(new_deaths,alive,cur_year)
    new_deaths[,age := cur_year - birth_year][.(TRUE,year_now),u := runif(.N)][.(TRUE,year_now),alive := (u <= population$surv_fun(age))]
    setkey(new_deaths,alive,mature,cur_year)
    new_deaths[.(TRUE,FALSE,year_now),w := runif(.N)][.(TRUE,FALSE,year_now),mature := (w <= population$maturity_fun(age))] 
    new_deaths[,u := NULL][,w := NULL]

    setkey(new_deaths,alive,mature,cur_year,Sex)


    mums = new_deaths[.(TRUE,TRUE,year_now,"female"),]
    mums = mums[,nkids := n_babies(age+1,.N,new_deaths)*2]

    faths = new_deaths[.(TRUE,TRUE,year_now,"male"),]
    IDmax = max(population$population$ID)+1

    glue = c("nkids","ID",population$namsmt)

    new_kids = mums[,..glue]
    setnames(new_kids,"ID","mother")
    new_kids = new_kids[rep(seq_len(nrow(new_kids)),nkids)]
    new_kids[, `:=`(birth_year=year_now+1,
                    father=NA,
                    cur_year=year_now+1,
                    alive = TRUE,
                    mature = FALSE,
                    Sex = sample(c("male","female"),.N,TRUE),
                    ID = seq(IDmax,length.out = nrow(new_kids)))]
    mums[,nkids := NULL]
    population$population[,age := NULL]
    new_kids[,nkids := NULL]
    setcolorder(new_kids,names(population$population))


    new_kids[,father := sample(faths$ID,.N,replace=TRUE)]
    kids_markers = samp_markers(population$smtlist,nrow(new_kids))

    marks = c("ID",population$namsmt)
    setkey(population$population,ID)
    DNAs = unique(population$population[.(ID = new_kids$father),..marks])
    allDNAs = as.data.frame(merge(new_kids,DNAs,by.x="father",by.y="ID",suffixes = c("",".m")))
    allDNAs2 = allDNAs


    gloop = system.time({    picksL = as.data.frame(as.table(kids_markers))
        picksL$Var1 = as.numeric(picksL$Var1)
        picksL$Var2 = as.numeric(picksL$Var2)+min(new_kids$ID)-1
        names(picksL) = c("whocares","ID","marker")

        new_kidsL = melt(new_kids,measure.vars=population$namsmt,value.name="geno",variable.name="marker")
        new_kidsL$marker = as.character(new_kidsL$marker)
        DNAsL = melt(DNAs,measure.vars=population$namsmt,value.name="geno",variable.name="marker")
        names(DNAsL)[1] = "father"
        setkey(new_kidsL,ID,marker)
        grab = new_kidsL[picksL[,2:3],.(father,marker)]
        setkey(DNAsL,father,marker)
        new_kidsL[picksL[,2:3],"geno"] = DNAsL[grab,3]
        DNAsL[grab]
        new_kids = dcast(new_kidsL,Sex+birth_year+mother+father+cur_year+alive+ID+mature ~ marker,
                         value.var = "geno")
    })

    population$population = rbind(population$population,new_kids)
    population


}

#' Generate kids randomly with overdispersion
#' 
#' Generate the list of kids born this year
#' @param IDs vector of IDs
#' @param ages vector of ages
#' @param sex vector of individuals sex
#' @param ave_fecun the average fecundity at age (matrix one column for each sex)
#' @param overdispersion overdispersion parameter (one for each sex)
#'
#' @export
#'
generate_kids_od <- function(IDs,ages,sex,ave_fecun,overdispersion){
  df = data.frame(ID=IDs,sex=sex,age=ages,sex_num=as.numeric(as.factor(sex)))
  ave_numkids = table(df$sex,df$age)[1,]*ave_fecun[,1]
  numkids = rpois(3,ave_numkids)
  ##Trying this?
  tot_kids = sum(2*ave_numkids)
  sdf = df |>
    dplyr::group_by(sex_num) |>
    dplyr::mutate(overd = overdispersion[sex_num]) |>
    dplyr::mutate(E = rgamma(dplyr::n(),overd)) |>
    dplyr::mutate(fecun = ave_fecun[cbind(age,sex_num)]*E) |>
    dplyr::mutate(nkids = rmultinom(1,tot_kids,fecun))

  sdf = split(sdf,sdf$sex)

  ##sample the labels to generate parents
  kidlabs = lapply(sdf,function(x){
    sample(rep(x$ID,x$nkids))
  })
  kids = do.call(cbind,kidlabs)
  kids


}

#' Breed one year with 2 sexes, overdispersion and sex ratios
  #'
  #' @param population the population to breed
  #' @param year_now the current year
  #' @param overdispersion the vector with the overdispersion parameter for each sex
  #' @param sex_fraction fraction of population that is male
  #'
  #' @export
  breed_one_year_od <- function(population,year_now,overdispersion,sex_fraction){
    setkey(population$population,alive,mature,cur_year)
    population$population[,age := cur_year-birth_year]

    mature = population$population[.(TRUE,TRUE,year_now),]
    ages = sort(unique(mature$age))
    ave_fecun_f = population$fecund_fun(ages,rep("female",length(ages)))
    ave_fecun_m = population$fecund_fun(ages,rep("male",length(ages)))
    ave_fecun = cbind(ave_fecun_f,ave_fecun_m)
    kids = as.data.frame(generate_kids_od(mature$ID,mature$age,mature$Sex,ave_fecun,overdispersion))

    new_kids = data.table(Sex = sample(c("male","female"),nrow(kids),TRUE,c(sex_fraction,1-sex_fraction)),
			    birth_year = year_now,
			    mother=kids$female,
			    father=kids$male,
			    cur_year = year_now+1,
			    alive=TRUE,
			    ID = (max(population$population$ID)+1):(max(population$population$ID)+nrow(kids)),
			    mature = TRUE)
    ##   namsmt = c(population$namsmt,"ID")
    ##   mom_marks = mature[ID %in% new_kids$mother,..namsmt]
    ##   setnames(mom_marks,"ID","mother")
    ## new_kids = merge(new_kids,mom_marks,all.x=TRUE)

    setkey(population$population,ID)
    mdna = as.data.frame(population$population[.(new_kids$mother),])
    fdna = as.data.frame(population$population[.(new_kids$father),])


    nmarks = length(population$smtlist)
msamp = replicate(nrow(new_kids),paste0("m_",1:nmarks,"_",sample(c("A","B"),nmarks,TRUE)))
fsamp = replicate(nrow(new_kids),paste0("m_",1:nmarks,"_",sample(c("A","B"),nmarks,TRUE)))
    msampcb = lapply(1:nrow(new_kids),function(x){
      cbind(x,match(msamp[,x],names(mdna)))})
    msampcb = do.call(rbind,msampcb)
  fsampcb = lapply(1:nrow(new_kids),function(x){
      cbind(x,match(fsamp[,x],names(fdna)))})
    fsampcb = do.call(rbind,fsampcb)

mA = as.data.frame(matrix(as.numeric(mdna[msampcb]),nrow=nrow(new_kids),byrow=TRUE))
names(mA) = paste0("m_",1:nmarks,"_A")
mB = as.data.frame(matrix(as.numeric(fdna[fsampcb]),nrow=nrow(new_kids),byrow=TRUE))
    names(mB) = paste0("m_",1:nmarks,"_B")
    new_kids = cbind(new_kids,mA,mB)

  ##kids_markers = samp_markers(population$smtlist,nrow(new_kids))


  ##     marks = c("ID",population$namsmt)
  ##     setkey(population$population,ID)
  ##     DNAs = unique(population$population[.(ID = new_kids$father),..marks])
  ##     allDNAs = as.data.frame(merge(new_kids,DNAs,by.x="father",by.y="ID",suffixes = c("",".m")))
  ##     allDNAs2 = allDNAs


  ## gloop = system.time({    picksL = as.data.frame(as.table(kids_markers))
  ##     picksL$Var1 = as.numeric(picksL$Var1)
  ##     picksL$Var2 = as.numeric(picksL$Var2)+min(new_kids$ID)-1
  ##     names(picksL) = c("whocares","ID","marker")

  ##     new_kidsL = melt(new_kids,measure.vars=population$namsmt,value.name="geno",variable.name="marker")
  ##     new_kidsL$marker = as.character(new_kidsL$marker)
  ##     DNAsL = melt(DNAs,measure.vars=population$namsmt,value.name="geno",variable.name="marker")
  ##     names(DNAsL)[1] = "father"
  ##     setkey(new_kidsL,ID,marker)
  ##     grab = new_kidsL[picksL[,2:3],.(father,marker)]
  ##     setkey(DNAsL,father,marker)
  ##     new_kidsL[picksL[,2:3],"geno"] = DNAsL[grab,3]
  ##     DNAsL[grab]
  ##     new_kids = dcast(new_kidsL,Sex+birth_year+mother+father+cur_year+alive+ID+mature ~ marker,
  ##       	       value.var = "geno")
  ## })

      ##death bit
      setkey(population$population,cur_year,alive)

      new_year = data.table::copy(population$population[.(year_now,TRUE)])

      setkey(new_year,alive,cur_year)
      new_year[,age := cur_year - birth_year][.(TRUE,year_now),u := runif(.N)][.(TRUE,year_now),alive := (u <= population$surv_fun(age,Sex))]
      setkey(new_year,alive,mature,cur_year)
      new_year[.(TRUE,FALSE,year_now),w := runif(.N)][.(TRUE,FALSE,year_now),mature := (w <= population$maturity_fun(age))] 
      new_year[,age := NULL][,u := NULL][,w := NULL]
	new_year[,cur_year := year_now+1]

      population$population[,age :=NULL]

      population$population = rbind(population$population,new_year)
      data.table::setcolorder(new_kids,names(population$population))

      population$population = rbind(population$population,new_kids)
      population


  }

#' Sample yearly non-fatally
#'
#' This samples yearly non-fatally. It's meant to be run AFTER
#' the simulation is complete as it randomly samples based on alive
#' individuals in each of the years.
#'
#' @param population the population to sample
#' @param sampquota the number at each age class to sample
#' @param samp_years the number of years to sample
#' @export
samp_non_fatally <- function(sim,sampquota,samp_years){
    sim$population[, age := cur_year - birth_year]
    setkey(sim$population,cur_year,alive)
    quts = expand.grid(sampquota,samp_years)

    selected = sim$population[.(samp_years,TRUE),]
    samps = selected[,.SD[sample(.N,sampquota[.BY$age+1])],by=.(age,cur_year)]
    
    sim$population[,age := NULL]
    sim$samps = samps
    sim
}

#' Sample yearly non-fatally with a fraction
#'
#' This samples yearly non-fatally. It's meant to be run AFTER
#' the simulation is complete as it randomly samples based on alive
#' individuals in each of the years. This version uses a
#' sampling fraction instead of a quota.
#'
#' @param population the population to sample
#' @param sampfrac the fraction at each age class to sample
#' @param samp_years the number of years to sample
#' 
#' @export
samp_non_fatally_frac <- function(sim,sampfrac,samp_years){
    sim$population[, age := cur_year - birth_year]
    setkey(sim$population,cur_year,alive)

    selected = sim$population[.(samp_years,TRUE),]
    samps = selected[,.SD[sample(.N,.N*sampfrac[.BY$age+1])],by=.(age,cur_year)]
    
    sim$population[,age := NULL]
    sim$samps = samps
    sim
}

#' make the pedigree of the indivduals
#'
#' @param indiv sim to get pedigree of
#' @export
make_ped <- function(indiv){
    indiv$pedigree = unique(indiv$population,by=c("ID","mother","father"))
    setorder(indiv$pedigree,ID)
    indiv
}

#' Sample markers from the population for breeding
#'
#' @param smtlist the list of markers used in the simulation
#' @param n the number to get back
#'
#' @export
samp_markers <- function(smtlist,n){
    inds = replicate(length(smtlist),sample(1:2,n,replace=TRUE))
    jammed = do.call(rbind,smtlist)
    nams = apply(inds,1,function(x){jammed[cbind(seq_along(x),x)]})
    nams
}

#' Generate rough sample size to target
#' @param tN matrix of true numbers at age
#' @param fun the function to use with 
#' @export
rough_sample_size <- function(tN,fun=function(y){10*sqrt(y)}){
    samp = apply(tN,2,function(x){
        s_est = fun(x)
        ssize = rpois(1,s_est)
        mprop = x/sum(x)
        samp = ssize*mprop
        floor(samp)
    })
    samp
}

#' Sample by age by year non-lethally
#'
#' @param population data.table of population from simulation
#' @param samp_size the matrix of the sample size each year at age
#' @param sample_years vector of years from the sim to sample from39
#' @param max_age the max age of individuals in the population
#' @export
sample_age_by_year <- function(population,samp_size,sample_years,max_age){
    samps = list()
    for(y in sample_years){
        y2 = y-min(sample_years)+1
        setkey(population,alive,cur_year)
        allive = population[.(TRUE,y)]
        allive[, `:=`(age,cur_year-birth_year)]
        ages = list()
        setkey(allive,age)
        for(a in 1:max_age){
            ages[[a]] = allive[.(a)]
        }
        sages = list()
        for(a in 1:max_age){
            sages[[a]] = sample(1:nrow(ages[[a]]),samp_size[a,y2])
        }
        ret = list()
        for(a in 1:max_age){
            ret[[a]] = ages[[a]][sages[[a]],]
        }
        samps[[y2]] = ret
    }
    samps
}
