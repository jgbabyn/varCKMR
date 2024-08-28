library(TMB)
library(tidyverse)
library(stringr)
library(varCKMR)
library(abind)
library(adegenet)
library(strataG)

pop_files = c("../populations/pop0020.rds","../populations/pop0002.rds","../populations/pop0003.rds","../populations/pop0004.rds","../populations/pop0005.rds","../populations/pop0006.rds")

for(i in 1:length(pop_files)){
    

    if(!dir.exists(paste0("../LDNEsamples/","LDNEPOP",str_pad(i,4,pad="0")))){
        dir.create(paste0("../LDNEsamples/","LDNEPOP",str_pad(i,4,pad="0")))
    }
    savepath = paste0("../LDNEsamples/","LDNEPOP",str_pad(i,4,pad="0"))

    indiv = readRDS(pop_files[[i]])
    ped = make_ped(indiv)
    tN = n_at_age(indiv)[,]

    popPar = indiv$popPars  

    dprob = c(1-popPar$surv1,popPar$surv1*(1-popPar$surv2),popPar$surv2*popPar$surv1)
    fec = c(popPar$fec1,popPar$fec2,popPar$fec3)
    LDNEs = list()

    for(j in 1:1000){
        print(paste("i: ",i,"j: ",j)) 
        sampsForLDNe = varCKMR::rough_sample_size(t(tN),function(y){0.25*y})
        sampForLDNe = varCKMR::sample_age_by_year(indiv$population,sampsForLDNe,40:50,1)
        samp2ForLDNe = lapply(sampForLDNe,function(x){
	    do.call(rbind,x)
        })

        getLDNE = function(ped){
	    pp = lapply(ped,function(by50){
                Agroup = select(by50,ends_with("_A"))[,1:100]
                Bgroup = select(by50,ends_with("_B"))[,1:100]

                ABpaste = mapply(function(x,y){
                    paste0(x,y)},Agroup,Bgroup)

                genindy = df2genind(ABpaste,ploidy=2,ncode=1)
                gi.g <- genind2gtypes(genindy)
                ldNeey = ldNe(gi.g,maf.threshold = 0.05,num.cores = 8)
                ldNeey
	    })
	    do.call(rbind,pp)
        }

        getLDNE2 = function(N){
	    ldNe = list()
	    for(i in 1:N){
                sampsForLDNe = varCKMR::rough_sample_size(t(tN),function(y){0.15*y})
                sampForLDNe = varCKMR::sample_age_by_year(indiv$population,sampsForLDNe,40:50,1)
                samp2ForLDNe = lapply(sampForLDNe,function(x){
                    do.call(rbind,x)
                })
                ldNeT = getLDNE(samp2ForLDNe)
                ldNeT$rawNb = 1/(3*(ldNeT[,4]-ldNeT[,5]))
                CVf = sd(fec)/mean(fec)
                AL = 3-1+1
                ldNeT$Nbadj = ldNeT$rawNb/(0.991-0.206*log10(AL)+0.256*log10(3)+0.137*CVf)
                ldNeT$Neadj = ldNeT$Nbadj/(0.833+0.637*log10(AL)-0.793*log10(3)-0.423*CVf)
                ldNe[[i]] = ldNeT
	    }
	    ldNe
        }

        
        ldNeEst = getLDNE(samp2ForLDNe)
        ldNeEst$rawNb = 1/(3*(ldNeEst[,4]-ldNeEst[,5]))
        CVf = sd(fec)/mean(fec)
        AL = 3-1+1
        ldNeEst$Nbadj = ldNeEst$rawNb/(0.991-0.206*log10(AL)+0.256*log10(3)+0.137*CVf)
        ldNeEst$Neadj = ldNeEst$Nbadj/(0.833+0.637*log10(AL)-0.793*log10(3)-0.423*CVf)
        ldNeEst$Nbadj2 = ldNeEst$rawNb/(0.991-0.206*log10(AL)+0.256*log10(0.0001)+0.137*CVf)
        ldNeEst$Neadj2 = ldNeEst$Nbadj/(0.833+0.637*log10(AL)-0.793*log10(0.0001)-0.423*CVf)
        AL2 = 4-1+1
        ldNeEst$Nbadj3 = ldNeEst$rawNb/(0.991-0.206*log10(AL2)+0.256*log10(0.001)+0.137*CVf)
        ldNeEst$Neadj3 = ldNeEst$Nbadj/(0.833+0.637*log10(AL2)-0.793*log10(0.001)-0.423*CVf)
        ldNeEst$sim = j
        ldNeEst$year = 1:11
        LDNEs[[j]] = ldNeEst
    }
LDNEscom = do.call(rbind,LDNEs)
saveRDS(LDNEscom,paste0(savepath,"/LDNEcom",str_pad(i,width=4,pad=0),".rds"))
}


les = make_les(fecun=c(popPar$fec1,popPar$fec2,popPar$fec3),surv=c(popPar$surv1,popPar$surv2))
esses = c(popPar$surv1,popPar$surv2)
fecun = c(popPar$fec1,popPar$fec2,popPar$fec3)

get_ell <- function(esses){
    ells = c(1)
    for(i in 2:(length(esses)+1)){
        ells = c(ells,(1-esses[i-1])*ells[i-1])
    }
    ells
}
ells = get_ell(esses)

get_T = function(lambda,ells,bs){
    T = 0
    for(i in 1:length(ells)){
        T = T+i*ells[i]*bs[i]*lambda^(-i)
    }
    T
}

get_vi <- function(lambda,ells,bs){
    vs = c(1)
    for(i in 2:length(ells)){
        vpart = 0
        for(j in i:length(ells)){
            vpart = vpart+ells[j]*bs[j]*lambda^(-j)
        }
        vs = c(vs,(lambda^(i-1)/ells[i])*vpart)
    }
    vs
}

get_K <- function(lambda,sigmas,ells,esses,vs){
    K1 = 0
    K2 = 0
    for(i in 1:length(ells)){
        K1 = K1 + sigmas[i]*ells[i]*lambda^(-i)
        K2 = K2 + ells[i]*lambda^(-i)*esses[i]*(1-esses[i])*vs[i]
    }
    K = K1 + K2
    K
}

Tt = get_T(popPar$theta,ells,2*fecun)
vs = get_vi(popPar$theta,ells,2*fecun)

sigmas = fecun+(2*fecun)/popPar$theta
Kk = get_K(popPar$theta,sigmas,ells,c(1,esses),vs)
(tN[41:51,1]*Tt)/Kk

  get_drift <- function(pedigree,stupidn=1000){
    pedigree = dplyr::mutate(pedigree,across(starts_with("m_"),~ .x-1))
    pop = as.data.frame(dplyr::select(pedigree,-Sex,-mother,-father,-cur_year,-alive,-ID,-mature))
      popn = dplyr::group_by(pop,birth_year) |>
	dplyr::count()
      pop = dplyr::group_by(pop,birth_year) |>
	dplyr::summarise_all(sum)

    popAs = dplyr::select(pop,ends_with("A"))
    popBs = dplyr::select(pop,ends_with("B"))

    pop = popAs+popBs
    pop = cbind(birth_year=popn$birth_year,pop)

    pedigreePar = pedigree |>
      dplyr::select(birth_year,mother,father) |>
      tidyr::gather(key=parent,value=ID,-birth_year) |>
      dplyr::filter(!is.na(ID)) |>
      dplyr::group_by(birth_year) |>
      dplyr::distinct(ID,.keep_all=TRUE)

    pedigreeParN = pedigreePar |>
      dplyr::group_by(birth_year) |>
      dplyr::count()

    popPar = pedigree[pedigreePar$ID,] |>
      dplyr::mutate(birth_year = pedigreePar$birth_year) |>
      dplyr::select(-Sex,-mother,-father,-cur_year,-alive,-ID,-mature) |>
      dplyr::group_by(birth_year) |>
      dplyr::summarise_all(sum)
    popParAs = dplyr::select(popPar,ends_with("A"))
    popParBs = dplyr::select(popPar,ends_with("B"))
    popPars = popParAs+popParBs
    parfreqs = popPars/(pedigreeParN$n*2)
    parfreq = parfreqs[,]

      childfreq = pop[,-1]/(popn$n*2)
      childfreq = cbind(birth_year=pop$birth_year,childfreq) |>
	dplyr::filter(birth_year >= 0) |>
	dplyr::select(-birth_year)

      differences = childfreq-parfreq
      Vdelta = apply(differences,1,var)

      ps = parfreq*(1-parfreq)
	Ne = ps/(2*Vdelta)

	list(Ne=Ne,Vdelta=Vdelta,parfreq=parfreq)

      }

(tN[41:51,1]*Tt)/Kk

1/(warp$Vdelta[41:50]*Tt)

library(tidyverse)
get_drift2 <-function(population){
    pop2 = dplyr::mutate(population,across(starts_with("m_"),~ .x-1))
    borned = unique(pop2,by=c("ID","mother","father"))


    
bornedAs = dplyr::select(borned,ends_with("A"))
    bornedBs = dplyr::select(borned,ends_with("B"))

    borneds = bornedAs+bornedBs
    borneds = cbind(birth_year=borned$birth_year,borneds)
    borneds = mutate(borneds,across(starts_with("m_"),~ .x/2))

    bornedsG = borneds |>
        group_by(birth_year) |>
        summarise_all(mean) |>
        rename(year=birth_year)

    bornedsG

    alives = pop2 |>
        filter(alive=TRUE,mature=TRUE)

    alivesAs = dplyr::select(alives,ends_with("A"))
    alivesBs = dplyr::select(alives,ends_with("B"))

    alivess = alivesAs+alivesBs
    alivess = cbind(cur_year=alives$cur_year,alivess)
    alivess = cbind(birth_year=alives$birth_year,alivess)
    alivess = mutate(alivess,across(starts_with("m_"),~ .x/2))

    alivessG = alivess |>
        group_by(cur_year) |>
        filter(cur_year != birth_year) |>
        summarise(across(starts_with("m_"),~ mean(.x))) |>
        rename(year=cur_year)

    ret = list()
    ret$alives = filter(alivessG,year %in% bornedsG$year)
    ret$borneds = filter(bornedsG,year %in% alivessG$year)
    ret
}

varp = (1/(apply(dorp$alives-dorp$borneds,2,var)))
(tN[41:51,1]*Tt)/varp[41:50]
