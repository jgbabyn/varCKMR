library(TMB)
library(tidyverse)
library(stringr)
library(varCKMR)
library(abind)
library(adegenet)
library(strataG)

pop_files = c("../populations/pop0022.rds","../populations/pop0023.rds","../populations/pop0024.rds")



for(i in 1:length(pop_files)){
    

    if(!dir.exists(paste0("../LDNEsamples/","LDNEPOP",str_pad(i+21,4,pad="0")))){
        dir.create(paste0("../LDNEsamples/","LDNEPOP",str_pad(i+21,4,pad="0")))
    }
    savepath = paste0("../LDNEsamples/","LDNEPOP",str_pad(i+21,4,pad="0"))

    indiv = readRDS(pop_files[[i]])
    ped = make_ped(indiv)
    tN = n_at_age(indiv)[,]

    popPar = indiv$popPars  

    dprob = c(1-popPar$surv1,popPar$surv1*(1-popPar$surv2),popPar$surv2*popPar$surv1)
    fec = c(popPar$fec1,popPar$fec2,popPar$fec3)
    LDNEs = list()

    for(j in 1:1000){
        print(paste("i: ",i,"j: ",j)) 
        sampsForLDNe = varCKMR::rough_sample_size(t(tN),function(y){0.15*y})
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
        ldNeEst$Nbadj2 = ldNeEst$rawNb/(0.991-0.206*log10(AL)+0.256*log10(1)+0.137*CVf)
        ldNeEst$Neadj2 = ldNeEst$Nbadj/(0.833+0.637*log10(AL)-0.793*log10(1)-0.423*CVf)
        AL2 = 4-1+1
        ldNeEst$Nbadj3 = ldNeEst$rawNb/(0.991-0.206*log10(AL2)+0.256*log10(1)+0.137*CVf)
        ldNeEst$Neadj3 = ldNeEst$Nbadj/(0.833+0.637*log10(AL2)-0.793*log10(1)-0.423*CVf)
        ldNeEst$sim = j
        ldNeEst$year = 1:11
        LDNEs[[j]] = ldNeEst
    }
LDNEscom = do.call(rbind,LDNEs)
saveRDS(LDNEscom,paste0(savepath,"/LDNEcom",str_pad(i+21,width=4,pad=0),".rds"))
}
