library(tidyverse)
library(stringr)
library(varCKMR)

sim_df = read.csv("../sims/sim_df.csv")
sim_loc = "../sims/"
files = list.files("../populations/","*.rds")
gfiles = gsub(".rds","",files)
for(i in 1:length(gfiles)){
dir.create(paste0(sim_loc,gfiles)[i])
}

library(reshape)
reps = data.frame(rep=1:1000)

popPar = list()
for(i in 1:length(files)){
  pop = readRDS(paste0("../populations/",files[i]))
  popPar[[i]] = pop$popPar
  popPar[[i]]$foldname = paste0("../sims/",gfiles[i])
}
popPar = do.call(rbind,popPar)
popPar$X = 1:nrow(popPar)

rep_df = expand.grid.df(popPar,reps)    
rep_df = split(rep_df,rep_df$X)
set.seed(46)
rep_df = lapply(rep_df,function(x){
  x$seed = sample(1:1000000,nrow(x))
  x})

nyears = 50
sampyears = 10
sim1y = nyears-sampyears

rep_df1 = rep_df[1:4]
rep_df2 = rep_df[5:8]
rep_df3 = rep_df[9:12]
rep_df4 = rep_df[13:16]
rep_df5 = rep_df[17:21]
rep_df6 = rep_df[22:24]
saveRDS(rep_df1,file="../sims/rep_df1.rds")
saveRDS(rep_df2,file="../sims/rep_df2.rds")
saveRDS(rep_df3,file="../sims/rep_df3.rds")
saveRDS(rep_df4,file="../sims/rep_df4.rds")
saveRDS(rep_df5,file="../sims/rep_df5.rds")
saveRDS(rep_df6,file="../sims/rep_df6.rds")
