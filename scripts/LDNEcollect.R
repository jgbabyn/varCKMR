LDNEfolders = list.dirs("../LDNEsamples",recursive = FALSE)
library(tidyverse)

LDNEQ = list()
for(i in 1:length(LDNEfolders)){
    file = list.files(LDNEfolders[i],".rds",full.names=TRUE)
    ldnes = readRDS(file[1])
    popbit =  regexec("POP[0-9][0-9][0-9][0-9]",LDNEfolders[i])
    mlength = sapply(popbit,function(x){attributes(x)$match.length})
    popID = unique(do.call(substring,list(text=LDNEfolders[i],first=popbit,last=unlist(popbit)+mlength-1)))
    popID = tolower(popID)
    
    LDNEQ[[i]] = ldnes |>
        group_by(year) |>
        summarize(Ne= quantile(Neadj3,probs=c(0.05,0.5,0.95)),q=c(0.05,0.5,0.95))
    LDNEQ[[i]]$population = popID
}

LDNEQC = do.call(rbind,LDNEQ)
saveRDS(LDNEQC,"../LDNEQuantilesC.rds")
