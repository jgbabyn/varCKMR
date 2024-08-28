library(tidyverse)
NESCom1 = readRDS("../NESCcom1B.rds")


get_sameys <- function(dat){
    dat = dat |>
    group_by(pop,year) |>
        summarise(simNE=unique(sim_Ne),LD5=unique(LD5),LD50=unique(LD50),LD95=unique(LD95))
    dat
}

get_mod_quan <- function(dat){
    dat = dat |>
        group_by(pop,year) |>
        summarize(mNe = quantile(modelNe,c(0.05,0.5,0.95)),q=c(0.05,0.5,0.95)) |>
        spread(key=q,value=mNe)
    dat
}    


combos = list()
for(i in 1:6){
    NEScom = readRDS(paste0("../NESCcom",i,"B.rds"))
    sames = get_sameys(NEScom)
    mod_quan = get_mod_quan(NEScom)
    combos[[i]] = cbind(sames,mod_quan)
}

CC = do.call(rbind,combos)
names(CC) = c("pop1","year1","simNe","LD5","LD50","LD95","pop2","year2","M5","M50","M95")

CCS = split(CC,CC$pop1)

library(plotly)
CC1 = CCS[[1]]

t2 <- list(size = 20)
m = list(pad = 50)

p1 = plot_ly(CC1,x=~40:50,y=~M50,type="scatter",mode="lines",name=TeX("\\hat{N}_e \\text{CKMR 50th Per.}"),color=I("blue"),line=list(dash="dash",width=4)) |>
    add_ribbons(ymin=~M5,ymax=~M95,name=TeX("\\hat{N}_e \\text{CKMR 5/95th Per.}"),color=I("blue"),line=list(width=2)) |>
    add_trace(x=~40:50,y=~LD50,name=TeX("\\hat{N}_{e(Adj3)} \\text{50th Per.}"),color=I("orange"),line=list(dash="dot",width=4)) |>
    add_ribbons(ymin=~LD5,ymax=~LD95,name=TeX("\\hat{N}_{e(Adj3)} \\text{5/95th Per.}"),color=I("orange"),line=list(dash="dot",width=2)) |>
    add_trace(x=~40:50,y=~simNe,name=TeX("\\text{Eq. 1 } N_{e}"),color=I("black"),line=list(dash="solid",width=4)) |>
    layout(title = list(text=TeX("N_{e} \\text{ for Pop. when }\\theta=0.1\\text{ and }\\lambda=0.95"),font=t2,yref="paper"),yaxis=list(title=TeX("N_{e}"),range=c(0,250)),xaxis=list(title="Year"),legend=list(font=list(size= 20))) |>
    config(mathjax = "cdn")
    
save_image(p1,file="/home/jonathan/Desktop/p1.png")

###Plot the sens analysis:

sensdat = readRDS("../senanalysisresults.rds")

Eests = lapply(1:length(sensdat),function(x){
    MBNe = sensdat[[x]]$N_eEstMB$estN_e
    MVP1Ne = sensdat[[x]]$N_eEstMVP1$estN_e
    MVP2Ne = sensdat[[x]]$N_eEstMVP2$estN_e
    df = data.frame(MBNe=MBNe,MVP1Ne=MVP1Ne,MVP2Ne=MVP2Ne,sim=x,year=40:50)
    df
})

rests = do.call(rbind,Eests)
library(tidyverse)

qrests = rests |>
    group_by(year) |>
    summarize(MBNe=quantile(MBNe,c(0.05,0.5,0.95)),MVP1Ne=quantile(MVP1Ne,c(0.05,0.5,0.95)),MVP2Ne=quantile(MVP2Ne,c(0.05,0.5,0.95)),q=c(0.05,0.5,0.95))

MVP2 = qrests |>
    select(year,MVP2Ne,q) |>
    spread(key=q,value=MVP2Ne)
names(MVP2) = c("year","P25","P250","P295")

MVP1 = qrests |>
    select(year,MVP1Ne,q) |>
    spread(key=q,value=MVP1Ne)
names(MVP1) = c("year","P15","P150","P195")

gorp = cbind(CC1,MVP2,MVP1)

p2 = plot_ly(data=gorp,x=~year,y=~simNe,name=TeX("\\text{Eq. 1 } N_{e}"),color=I("black"),line=list(dash="solid",width=4),mode="lines",type="scatter") |>
    add_trace(x=~year,y=~M50,name=TeX("\\text{Neg. Bin 50th Per.}"),color=I("blue"),line=list(dash="dash",width=4)) |>
    add_trace(x=~year,y=~M5,name=TeX("\\text{Neg. Bin 5/95th Per.}"),legendgroup=TeX("\\text{Neg. Bin 5/95th Per.}"),color=I("blue"),line=list(dash="dash",width=2)) |>
    add_trace(x=~year,y=~M95,name=TeX("\\text{Neg. Bin 5/95th Per.}"),showlegend=FALSE,color=I("blue"),line=list(dash="dash",width=2)) |>                                     
    add_trace(x=~year,y=~P250,name=TeX("\\gamma = 2, \\text{50th Per.}"),line=list(dash="dot",width=4),color=I("orange")) |>
    add_trace(x=~year,y=~P25,name=TeX("\\gamma = 2, \\text{5/95th Per.}"),color=I("orange"),line=list(dash="dot",width=2)) |>
    add_trace(x=~year,y=~P295,name=TeX("\\gamma = 2, \\text{5/95th Per.}"),showlegend=FALSE,color=I("orange"),line=list(dash="dot",width=2)) |>
    add_trace(x=~year,y=~P150,name=TeX("\\gamma = 1, \\text{50th Per.}"),line=list(dash="dashdot",width=4),color=I("red")) |>
    add_trace(x=~year,y=~P15,name=TeX("\\gamma = 1, \\text{5/95th Per.}"),color=I("red"),line=list(dash="dashdot",width=2)) |>
    add_trace(x=~year,y=~P195,name=TeX("\\gamma = 1, \\text{5/95th Per.}"),color=I("red"),line=list(dash="dashdot",width=2)) |>
    layout(xaxis=list(title="Year"),yaxis=list(title=TeX("\\text{Effective Population Size,}, N_{e}"),range=c(0,350)),legend=list(font=list(size= 20))) |>
    config(mathjax="cdn")

