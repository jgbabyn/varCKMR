#' Make Leslie Matrix
#'
#' Creates a simple Leslie matrix from fecundity and survival vectors
#' @param fecun the fecundity vector
#' @param surv the survival vector
#'
#' @export
make_les <- function(fecun,surv){
    Les = matrix(0,nrow=length(fecun),ncol=length(fecun))
    Les[1,] = fecun
    for(i in 1:length(surv)){
        Les[i+1,i] = surv[i]
    }
    Les
}

#' Find deterministic growth rate of Leslie Matrix
#'
#' @param Les the Leslie matrix to find the growth rate of
#'
#' @export
get_growth <- function(Les){
    ev <- eigen(Les,only.values=TRUE)$values
    Re(ev[1])
}

#' Simple fecundity at age function
#'
#' @param age the age of the individual
#' @param fpar the fecundity parameter
#' @export
fecund_inv_fun <- function(age,fpar){
    exp(fpar)*age^2
}

#' Keeps survival within a specified range
#'
#' @param y a value
#' @param a the lower part of the range
#' @param b the upper part of the range
#' @export
surv_inv_fun <- function(y,a,b){
    x = a+(b-a)*plogis(y)
    x
}

#' Find missing parameter for target growth rate
#'
#' @param target the target deterministic growth rate to hit
#' @export
get_pars_with_target <- function(target){
    pars = c(runif(1,-2,2))
    fpar = runif(1,-1,-.5)
    y = runif(1,0.1,0.4)

    optomax = function(pars){
        fecun = fecund_inv_fun(1:3,fpar)
        survt = surv_inv_fun(pars[1],0.1,0.9)
        surv = c(survt,y)
        Les <- make_les(fecun,surv)
        gw = get_growth(Les)
        abs(target-gw)
    }
    opt = optim(pars,optomax,method="Brent",lower=-5,upper=5)
    print(opt)

    fex = fecund_inv_fun(1:3,fpar)
    ##surv = c(surv_inv_fun(opt$root,0.1,0.7),y)
    surv = c(surv_inv_fun(opt$par,0.1,0.9),y)
    Les = make_les(fex,surv)
    gw = get_growth(Les)

    ret = list(growth_rate=gw,fec=fex,surv=surv,fpar=fpar)
    ret

}

#' Get numbers at age matrix from population
#'
#' @param pop the population to get numbers at age for
#' @export
n_at_age <- function(pop){
    population = pop$population[pop$population$alive == TRUE]
    population$age = population$cur_year - population$birth_year
    table(population$cur_year,population$age)
}
