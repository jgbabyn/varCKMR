##See README
##here is a small example using the sim parameters from one of the populations

surv = c(0.215,0.280)
death_probs = c(1-surv[1],surv[1]*(1-surv[2]),surv[1]*surv[2])
fecundity = 2*c(0.380,1.518,3.416)
theta = 0.1

EXD = cumsum(fecundity)
V1 = fecundity+fecundity^2/theta
varXD = cumsum(V1)


oldform <- function(varXD,EXD,Delta){
  firstsum <- sum(varXD*Delta)
  secondsum <- sum(EXD^2*(1-Delta)*(Delta))

  thirdsum <- 0
  for(i in 2:length(EXD)){
    for(j in 1:(i-1)){
      thirdsum = thirdsum + EXD[i]*Delta[i]*EXD[j]*Delta[j]
    }
  }

  total <- firstsum + secondsum -2*thirdsum
  total


}

newform <- function(varXD,EXD,Delta){
  mu_tot = sum(EXD*Delta)
  firstsum = sum(varXD*Delta)
  secondsum = sum(EXD^2*Delta)
  total = firstsum+secondsum-mu_tot^2
  total
}

oldv = oldform(varXD,EXD,death_probs)
newv = newform(varXD,EXD,death_probs)
print(paste0("old: ",round(oldv,3)," new: ", round(newv,3)))
