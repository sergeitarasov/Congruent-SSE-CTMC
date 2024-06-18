# See
# James Ledoux. A Necessary Condition
# for Weak Lumpability in Finite Markov Processes


library("corHMM")
library("phytools")
library(plyr)
library(dplyr)
source('R/mod_ray_disc.R')

# Weak lumpability


# make rate matrix Q
Q <- matrix(
  c(
    -18, 6, 12,
    0, -20, 20,
    21, 3, -24
  ), 3,3, byrow = T)

# pi should of the form (1-3a, a, 2a)
pi <- c(1,0,0)
# equibrium pi is
pi.eq = c(7/16, 3/16, 3/8)

# lumping as Q.lum=MAB
M <- matrix(
  c(
    1,0,0,
    0, 1/3, 2/3
  ), 2,3, byrow = T)

B <- matrix(
  c(
    1,0,
    0,1,
    0,1
  ), 3,2, byrow = T)

Q.lump <- M%*%Q%*%B


library(expm)
pi%*%expm(Q*.1)
c(1,0)%*%expm(Q.lump*.1)

c(.7,.1,.2)%*%expm(Q*.1)
c(.7,.3)%*%expm(Q.lump*.1)

c(1/3,1/3,1/3)%*%expm(Q*.1)
c(1/3,2/3)%*%expm(Q.lump*.1)

tree<-pbtree(n=200, scale=100, b=1, d=0)
Qsim <- Q.lump*.1
hist <- sim.history(tree, Qsim, nsim=10)

# Inference under original model
Q2.inf <- Q2model(Qsim, diag.as = NA)
out.org <- list()
for (i in 1:10){
  taxa <- cbind(hist[[i]]$tip.label, hist[[i]]$states)
  out.org[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=Q2.inf, root.p=c(.7, 0.3), p=c(1.4, 1.8))
}


# Inference under EHE transfrom model
Qext <- Q*.1
Qext.inf <- matrix(
  c(
    NA, 1, 2,
    NA, NA, 3,
    4, 5, NA
  ), 3,3, byrow = T)
Qext.inf
Qext
rates <- c(.6, 1.2, 2, 2.1, 0.3)
rates
out.ext <- list()
for (i in 1:10){
  states.ehe<- mapvalues(hist[[i]]$states, from = c("1", "2"), to=c("1", "2&3") )
  taxa<-cbind( names(states.ehe), states.ehe)
  out.ext[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=Qext.inf, root.p=c(.7,0.1,0.2), p=rates)
}

lapply(out.org, function(x) x$loglik) %>% unlist()
lapply(out.ext, function(x) x$loglik) %>% unlist()

out.org[1]
out.ext[1]

