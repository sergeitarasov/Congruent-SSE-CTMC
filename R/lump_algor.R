setwd("~/Library/CloudStorage/OneDrive-UniversityofHelsinki/My_papers/State_Curvature/Github/Congruent-SSE-CTMC")

library(DT)
library(ggplot2)
library(dplyr)
library(knitr)
library(kableExtra)
library(here)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ installations and dependencies  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source('R/utils/dependencies.R')
source('R/hiclasse/HiClaSSE-R.R') # pure R implementation of HiCLaSSE
source('R/hiclasse/HiClaSSE_cpp.R') # fast implementation



# Create the matrix as a data frame
Q_ehe8_C.t <- data.frame(
  a0 = c(-4, 0, 1, 1, 1, 1, 0, 0),
  b0 = c(0, -4, 1, 1, 1, 1, 0, 0),
  c0 = c(1, 1, -4, 0, 0, 0, 1, 1),
  d0 = c(1, 1, 0, -4, 0, 0, 1, 1),
  a1 = c(1, 1, 0, 0, -4, 0, 1, 1),
  b1 = c(1, 1, 0, 0, 0, -4, 1, 1),
  c1 = c(0, 0, 1, 1, 1, 1, -4, 0),
  d1 = c(0, 0, 1, 1, 1, 1, 0, -4)
)

# Set row names and column names
rownames(Q_ehe8_C.t) <- colnames(Q_ehe8_C.t)
Q_ehe8_C.t <-as.matrix(Q_ehe8_C.t)


v=c(1,5, 2,6, 3,7, 4,8)
Q_ehe8_C.r <- Q_ehe8_C.t[v,v]

print(Q_ehe8_C.r)
#print(La8)
print(Q_ehe8_C.t)
#print(La8[v,v])

library(rphenoscate)


# Q <- initQ(c(1, 2), c(.3,.2))
# Q.h <- EHEtransform(Q)
# 
# part_scheme=list(c(1, 2), c(3,4,5))
# is_slumpable(Q.h, part_scheme)

part_scheme=list(c(1, 2), c(3,4), c(5,6), c(7,8))
is_slumpable(Q_ehe8_C.t, part_scheme)

#----- check irreducibility

# 100 tips
Qs <- readRDS( file='R/data/out/SEM/best_Qs-SEM8_phy_CID4-100tips-100tr.RDS')
vector_list <-lapply(Qs, function(x) c(x))
vector_matrix <- do.call(rbind, vector_list)
unique_vectors <- unique(vector_matrix)
nrow(unique_vectors)
# [1] 95

part_scheme=list(c(1, 2), c(3,4), c(5,6), c(7,8))
lapply(Qs, function(x) is_slumpable(x, part_scheme)) %>% unlist %>% any()

#----- check correlation
Q <- initQ(c('a', 'b', 'c', 'd'), c(1:12))
Q2<- initQ(c(0,1), c(13, 14))
smm=amaSMM(Q2, Q)

Qin <- Qs[[1]]

is_trait_INdependent <- function(Qin){
  Q <- initQ(c('a', 'b', 'c', 'd'), c(1:12))
  Q2<- initQ(c(0,1), c(13, 14))
  smm=amaSMM(Q2, Q)
  
  i=3
  # correlation according to Pagel
  out <- c()
  for (i in 1:12){
    cells <- which(smm==i)
    rates=Qin[cells]
    res <- all(rates==rates[1])
    out <- c(out, res)
  }
  
  # dual transitions
  #which((Q_ehe8_C.t + smm)==0)
  duals = c(7,  8, 15, 16, 21, 22, 29, 30, 35, 36, 43, 44, 49, 50, 57, 58)
  drates=Qin[duals]
  res <- all(drates==0)
  out <- c(out, res)
  
  return(all(out))
}


is_trait_INdependent(Q_ehe8_C.t)
lapply(Qs, function(x) is_trait_INdependent(x)) %>% unlist %>% any()
is_indep=lapply(Qs, function(x) is_trait_INdependent(x)) %>% unlist
Qin=Qs[is_indep][[1]]
is_trait_INdependent(Qin)

#_-----------------------------


cid.aic=get_item(readRDS(file='R/data/out/lik-CID4-phy_CID4-100tips-100tr.RDS'), 'lnLik') %>% get_aic(., 4)
sem.aic=get_item(readRDS(file='R/data/out/SEM/best_Ln-SEM8_phy_CID4-100tips-100tr.RDS'), 'lnLik') %>% get_aic(., 2)
daic=cid.aic-sem.aic
hist(daic)

daic[is_indep]
hist(daic, probability =T)
daic[!is_indep] %>% hist(., probability = T)




