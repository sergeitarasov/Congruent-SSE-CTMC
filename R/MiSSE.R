## GeoSSE equivalence
## Same tree simulated in ?make.geosse
pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
names(pars) <- diversitree:::default.argnames.geosse()
set.seed(5)
phy <- tree.geosse(pars, max.t=4, x0=0)

lik.g <- make.geosse(phy, phy$tip.state)
pars.g <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
names(pars.g) <- argnames(lik.g)

make.classe(phy, phy$tip.state+1, 3)

lik.c <- make.classe(phy, phy$tip.state+1, 3)
pars.c <- 0 * starting.point.classe(phy, 3)
pars.c['lambda222'] <- pars.c['lambda112'] <- pars.g['sA']
pars.c['lambda333'] <- pars.c['lambda113'] <- pars.g['sB']
pars.c['lambda123'] <-  pars.g['sAB']
pars.c['mu2'] <- pars.c['q13'] <- pars.g['xA']
pars.c['mu3'] <- pars.c['q12'] <- pars.g['xB']
pars.c['q21'] <- pars.g['dA']
pars.c['q31'] <- pars.g['dB']

lik.g(pars.g)   # -175.7685
lik.c(pars.c)   # -175.7685


#-----
Q_t <- initQ(c(0,1), c(1,2)) # trait
Q_r <- initQ(c('A','B'), c(3,4)) # diversification regimes
Q_j <- amaSMM(Q_r, Q_t, diag.as = NA)
v=c(1,3,2,4)
Q_j[v,v]

lik.t <- make.bisse.td(phy, phy$tip.state, 2)

#-----------
bookdown::render_book("Rmd", output_format = "bookdown::gitbook", output_dir = '../Report')
bookdown::render_book("Rmd", output_format = "bookdown::pdf_book", output_dir = '../Report-pdf')

Qnew <- edit(Q_cid8.t)
mean <- edit(mean, editor = "xedit")

# Create the matrix as a data frame
matrix_data <- data.frame(
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
rownames(matrix_data) <- colnames(matrix_data)
matrix_data <-as.matrix(matrix_data)


# Create a sample matrix
my_matrix <- matrix(1:12, nrow = 4, ncol = 3)
colnames(my_matrix) <- c("A", "B", "C")
rownames(my_matrix) <- c("Row1", "Row2", "Row3", "Row4")

edited_matrix <- edit(my_matrix)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ installations and dependencies  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# source dependencies and install them if they are not
source('R/utils/dependencies.R')
#source('R-models/my-ClaSSE.R') # pure R implementation of HiCLaSSE
source('R/hiclasse/HiClaSSE-R.R') # pure R implementation of HiCLaSSE
source('R/hiclasse/HiClaSSE_cpp.R') # fast implementation



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##------------------------- SIMULATE DATA: CID-4 MODEL--------------------------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NSIM=10 # N of trees to simulate
max.taxa = 100

Q_hisse <- TransMatMakerHiSSE(hidden.traits=1)
Q_t <- initQ(c(0,1), c(0.05,0.2)) # trait
Q_r <- initQ(c('A','B'), c(0.2,0.05)) # diversification regimes
Q_cid <- amaSMM(Q_r, Q_t, diag.as = NA)
Q_cid
# > Q_cid # same state order as in Q_hisse
#      A0   A1  B0   B1
# A0   NA 0.05 0.2 0.00
# A1 0.20   NA 0.0 0.20
# B0 0.05 0.00  NA 0.05
# B1 0.00 0.05 0.2   NA
v = c(1,3, 2, 4)
Q_cid[v,v] # lumpable with respect to trait

# diversification pars
#lam <- c(0.05, 0.05, 0.1, 0.1)
lam <- c(0.1, 0.1, 0.05, 0.05)
mu <- rep(0, 4)
div.pars <- convert2ratesHisse(lam, mu)
div.pars

phy<- list()
i=1
for (i in 1:NSIM){
  simulated.result <- SimulateHisse(div.pars$turnover, div.pars$eps, Q_cid, max.taxa=max.taxa, x0 = sample(0:3, 1), nstart = 2)
  simulated.result
  phy[[i]] <- SimToPhylo(simulated.result$results, include.extinct=FALSE)
  
}

#saveRDS(phy, file='R/data/phy_CID4.RDS')
phy <- readRDS(file='R/data/phy_CID4.RDS')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      SET-UP INFERENCE                                    ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Q3 <- initQ(c(0,1,2), c(0.01,0.01, 0.01))
Q3 <- make_Qcol(4, c(0.01, 0.03, 0.05, 0))
Q3[1,3] <- Q3[2,1] <- Q3[3,2] <- 0
Q3[4,] <- Q3[4,] <- rep(0, 4)
v=c(1,3,2)
Q3[v,v]
diag(Q3) <- NA
Q3
# > Q3
#       [,1] [,2] [,3] [,4]
# [1,]   NA 0.01 0.00    0
# [2,] 0.00   NA 0.05    0
# [3,] 0.03 0.00   NA    0
# [4,] 0.00 0.00 0.00   NA


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      MiSSE-2                                            ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Args <- list(
  Nstates = 2L,
  y = list(
    c(0,0, 1,1),
    c(0,0, 1,1)
  ))

newArgs <- makeArgs(Args)
printArgsGlobal()
args <- argnames_HiClaSSE(2)
args
args$pars
length(args$pars)

Q3

pars.hc <- rep(0, length(args$pars))
names(pars.hc) <- args$pars
pars.hc['lam000'] <- 0.05
pars.hc['lam111'] <- 0.02
pars.hc['mu0'] <-pars.hc['mu1'] <-0.1
pars.hc['q01'] <- Q3[1,2]
pars.hc['q10'] <- Q3[3,1]

pars.hc
pars_to_arrays(pars.hc, 2)


Misse2 <- list()
i=1
for (i in 1:NSIM){
  print(paste0('Working on: ', i))
  tree <- phy[[i]]
  states<- tree$tip.state
  # states<- mapvalues(states, from = c("0", "1", "2"), to=c(0, 1, 0) )
  states<- mapvalues(states, from = c("0", "1", "2", "3"), to=c(0, 1, 0, 1) )
  root <- rep(1/2, 2)
  
  lik.c <- make.HiClasse_cpp(tree, states,  sampling.f=NULL,  strict=TRUE, control=list(backend = "gslode"), newArgs)
  #res <- lik.c(pars.hc, intermediates=F, root=ROOT.GIVEN, root.p=root, condition.surv=T)
  lik.const <- constrain(lik.c, lam001~0, lam011~0, lam100~0, lam101~0, mu0 ~ mu1)
  arg.const <- argnames(lik.const)
  starting.point <- pars.hc[arg.const]
  Misse2[[i]] <- find.mle(lik.const, starting.point, method="subplex", keep.func=F, root=ROOT.GIVEN, root.p=root, condition.surv=TRUE)
}

get_item(Misse2, 'lnLik')
pars_to_arrays(Misse2[[1]]$par.full, 2)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      MiSSE-3     Congruent                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Args <- list(
  Nstates = 3L,
  y = list(
    c(0,0,0, 1,1,1),
    c(0,0,0, 1,1,1)
  ))

newArgs <- makeArgs(Args)
printArgsGlobal()
args <- argnames_HiClaSSE(3)
args
args$pars
length(args$pars)

Q30 <- Q3
Q30 <- Q30[-4,-4]
Q30[2,3] <- 0.01
Q30[2,1] <- 0.01
diag(Q30) <- 0
Q30

pars.hc <- rep(0, length(args$pars))
names(pars.hc) <- args$pars
pars.hc['lam000'] <- 0.05
pars.hc['lam111'] <- 0.02
pars.hc['lam222'] <- 0.02
pars.hc['mu0'] <-pars.hc['mu1'] <-pars.hc['mu2'] <-0.1
qs <- extract_off_diagonal(Q30)
pars.hc[22:27] <- qs

pars.hc
pars_to_arrays(pars.hc, 3)

zero.constr <- formulas_zero_pars(pars.hc)
#f.qs <- c(q10 ~ q01, q02 ~ q01, q03 ~ q01, q20 ~ q01, q23 ~ q01, q31 ~ q01, q32 ~ q01, q12 ~ q13)
#f.qs <- c(q02 ~ q03, q12 ~ q13, q20 ~ q13, q21 ~ q13, q30 ~ q03, q31 ~ q03)
f.qs <- c(q12 ~ q20, q10 ~ q20)
f.mu <- c(mu1 ~ mu0, mu2 ~ mu0)
f.lams <- c(lam222 ~  lam111)
f.list<- c(zero.constr,  f.lams, f.mu, f.qs)



Misse3C <- list()
i=1
for (i in 1:NSIM){
  print(paste0('Working on: ', i))
  tree <- phy[[i]]
  states<- tree$tip.state
  states<- mapvalues(states, from = c("0", "1", "2", "3"), to=c(0, 1, 0, 1) )
  root <- c(1/2, 1/4, 1/4)
  
  lik.c <- make.HiClasse_cpp(tree, states,  sampling.f=NULL,  strict=TRUE, control=list(backend = "gslode"), newArgs)
  #res <- lik.c(pars.hc, intermediates=F, root=ROOT.GIVEN, root.p=root, condition.surv=T)
  lik.const <- constrain(lik.c, formulae = f.list)
  arg.const <- argnames(lik.const)
  starting.point <- pars.hc[arg.const]
  Misse3C[[i]] <- find.mle(lik.const, starting.point, method="subplex", keep.func=F, root=ROOT.GIVEN, root.p=root, condition.surv=TRUE)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      MiSSE-3     Inferred Pars                          ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pars_to_arrays(Misse2[[1]]$par.full, 2)
pars_to_arrays(Misse3C[[1]]$par.full, 3)







Misse3C_pars <- list()
i=1
for (i in 1:NSIM){
  print(paste0('Working on: ', i))
  #---
  mi2_est <- get_item(Misse2[i], 'par')
  pars.hc <- rep(0, length(args$pars))
  names(pars.hc) <- args$pars
  pars.hc['lam000'] <- mi2_est['lam000']
  pars.hc['lam111'] <- mi2_est['lam111']
  pars.hc['lam222'] <- mi2_est['lam111']
  pars.hc['mu0'] <-pars.hc['mu1'] <-pars.hc['mu2'] <- mi2_est['mu1']
  qs <- extract_off_diagonal(Q30)
  pars.hc[22:27] <- c(mi2_est['q01'], 0, mi2_est['q10'], mi2_est['q10'], mi2_est['q10'], 0)
  #pars.hc
  #pars_to_arrays(pars.hc, 3)
  #---
  tree <- phy[[i]]
  states<- tree$tip.state
  states<- mapvalues(states, from = c("0", "1", "2", "3"), to=c(0, 1, 0, 1) )
  root <- c(1/2, 1/4, 1/4)
  
  lik.c <- make.HiClasse_cpp(tree, states,  sampling.f=NULL,  strict=TRUE, control=list(backend = "gslode"), newArgs)
  Misse3C_pars[[i]] <- lik.c(pars.hc, intermediates=F, root=ROOT.GIVEN, root.p=root, condition.surv=T)
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      MiSSE-3     Non-congruent                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Args <- list(
  Nstates = 3L,
  y = list(
    c(0,0,0, 1,1,1),
    c(0,0,0, 1,1,1)
  ))

newArgs <- makeArgs(Args)
printArgsGlobal()
args <- argnames_HiClaSSE(3)
args
args$pars
length(args$pars)

Q30 <- Q3
Q30 <- Q30[-4,-4]
Q30[2,3] <- 0.01
diag(Q30) <- 0
Q30

pars.hc <- rep(0, length(args$pars))
names(pars.hc) <- args$pars
pars.hc['lam000'] <- 0.05
pars.hc['lam111'] <- 0.02
pars.hc['lam222'] <- 0.02
pars.hc['mu0'] <-pars.hc['mu1'] <-pars.hc['mu2'] <-0.1
qs <- extract_off_diagonal(Q30)
pars.hc[22:27] <- qs

pars.hc
pars_to_arrays(pars.hc, 3)

zero.constr <- formulas_zero_pars(pars.hc)
#f.qs <- c(q10 ~ q01, q02 ~ q01, q03 ~ q01, q20 ~ q01, q23 ~ q01, q31 ~ q01, q32 ~ q01, q12 ~ q13)
#f.qs <- c(q02 ~ q03, q12 ~ q13, q20 ~ q13, q21 ~ q13, q30 ~ q03, q31 ~ q03)
f.qs <- c(q12 ~ q20)
f.mu <- c(mu1 ~ mu0, mu2 ~ mu0)
f.lams <- c(lam222 ~  lam111)
f.list<- c(zero.constr,  f.lams, f.mu, f.qs)



Misse3 <- list()
i=1
for (i in 1:NSIM){
  print(paste0('Working on: ', i))
  tree <- phy[[i]]
  states<- tree$tip.state
  states<- mapvalues(states, from = c("0", "1", "2", "3"), to=c(0, 1, 0, 1) )
  root <- c(1/2, 1/4, 1/4)
  
  lik.c <- make.HiClasse_cpp(tree, states,  sampling.f=NULL,  strict=TRUE, control=list(backend = "gslode"), newArgs)
  #res <- lik.c(pars.hc, intermediates=F, root=ROOT.GIVEN, root.p=root, condition.surv=T)
  lik.const <- constrain(lik.c, formulae = f.list)
  arg.const <- argnames(lik.const)
  starting.point <- pars.hc[arg.const]
  Misse3[[i]] <- find.mle(lik.const, starting.point, method="subplex", keep.func=F, root=ROOT.GIVEN, root.p=root, condition.surv=TRUE)
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      Results                                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mi2 = get_item(Misse2, 'lnLik')
mi3c = get_item(Misse3C, 'lnLik')
mi3c_pars = unlist(Misse3C_pars)
mi3 = get_item(Misse3, 'lnLik')

mi2
mi3c
mi3c_pars
mi3

abs(mi2-mi3c) %>% mean()
abs(mi2-mi3c_pars) %>% mean()
abs(mi2-mi3) %>% mean()

get_item(Misse2[1], 'par')
get_item(Misse3C[1], 'par')
get_item(Misse3[1], 'par')


pars_to_arrays(Misse2[[1]]$par.full, 2)
pars_to_arrays(Misse3C[[1]]$par.full, 3)
pars_to_arrays(Misse3[[1]]$par.full, 3)


