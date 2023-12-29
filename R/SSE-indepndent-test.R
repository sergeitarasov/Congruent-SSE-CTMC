##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##----------------------- SEPARABLE BUT DEPENDENT MODEL-------------------------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
##                                                                            --
##-------------------------------- INFERENCES-----------------------------------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                     CID-4 (i.e. 4 states in total) (as HiClaSSE-4)       ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Args <- list(
  Nstates = 4L,
  y = list(
    c(0,0,0,0, 1,0,1,0),
    c(0,0,0,0, 0,1,0,1)
  ))



newArgs <- makeArgs(Args)
printArgsGlobal()
args <- argnames_HiClaSSE(4)
args
args$pars
length(args$pars)

pars.hc <- rep(0, length(args$pars))
names(pars.hc) <- args$pars
pars.hc['lam000'] <- 0.1
pars.hc['lam111'] <- 0.1
pars.hc['lam222'] <- 0.05
pars.hc['lam333'] <- 0.05
pars.hc['mu0'] <-pars.hc['mu1'] <- pars.hc['mu2'] <- pars.hc['mu3'] <-0
qs <- extract_off_diagonal(Q_cid)
qsl=length(qs)
pars.hc[c(length(pars.hc)-qsl+1):length(pars.hc)] <- qs

pars.hc
pars_to_arrays(pars.hc,4)
args

zero.constr <- formulas_zero_pars(pars.hc)
f.qs <- assign_classes_pairwise(Q_cid, args$arrays$Q)
f.lams <- c(lam111 ~ lam000, lam333 ~ lam222)
f.list<- c(zero.constr,  f.qs, f.lams)

#------------------
cid4 <- list()
i=1
for (i in 1:NSIM){

  print(paste0('Working on: ', i))

  tree <- phy[[i]]
  states<- tree$tip.state
  states<- mapvalues(states, from = c("0", "1", "2", "3"), to=c(0, 1, 0, 1) )
  root <- rep(1/4, 4)

  lik.c <- make.HiClasse_cpp(tree, states,  sampling.f=NULL,  strict=T, control=list(backend = "gslode"), newArgs)
  #res <- lik.c(pars.hc, intermediates=T, root=ROOT.GIVEN, root.p=root, condition.surv=T)
  lik.const <- constrain(lik.c, formulae=f.list)
  arg.const <- argnames(lik.const)
  starting.point <- pars.hc[arg.const]
  cid4[[i]] <- find.mle(lik.const, starting.point, method="subplex", keep.func=F, root=ROOT.GIVEN, root.p=root, condition.surv=TRUE)

}

# Results
get_item(cid4, 'lnLik')
pars_to_arrays(cid4[[2]]$par.full, 4)

