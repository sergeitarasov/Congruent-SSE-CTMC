source('R/hiclasse/HiClaSSE-utils.R')
sourceCpp("R/hiclasse/HiClaSSE_func.cpp")

#  R CMD SHLIB hiclasse_c.c
#dyn.unload("R-models/hiclasse_c.so")
dyn.load("R/hiclasse/hiclasse_c.so")
is.loaded("do_derivs_hiclasse2", "hiclasse_c") %>% print()
is.loaded("derivs_hiclasse2_gslode", "hiclasse_c") %>% print()
is.loaded("derivs_hiclasse3_gslode", "hiclasse_c") %>% print()
is.loaded("derivs_hiclasse8_gslode", "hiclasse_c") %>% print()

derivs.HiClasse_cpp <- function (t, y, pars)
{
  array <- cpp_pars_to_arrays(pars, MY_Nstates)
  lam.tensor <- array$lam.tensor
  mu <- array$mu
  Q <- array$Q

  E <- y[1:MY_Nstates]
  D <- y[-1:-MY_Nstates]
  DE <- D %*% t(E)
  EE <- E %*% t(E)

  # Equation for Lambdas
  #matrix_times_2 <- make_matrix_diag2(MY_Nstates)
  #matrix_times_half <- OffDiagonalHalve(MY_Nstates)
  #flat_vec <- rep(1,MY_Nstates)

  Lambdas <-lapply(lam.tensor, function(x) flat_vec %*% (makeSymmetricMatrix(x)*matrix_times_2*DE) %*% flat_vec )
  Lambdas <- unlist(Lambdas)
  # Equation for Extinctions
  Mus <- lapply(lam.tensor, function(x) flat_vec %*% (makeSymmetricMatrix(x)*matrix_times_half*EE) %*% flat_vec )
  Mus <- unlist(Mus)

  NoChnageLambdas <-  lapply(lam.tensor, function(x) -1*(flat_vec %*% x %*% flat_vec) )
  NoChnageLambdas <- unlist(NoChnageLambdas)

  # Eeq Deq
  c(mu + (NoChnageLambdas-mu)*E + Q %*% E + Mus, (NoChnageLambdas-mu)*D + Q %*% D + Lambdas)
}

# intit for BiSSE
#       [,1]       [,2]
# [1,] 0.06884948 0.06884948
# [2,] 0.08884444 0.08884444
# [3,] 0.72557274 0.56003710
# [4,] 0.27442726 0.43996290

# init <- cbind(c(1:4), 6:9)
# pars.bisse <- c(1,  2,  0.03, 0.03,  5, 6)
# pars <- c(1,0,0,  0,0,2,  0.03, 0.03,  5, 6)
# diversitree:::initial.conditions.bisse(init, pars.bisse, t, idx)
# initial.conditions.classeNew(init, pars, t, idx)

initial.conditions.HiClasse_cpp <- function(init, pars, t, idx) {
  #print(init)
  E       <- init[c(1:MY_Nstates),1] # take the first branch arbitrarily
  D.left  <- init[c(-1:-MY_Nstates),1]
  D.right <- init[c(-1:-MY_Nstates),2]

  array <- cpp_pars_to_arrays(pars, MY_Nstates)
  lam.tensor <- array$lam.tensor

  Dmn = D.left %*% t(D.right)
  #flat_vec <- rep(1,MY_Nstates)
  #matrix_times_half <- OffDiagonalHalve(MY_Nstates)
  nodes <- lapply(lam.tensor, function(x) flat_vec %*% (makeSymmetricMatrix(x)*matrix_times_half*Dmn) %*% flat_vec )

  c(E, unlist(nodes))
}

# Root ClaSSE
# diversitree:::rootfunc.bisseness()
# pars <- c(1,0,0,  0,0,2,  0.03, 0.03,  5, 6)
# lik.new <- make.classeNew(tree, states,  sampling.f=NULL,  strict=TRUE, control=list(), newArgs)
# xx <- lik.new(pars, intermediates=T, root=ROOT.GIVEN, root.p=c(.5, .5))
# res <- attr(xx,"intermediates")
# root=ROOT.GIVEN
# root.p=c(0.5, 0.5)

rootfunc.HiClasse_cpp <- function (res, pars, condition.surv, root, root.p, intermediates){
  vals <- res$vals
  lq <- res$lq
  #d.root <- vals[3:4]
  d.root <- vals[-1:-MY_Nstates]
  root.p <- diversitree:::root_p_calc(d.root, pars, root, root.p, stationary.freq.bisseness)
  if (condition.surv) {
    #e.root <- vals[1:2]
    # lambda <- pars[1:2]
    e.root <- vals[1:MY_Nstates]
    #array <- pars_to_arrays(pars, MY_Nstates)
    array <- cpp_pars_to_arrays(pars, MY_Nstates)
    lam.tensor <- array$lam.tensor

    EE_root <- (1-e.root) %*% t(1-e.root)
    #matrix_times_half <- OffDiagonalHalve(MY_Nstates)
    #flat_vec <- rep(1,MY_Nstates)
    Nonextinct <- lapply(lam.tensor, function(x) flat_vec %*% (makeSymmetricMatrix(x)*matrix_times_half*EE_root) %*% flat_vec )
    Nonextinct <- unlist(Nonextinct)
    d.root <- d.root/sum(root.p * Nonextinct)
  }
  if (root == ROOT.ALL)
    loglik <- log(d.root) + sum(lq)
  else loglik <- log(sum(root.p * d.root)) + sum(lq)
  if (intermediates) {
    res$root.p <- root.p
    attr(loglik, "intermediates") <- res
    attr(loglik, "vals") <- vals
  }
  loglik
}


make.info.HiClasse_cpp <- function (phy, newArgs)
{
  list(name = newArgs$name, dll="hiclasse_c", name.pretty = "HiClaSSE", np = as.integer(length(argnames_HiClaSSE(newArgs$Nstates)$pars)),
       argnames = argnames_HiClaSSE(newArgs$Nstates)$pars,
       ny = newArgs$ny, k = newArgs$Nstates, idx.e = newArgs$idx.e, idx.d = newArgs$idx.d, derivs = derivs.HiClasse_cpp,

       phy = phy, ml.default = "subplex", mcmc.lowerzero = TRUE,
       doc = NULL, reference = c("no ref"))
}


make.cache.HiClasse_cpp <- function (tree, states, sampling.f = NULL, strict = TRUE, newArgs)
{
  cache <- diversitree:::make.cache.musse(tree, states + 1L, 2L, sampling.f, strict)
  cache$info <- make.info.HiClasse_cpp(tree, newArgs)
  cache
}


make.HiClasse_cpp <- function (tree, states, sampling.f = NULL, strict = TRUE, control = list(), newArgs=list()) {
  cache <- make.cache.HiClasse_cpp(tree, states, sampling.f, strict, newArgs)

  #-- custom states at tips
  # cache$y[[1]]$y <- c(0,0,0,0, 1,1,0,0)
  for (i in 1:length(cache$y)){
    # print(cache$y[[i]]$y)
    ca <- paste0(cache$y[[i]]$y, collapse = '')
    new <- paste0(newArgs$y[[i]], collapse = '')
    print(paste0('Recoding: ', ca, ' -> ', new) )
    cache$y[[i]]$y <- newArgs$y[[i]]
  }
  #control <- list(backend = "deSolve")
  all_branches <- diversitree:::make.all_branches.dtlik(cache, control, initial.conditions.HiClasse_cpp)
  #rootfunc <- diversitree:::rootfunc.geosse
  #rootfunc <- diversitree:::rootfunc.musse
  rootfunc <- rootfunc.HiClasse_cpp
  ll <- function(pars, condition.surv = TRUE, root = ROOT.OBS,
                 root.p = NULL, intermediates = FALSE) {
    #check.pars.geosse(pars)
    diversitree:::check.pars.nonnegative(pars, cache$info$np)
    ans <- all_branches(pars, intermediates)
    #diversitree:::all_branches(pars)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("dtlik", "function")
  ll
}

