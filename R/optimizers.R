### NSGA-II -----------------------------------------------------------------------------
run_nsga2 <- function(config, control) {
  
  pop.size   <- control$pop.size
  no.iters   <- control$no.iters
  cross.prob <- control$cross.prob
  mut.prob   <- control$mut.prob
  no.cycles  <- control$no.cycles
  no.weights <- 4 * (config$no.floors - 1)
  max.skips  <- 4 
  
  if (requireNamespace("memoise", quietly = TRUE)){
    tryCatch({
      sring <- memoise::memoise(sring)
    }, error = function(err) {
    })
  }
  
  objs <- function(w) sring(config, w, no.cycles)[c(1, 2)]
  cons <- function(w) max.skips - sring(config, w, no.cycles)[3]
  both <- function(w) sring(config, w, no.cycles)
  
  if (requireNamespace("mco", quietly = TRUE)){
    res <- mco::nsga2(fn           = objs, 
                      idim         = no.weights, 
                      odim         = 2, 
                      lower.bounds = rep(-1, no.weights), 
                      upper.bounds = rep(1, no.weights), 
                      constraints  = cons,
                      cdim         = 1, 
                      popsize      = pop.size,
                      cprob        = cross.prob,
                      mprob        = mut.prob,
                      generations  = 1:no.iters)
    pf <- find_pf(res, both)
  } else {
    pf <- NULL
  }
  
  
  return(pf)
}

### DEMO --------------------------------------------------------------------------------
run_demo <- function(config, control) {
  
  pop.size   <- control$pop.size
  no.iters   <- control$no.iters
  cross.prob <- control$cross.prob
  mut.prob   <- control$mut.prob
  no.cycles  <- control$no.cycles
  no.weights <- 4 * (config$no.floors - 1)
  max.skips  <- 4 
  
  if (requireNamespace("memoise", quietly = TRUE)){
    tryCatch({
      sring <- memoise::memoise(sring)
    }, error = function(err) {
    })
  }
  
  objs <- function(w) {
    out <- sring(config, w, no.cycles)
    out[3] <- max(out[3] - max.skips, 0)
    return(out)
  }
  
  both <- function(w) sring(config, w, no.cycles)
  
  res <- demo(fn         = objs, 
              lower      = rep(-1, no.weights), 
              upper      = rep(1, no.weights), 
              pop.size   = pop.size,
              no.iter    = no.iters,
              cross.prob = cross.prob,
              scal.fac   = mut.prob,
              no.cons    = 1)
  
  pf <- find_pf(res, both)
  
  return(pf)
}

### MOEA/D ------------------------------------------------------------------------------
run_moead <- function(config, control, seed) {
  
  pop.size   <- control$pop.size
  no.iters   <- control$no.iters
  cross.prob <- control$cross.prob
  mut.prob   <- control$mut.prob
  no.cycles  <- control$no.cycles
  no.weights <- 4 * (config$no.floors - 1)
  max.skips  <- 4 
  
  if (requireNamespace("memoise", quietly = TRUE)){
    tryCatch({
      sring <- memoise::memoise(sring)
    }, error = function(err) {
    })
  }
  
  fn.objs <- function(w) sring(config, w, no.cycles)[c(1, 2)]
  fn.cons <- function(w) max.skips - sring(config, w, no.cycles)[3] - max.skips
  both <- function(w) sring(config, w, no.cycles)
  
  pop.x.archive <<- lapply(1:no.iters, function(x) list(par = NULL))
  i <<- 1
  
  objs <<- function(X) {
    pop.x.archive[[i]]$par <<- X
    i <<- i + 1
    out <- t(apply(X, 1, fn.objs))
    return(out)
  }
  
  cons <<- function(X) {
    Cmatrix <- matrix(numeric(), 
                      nrow = nrow(X),
                      ncol = 2 * no.weights + 1)
    Xmin <- matrix(-1, nrow = nrow(X), ncol = no.weights)
    Xmax <- matrix(1, nrow = nrow(X), ncol = no.weights)
    Cmatrix[, 1:no.weights] <- Xmin - X
    Cmatrix[, (no.weights + 1):(2 * no.weights)] <- X - Xmax
    Cmatrix[, (2 * no.weights + 1)] <- as.matrix(apply(X, 1, fn.cons))
    Vmatrix <- pmax(Cmatrix, 0)  
    return(list(Cmatrix = Cmatrix,
                Vmatrix = Vmatrix,
                v       = rowSums(Vmatrix)))
  }
  
  problem    <- list(name = "objs", xmin = rep(-1, no.weights), 
                     xmax = rep(1, no.weights), m = 2, 
                     constraints = list(name = "cons"))
  decomp     <- list(name = "SLD", H = (pop.size - 1))
  neighbors  <- list(name = "lambda", T = 20, delta.p = 1)
  aggfun     <- list(name = "wt")
  variation  <- list(list(name = "sbx", etax = 20, pc = cross.prob),
                     list(name = "polymut", etam = 20, pm = mut.prob),
                     list(name = "truncate"))
  update     <- list(name = "standard", UseArchive = TRUE)
  scaling    <- list(name = "none")
  constraint <- list(name = "penalty", beta = 1)
  stopcrit   <- list(list(name = "maxiter", maxiter = no.iters - 1)) 
  showpars   <- list(show.iters = "none")
  
  if (requireNamespace("MOEADr", quietly = TRUE)){
    res <- MOEADr::moead(preset     = NULL,
                         problem    = problem,
                         decomp     = decomp,
                         neighbors  = neighbors,
                         aggfun     = aggfun,
                         variation  = variation,
                         update     = update,
                         constraint = constraint,
                         scaling    = scaling,
                         stopcrit   = stopcrit,
                         showpars   = showpars,
                         seed       = seed)
    
    pf <- find_pf(pop.x.archive, both)
  } else {
    pf <- NULL
  }
  
  return(pf)
}

### Helpers -----------------------------------------------------------------------------

find_pf <- function(res, fn) {
  res <- lapply(1:length(res), function(i) res[[i]]$par)
  res <- do.call("rbind", res)
  res <- t(apply(res, 1, fn)) 
  res <- subset(res, res[, 3] < 5)[, 1:2]
  res <- t(emoa::nondominated_points(t(as.matrix(res))))
  res <- data.frame(x = res[, 1], y = res[, 2])
  return(res)
}
