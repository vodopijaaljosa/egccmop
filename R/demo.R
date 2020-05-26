### DEMO algorithm ----------------------------------------------------------------------

demo <- function(fn, 
                 lower, 
                 upper, 
                 pop.size = 100, 
                 no.iter = 100, 
                 cross.prob = 0.7, 
                 scal.fac = 0.5,
                 no.cons = 0) {
  
  # population initialization, have to be modifed for other problem
  pop.x <- popInit(pop.size, lower, upper)
  
  # population evaluation
  pop.y <- t(apply(pop.x, 1, fn))
  
  # initialize archive
  pop.x.archive <- lapply(1:no.iter, function(x) list(par = NULL))
  pop.x.archive[[1]]$par <- pop.x
  
  # main loop  
  for (i in 2:no.iter) {
    pop.x.new <- deRand1Bin(pop.x, lower, upper, scal.fac)
    pop.x.new <- crossover(pop.x.new, pop.x, cross.prob)
    pop.y.new <- t(apply(pop.x.new, 1, fn))
    
    sort <- nsga2CdSelection(pop.x, pop.y, pop.x.new, pop.y.new, no.cons = no.cons)
    pop.x <- sort$x
    pop.y <- sort$y
    pop.x.archive[[i]]$par <- pop.x
  }
  
  return(pop.x.archive)
}

### Helpers -----------------------------------------------------------------------------

# Random initialization
popInit <- function(pop.size, lower, upper) {
  d <- length(lower)
  pop.x <- t(replicate(pop.size, runif(d, lower, upper)))
  return(pop.x)
}

# DE/rand/1/bin
deRand1Bin <- function(pop.x, lower, upper, scaling.factor) {
  pop.x.new <- t(sapply(1:nrow(pop.x), function(ind.i) mutDE(ind.i, pop.x, scaling.factor)))
  pop.x.new <- t(apply(pop.x.new, 1, function(ind.x) ensCons(ind.x, lower, upper)))
  return(pop.x.new)  
}

mutDE <- function(ind.i, pop.x, scaling.factor) {
  inds <- sample(1:nrow(pop.x), 4)
  inds <- inds[!inds %in% ind.i]
  ind.x.new <- pop.x[inds[1], ] + scaling.factor * (pop.x[inds[2], ] - pop.x[inds[3], ])  
  return(ind.x.new)
}

ensCons <- function(ind.x, lower, upper) {
  d <- length(lower)
  ind.x.new <- runif(d, lower, upper)
  ind.x[ind.x > upper] <- ind.x.new[ind.x > upper]
  ind.x[ind.x < lower] <- ind.x.new[ind.x < lower]
  return(ind.x)
}

# Crossover
crossover <- function(pop.x.new, pop.x, cross.prob) {
  pop.x.new <- t(sapply(1:nrow(pop.x.new), function(ind.i) crossDE(ind.i, pop.x.new, pop.x, cross.prob)))
  return(pop.x.new)
}

crossDE <- function(ind.i, pop.x.new, pop.x, cross.prob) {
  d <- ncol(pop.x)
  j.rand <- sample(1:d, 1)
  j.new  <- runif(d, 0, 1) <= cross.prob
  ind.x  <- pop.x[ind.i, ]
  ind.x.new <- pop.x.new[ind.i, ]
  ind.x[j.new] <- ind.x.new[j.new]
  ind.x[j.rand] <- ind.x.new[j.rand]
  return(ind.x)
}

# Survivor selection
nsga2CdSelection <- function(pop.x, 
                             pop.y, 
                             pop.x.new, 
                             pop.y.new,
                             no.cons = 0) {
  
  pop.x.all <- rbind(pop.x, pop.x.new)
  pop.y.all <- rbind(pop.y, pop.y.new)
  
  list.1 <- sample(1:nrow(pop.x.all), nrow(pop.x))
  list.2 <- sample(1:nrow(pop.x.all), nrow(pop.x))
  keep <- c() 
  
  for (i in 1:nrow(pop.x)) {
    i1 <- list.1[i]
    i2 <- list.2[i]
    y1 <- pop.y.all[i1, ]
    y2 <- pop.y.all[i2, ]
    y1.obj <- head(y1, length(y1) - no.cons)
    y2.obj <- head(y2, length(y2) - no.cons)
    y1.phi <- sum(tail(y1, no.cons))
    y2.phi <- sum(tail(y2, no.cons))
    
    if (y1.phi < y2.phi) {
      keep <- c(keep, i1)
    } else if (y2.phi < y1.phi) {
      keep <- c(keep, i2)
    } else {
      dominated <- emoa::is_dominated(matrix(c(y1.obj, y2.obj), length(y1.obj)))
      if (dominated[1]) {
        keep <- c(keep, i2)
      } else {
        keep <- c(keep, i1)
      }
    }
  }
  
  pop.x.next <- pop.x.all[keep, ]
  pop.y.next <- pop.y.all[keep, ]
  
  return(list(x = pop.x.next, y = pop.y.next))
}

