# Decoding results
spotDecoding <- function(algpar, max.eval) {
  pop.size   <- algpar[1] 
  pop.size   <- round((1 - pop.size) * 50 + pop.size * 500)
  pop.size   <- floor(pop.size / 4) * 4
  no.iters   <- floor(max.eval / pop.size)
  cross.prob <- algpar[2] 
  mut.prob   <- algpar[3]
  return(list(pop.size = pop.size, no.iters = no.iters, cross.prob = cross.prob, mut.prob = mut.prob))
}

# SPOT interface
methodToSpot <- function(method, config, algpar, max.eval) {
  performance <- NULL
  
  for (i in 1:nrow(algpar)) {
    pop.size <- algpar[i, 1] 
    pop.size <- round((1 - pop.size) * 50 + pop.size * 500)
    pop.size <- floor(pop.size / 4) * 4 
    no.iters  <- floor(max.eval / pop.size)
    
    control <- list(pop.size   = pop.size, 
                    no.iters   = no.iters,
                    cross.prob = algpar[i, 2],
                    mut.prob   = algpar[i, 3],
                    no.cycles  = 1e3)
      
    res <- opt(config, 
               method  = method, 
               control = control, 
               no.runs = 1)[[1]]
    
    res <- emoa::dominated_hypervolume(t(as.matrix(res)), rep(1.1, 2))
    
    performance <- c(performance, -res)
  }
  return(matrix(performance, ncol = 1))
}