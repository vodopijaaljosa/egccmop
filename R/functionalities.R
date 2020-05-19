### 0. Prerequisits ---------------------------------------------------------------------

### 1. S-Ring ---------------------------------------------------------------------------

#' S-Ring
#' 
#' This function run an S-Ring simulation.
#' 
#' @param config A list representing an elevator group control configuration:
#' \describe{
#'		\item{\code{no.floors}}{An integer denoting the number of floors. 
#'		It is reasonable to pick up to 50 floors. For larger configurations 
#'		the simulation becomes time consuming.}
#'		\item{\code{no.elevators}}{An integer denoting the number of elevators. 
#'		It should be smaller than the number of floors.}
#'		\item{\code{prob}}{A float from [0,0.3] representing the probability of 
#'		a newly arriving customer.}
#' }
#' @param weights A vector of perceptron weights. Vector components take values from [-1,1].
#' @param no.cycles An integer denoting the number of cycles used in the S-Ring simulation. 
#' It should take values from  [1e3,1e5]. Default is 1e4.
#' 
#' @return An vector of normalized objective values: 1) Proportion of state 
#' with waiting customers, 2) Proportion of elevator stops, 3) Maximal number 
#' of elevator skips.
#' 
#' @export
sring <- function(config, weights, no.cycles = 1e4) {
  
  if ((!no.cycles%%1 == 0) || (no.cycles < 1e3) || (no.cycles > 1e5)) {
    stop("'no.cycles' is not an integer or no.cycles not in [1e3,1e5]")
  }
  
  if ((config$prob < 0) || (config$prob > 0.3)) {
    stop("'prob' not in [1e3,1e5]")
  }
  
  return(s_ring(weights      = weights, 
                no_elevators = config$no.elevators, 
                prob         = config$prob, 
                alpha        = 1.5, 
                beta         = 1.5, 
                no_cycles    = no.cycles))
}

### 2. Optimization ---------------------------------------------------------------------

#' Optimization
#' 
#' This function solves the elevator group control problem using a 
#' Multiobjective Optimization Evolutionary Algorithm (MOEA).
#' 
#' @param config A list representing an elevator group control configuration:
#' \describe{
#'		\item{\code{no.floors}}{An integer denoting the number of floors. 
#'		It is reasonable to pick up to 50 floors. For larger configurations 
#'		the simulation becomes time consuming.}
#'		\item{\code{no.elevators}}{An integer denoting the number of elevators. 
#'		It should be smaller than the number of floors.}
#'		\item{\code{prob}}{A float from [0,0.3] representing the probability of 
#'		a newly arriving customer.}
#' }
#' @param method A string indicating the used optimization algorithm. 
#' It can be NSGA-II ("nsga2"), DEMO ("demo") or MOEA/D ("moead"). Default is DEMO.
#' @param control A list of control parameters: 
#' #' \describe{
#'		\item{\code{pop.size}}{An integer denoting the number of solutions 
#'		used by the optimization algorithm. Default is 100.}
#'		\item{\code{no.iters}}{An integer denoting the number of iterations 
#'		(generations) used by the optimization algorithm. Default is 100.}
#'		\item{\code{cross.prob}}{A float from [0,1] denoting the crossover 
#'		probability. Default is 0.9.}
#'		\item{\code{mut.prob}}{A float from [0,1] denoting the mutation probability. 
#'		For DEMO this is the scaling factor. Default is 0.1 for NSGA-II and MOEA/D, 
#'		while 0.7 for DEMO.}
#'		\item{\code{no.cycles}}{An integer denoting the number of cycles used in 
#'		the S-Ring simulation. 
#'		It should take values from  [1e3,1e5]. Default is 1e4.}
#' }
#' Default values are used for unspecified parameters. 
#' @param no.runs An integer denoting the number of runs. Default is 1. 
#' 
#' @return List of data frames storing the objective values of all 
#' nonodominated feasible solutions found during each run.
#' 
#' @export
opt <- function(config, 
                method  = c("nsga2", "demo", "moead"), 
                control = NULL, 
                no.runs = 1) {
  
  method <- match.arg(method)
  
  if ((!no.runs%%1 == 0) || (no.runs < 1) || (no.runs > 30)) {
    stop("'no.runs' not an integer or no.runs not in [1,30]")
  }
  
  control <- list(no.cycles  = ifelse("no.cycles" %in% names(control), control$no.cycles, 1e4), 
                  pop.size   = ifelse("pop.size" %in% names(control), control$pop.size, 100), 
                  no.iters   = ifelse("no.iters" %in% names(control), control$no.iters, 100), 
                  cross.prob = ifelse("cross.prob" %in% names(control), control$cross.prob, 0.9), 
                  mut.prob   = ifelse("mut.prob" %in% names(control), control$mut.prob, 0.1))
  
  if ((control$mut.prob < 0) || (control$mut.prob > 1)) {
    stop("'mut.prob' not in [0,1]")
  }
  
  if ((control$cross.prob < 0) || (control$cross.prob > 1)) {
    stop("'cross.prob' not in [0,1]")
  }
  
  pfs <- list()
  
  if (method == "nsga2") {
    for (i in 1:no.runs) {
      if (!requireNamespace("mco", quietly = TRUE)) stop("'mco' package is required", call. = FALSE)
      pfs[[paste0("run.", i)]] <- run_nsga2(config  = config,
                                            control = control)
    }  
    
  } else if (method == "demo") {
    for (i in 1:no.runs) {
      pfs[[paste0("run.", i)]] <- run_demo(config  = config,
                                           control = control)
    } 
    
  } else if (method == "moead") {
    if (!requireNamespace("MOEADr", quietly = TRUE)) stop("'MOEADr' package is required", call. = FALSE)
    for (i in 1:no.runs) {
      pfs[[paste0("run.", i)]] <- run_moead(config  = config,
                                            control = control,
                                            seed    = i)
    } 
  }
  
  return(pfs)
}

### 3. Statistics -----------------------------------------------------------------------

#' Computing hypervolume statisitcs
#' 
#' This function compute hypervolume statistics from runs.
#' 
#' @param pfs A list of Pareto front approximations obtained by \code{opt} function.
#' 
#' @return A vector of Min, 1st Qu., Median, Mean, 3rd Qu., Max. of the hypervulume values. 
#' 
#' @export
compute_hv_stats <- function(pfs) {
  stats <- sapply(pfs, function(pf) emoa::dominated_hypervolume(t(as.matrix(pf)), rep(1.1, 2)))
  return(summary(stats))
}

### 4. Plots ----------------------------------------------------------------------------

#' Making plots
#' 
#' This function meka a plot of the selected Pareto front approximation.
#' 
#' @param pfs A list of Pareto front approximations obtained by \code{opt} function.
#' @param run An integer denoting which run to depict in the plot.  
#' 
#' @return A plot depicting Pareto front approximation of the selected run. 
#' 
#' @export
make_plot <- function(pfs, run = 1) {
  if ((!run%%1 == 0) || (run < 1) || (run > length(pfs))) {
    stop("'run' not an integer or run not in [1,length(pfs)]")
  }
  pf <- pfs[[paste0("run.", run)]]
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(pf, ggplot2::aes(x = x, y = y)) +
      ggplot2::theme_bw() +
      ggplot2::xlab("f1: Proportion of states with waiting customers") +
      ggplot2::ylab("f2: Proportion of elevator stops") +
      ggplot2::geom_point()
  } else {
    plot(pf,
      xlab = "f1: Proportion of states with waiting customers",
      ylab = "f2: Proportion of elevator stops")
  }
}

### 5. Parameter tuning -----------------------------------------------------------------

#' Parameter tuning
#' 
#' This function is used to tune algorithm parameters. The SPO Toolbox is 
#' used for this purpose.
#' 
#' @param config A list representing an elevator group control configuration:
#' \describe{
#'		\item{\code{no.floors}}{An integer denoting the number of floors. 
#'		It is reasonable to pick up to 50 floors. For larger configurations 
#'		the simulation becomes time consuming.}
#'		\item{\code{no.elevators}}{An integer denoting the number of elevators. 
#'		It should be smaller than the number of floors.}
#'		\item{\code{prob}}{A float from [0,0.3] representing the probability of 
#'		a newly arriving customer.}
#' }
#' @param method A string indicating the used optimization algorithm. 
#' It can be NSGA-II ("nsga2"), DEMO ("demo") or MOEA/D ("moead"). Default is DEMO.
#' @param no.cycles An integer denoting the number of cycles used in the S-Ring simulation. 
#' It should take values from  [1e3,1e5]. Default is 1e4.
#' @param spot.config A list of SPOT settings. See the SPOT documentaion for more details.
#' 
#' @return The control list of tuned parameters. 
#' 
#' @export
param_tuning <- function(config, 
                         method      = c("nsga2", "demo", "moead"), 
                         no.cycles   = 1e4, 
                         spot.config = NULL) {
  
  method <- match.arg(method)
  
  if (requireNamespace("MOEADr", quietly = TRUE)) {
  
    if (is.null(spot.config)) {
      spot.config <- list(types            = c("numeric", "numeric", "numeric"),
                          funEvals         = 15,
                          noise            = TRUE,
                          seedFun          = 1,
                          seedSPOT         = 1,
                          replicates       = 1,
                          design           = SPOT::designLHD,
                          model            = SPOT::buildKriging,
                          optimizer        = SPOT::optimDE,
                          optimizerControl = list(funEvals = 1000))
    }
    
    max.eval <- 2000 * config$no.floors 
    fun <- function(algpar) methodToSpot(method, config, algpar, max.eval)
    
    res <- SPOT::spot(x       = NULL,
                      fun     = fun,
                      lower   = rep(0, 3),
                      upper   = rep(1, 3),  
                      control = spot.config)
    
    control.tuned <- spotDecoding(res$xbest, max.eval)
    control.tuned$no.cycles <- no.cycles
    
    return(control.tuned)
  } else {
    stop("'SPOT' package is required", call. = FALSE)
  }
} 
