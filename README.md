## Elevator Group Control as a Constrained Multiobjective Optimization Problem

**NOTE:** This repo contains the supplementary material to the paper:

> [1] A. Vodopija, J. Stork, T. Bartz-Beielstein, B. Filipič. Elevator group control as a constrained multiobjective optimization problem, *Applied Soft Computing*.

It's recommended to read this paper for a better understanding of the following content.

### 1. SETUP
#### 1.1. Prerequisites

The `devtools` and `Rcpp` packages, as well as an appropriate `C++` compiler, are required to install the package. In addition, the `emoa` package is mandatory to run the S-Ring simulations and to experiment with Differential Evolution for Multiobjective Optimization (DEMO).    

For other tasks you need: 
- `mco`: To experiment with Nondominated Sorting Genetic Algorithm II (NSGA-II).
- `MOEADr`: To experiment with Multiobjective Optimization Evolutionary Algorithm based on Decomposition (MOEA/D). 
- `SPOT`: For parameter tuning with Sequential Parameter Optimization Toolbox (SPOT).

For more efficent computations and nicer plots:
- `memoise`: To efficiently execute the S-Ring simulations.
- `ggplot2`: For nicer plots. 

#### 1.2. Installation

The package can be installed by running the following cell (`devtools` is required). 


```R
devtools::install_github("https://github.com/vodopijaaljosa/egccmop")
```

The installation procedure will automatically install `Rcpp` and `emoa`. Other packages are optional and can be installed using the command `install.packages("pkg")`. See Section 1.1 to select the packages you are interested in.

Once the installation is completed, the package functionalities can be imported by running the following cell. 


```R
library(egccmop)
```

You're done! You can start experimenting.

### 2. S-RING

In this part, you'll run an S-Ring simulation to compute the objective and constraint values, but first, let's see how to define a new elevator group control (EGC) configuration. For more details on S-Ring see Section 3 in [1].

#### 2.1. Elevator group control configurations

You can define a custom EGC configuration as shown below. The following parameters are mandatory:
- `no.floors`: An integer denoting the number of floors. It is reasonable to pick up to 50 floors. For larger configurations the simulation becomes time consuming. 
- `no.elevators`: An integer denoting the number of elevators. It should be smaller than the number of floors.
- `prob`: A float from [0,0.3] representing the probability of a newly arriving customer.


```R
config <- list(no.floors    = 5,
               no.elevators = 3,  
               prob         = 0.1)
```

#### 2.2. Computation of objective and constraint values

Let's compute the objective and constraint values for a given solution (perceptron weights) of the EGC problem defined above. 

The input parameters for the function performing the S-Ring simulations are:
- `config`: A list representing an elevator group control configuration (see Section 2.1 of this notebook).
- `weights`: A vector of perceptron weights. Vector components take values from [-1,1].
- `no.cycles`: An integer denoting the number of cycles used in the S-Ring simulation. It should take values from  [1e3,1e5]. Default is 1e4.

The first element of the output is the proportion of states with waiting customers ($f_1$), the second the proportion of elevator stops ($f_2$), and the last one the maximal number of skips ($c$). For more details on objectives and constraints see Section 4 in [1]. 


```R
weights <- c(-0.8, -0.9, 0.6, -0.3, 0.9, 0.1, 0.6, 0.3, 1.0, -0.8, -0.5, 0.0, -0.1, 0.5, 0.7, 0.3)
objs    <- sring(config    = config, 
                 weights   = weights,
                 no.cycles = 1e4)
print(objs)
```

### 3. OPTIMIZATION

Now, let's see how to run an optimization procedure. See [1], Section 5.2 for more details on the optimization algorithms and their settings.

#### 3.1. Parameter settings

First, let's select the values of control parameters, i.e., the algorithm parameters and the number of cycles used in the S-Ring simulations.
- `pop.size`: An integer denoting the number of solutions used by the optimization algorithm. Default is 100.
- `no.iters`: An integer denoting the number of iterations (generations) used by the optimization algorithm. Default is 100.
- `cross.prob`: A float from [0,1] denoting the crossover probability. Default is 0.9. 
- `mut.prob`: A float from [0,1] denoting the mutation probability. For DEMO this is the scaling factor ($F$). Default is 0.1 for NSGA-II and MOEA/D, while 0.7 for DEMO.
- `no.cycles`: An integer denoting the number of cycles used in the S-Ring simulation. It should take values from  [1e3,1e5]. Default is 1e4.

**NOTE:** Default values are used for unspecified parameters.  


```R
control <- list(pop.size   = 100, 
                no.iters   = 50,
                cross.prob = 0.9,
                mut.prob   = 0.1,
                no.cycles  = 1e3)
```

#### 3.2. Run the optimization

The function `opt` performs the optimization procedure.

Its input parameters are:
- `config`: A list representing an elevator group control configuration (see Section 2.1 of this notebook).
- `method`: A string indicating the used optimization algorithm. It can be NSGA-II ("nsga2"), DEMO ("demo") or MOEA/D ("moead"). Default is DEMO.
- `control`: A list of control parameters (see Section 3.1 of this notebook).
- `no.runs`: An integer denoting the number of runs. Default is 1. 

The output of the function is a list of data frames storing the objective values of all nonodominated feasible solutions found during each run.


```R
pfs <- opt(config  = config, 
           method  = "nsga2", 
           control = control,
           no.runs = 3)
```

### 4. HYPERVOLUME STATISTICS

The following function summarizes the results. First, it computes the hypervolume values of all runs (all nondominated feasible solutions found during the entire run are used). Then, it derives six statistics: Minimum, 1st Quartile, Median, Mean, 3rd Quartile, Maximum. 


```R
compute_hv_stats(pfs)
```

### 5. PARETO FRONT APPROXIMATIONS

You can plot the obtained Pareto front approximations by running the following cell. The `run` parameter is used to select which run to depict in the plot. 


```R
make_plot(pfs, run = 2)
```

### 6. PARAMETER TUNING

**NOTE:** The `SPOT` package is required for this section. Let's install it by running the following cell.


```R
install.packages("SPOT")
```

#### 6.1. Tuning algorithm parameters

This section shows how to tune the algorithm parameters. See [1], Section 5.2 for more details on parameter tuning.

The input parameters of the tuning process are:
- `config`: A list representing an elevator group control configuration (see Section 2.1 of this notebook).
- `method`: A string indicating the used optimization algorithm. It can be NSGA-II ("nsga2"), DEMO ("demo") or MOEA/D ("moead"). Default is DEMO.
- `no.cycles`: An integer denoting the number of cycles used in the S-Ring simulation. It should take values from  [1e3,1e5]. Default is 1e4.
- `spot.config`: A list of SPOT settings. See the SPOT documentation for more details (<https://cran.r-project.org/web/packages/SPOT/index.html>). 

The following SPOT settings were used in the experiments reported in [1]. 

**NOTE:** Here, the parameter `funEvals` is reduced for faster computations!


```R
spot.config <- list(types            = c("numeric", "numeric", "numeric"),
                    funEvals         = 15, # In the paper experiments 'funEvals' was set to 50 
                    noise            = TRUE,
                    seedFun          = 1,
                    seedSPOT         = 1,
                    replicates       = 1,
                    design           = designLHD,
                    model            = buildKriging,
                    optimizer        = optimDE,
                    optimizerControl = list(funEvals = 1000))
```

Run the following cell to tune the parameters.


```R
control.tuned <- param_tuning(config      = config, 
                              method      = "demo",
                              no.cycles   = 1e3,
                              spot.config = spot.config)
```

#### 6.2. Re-run the optimization using the tuned parameters

Let's re-run the optimization with the tuned parameters. 


```R
pfs.tuned <- opt(config  = config, 
                 method  = "demo", 
                 control = control.tuned,
                 no.runs = 3)
```

Finally, let's compare the results obtained before and after the tuning. 


```R
hv.stats           <- rbind(compute_hv_stats(pfs), compute_hv_stats(pfs.tuned)) 
rownames(hv.stats) <- c("before", "after")
print(hv.stats)
```
