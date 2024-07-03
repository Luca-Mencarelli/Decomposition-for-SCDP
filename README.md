# Decomposition-for-SCDP

This repository contains the Julia scripts to replicate results in L. Mencarelli and J. Floquet, Decomposition based heuristic approaches for the
satellite constellation design problem, working paper, 2024. 

## Requirements

The packages required to run the code are ``JuMP, SCDP, StaticArrays, SatelliteToolbox, Random, JLD2, CPLEX, DataStructures, Ipopt, Distributions``. To install requirements:

```
import Pkg; Pkg.add("<package_name>")
```

where <package_name> is one of the previous packages.

## Code

``main-random*.jl`` files contain the main Julia scripts and to generate instances and obtain the corresponding results (contains in the directory ``Results/Global``). In order to run the code, it is sufficient to type in Terminal:

```
julia main*.jl
```

In particular the strategies implemented to generate and add new orbits to the master problem, are:
* the pure random strategy in ``main-random-random.jl``;
* the DS-1 strategy ``main-random-v2iter.jl``;
* the DS-2 strategy ``main-random-v3iter.jl``;
