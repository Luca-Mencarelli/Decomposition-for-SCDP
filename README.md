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

In particular:
* ``main-random-random.jl`` implements the pure random strategy to generate and adding new orbits to the master problem;
* ``main-random-v2iter.jl`` implements the DS-1 strategy to generate and adding new orbits to the master problem;
* ``main-random-v3iter.jl`` implements the DS-2 strategy to generate and adding new orbits to the master problem; 
