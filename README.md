# Car crashworthiness global sensitivity analysis for USD 2024 


## Prerequisites

This code was developped using Python 3.6.8 and the following modules:

- Numpy 1.19.4
- Matplotlib 3.36.2
- SALib 1.3.11
- openturns 1.19.post1

LS-DYNA executable is `ls-dyna_smp_s_r1010_x64_redhat5_ifort160`. 


## LS-DYNA car crash model

The model is located in `lsopt_car6_v3` folder.
There are 3 LS-DYNA files:

- the main file is `main_v223.k`. It contains parameters set at their initial values;
- the main file `main_v223_param.k` is used in case the values of the parameters are changed;
- the rest of the definition of the model is included in `car6_crash_v223.k`.


## Running the model from Python

Running the model from Python and performing the Global Sensitivity analysis is 
done with `LSDYNAmodel.py` module. 


## Post-treatment: plot all the Morris analyses

To plot the various Morris analyses, use `GSAutils.py` module. 
It reads and plot data stored manually in `GSA/car_v223_right-impact_v30/morris_n10_output.md`.


## Polynomial Chaos metamodel and Sobol indices

Module `PCEsobol.py` builds a PCE metamodel and computes the Sobol indices. 
It relies on [OpenTURNS](https://openturns.github.io/openturns/latest/index.html). 

The input and output of the car crash model are generated with `LSDYNAmodel.py`,
they are stored in the `LHS4Eric_X.csv` and `LHS4Eric_Y_xxxx.csv` files. 