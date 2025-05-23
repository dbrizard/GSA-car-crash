# Car crashworthiness global sensitivity analysis: Morris Analysis vs PCE and Sobol indices


## Prerequisites

This code was developped using Python 3.6.8 and the following modules:

- Numpy 1.19.4
- Matplotlib 3.36.2
- SALib 1.3.11
- openturns 1.19.post1
- dynareadout 23.10.2 **Not available any more. Can be replaced by *lasso-python* with just a few modifactions (not yet)**

LS-DYNA executable is `ls-dyna_smp_s_r1010_x64_redhat5_ifort160`. 


## LS-DYNA car crash model

The model is located in `lsopt_car6_v3` folder. 
The other `lsopt_xxx`folders contain older/developement versions of the model. 
There are 3 LS-DYNA files:

- the main file is `main_v223.k`. It contains parameters set at their initial values;
- the main file `main_v223_param.k` is used in case the values of the parameters are changed;
- the rest of the definition of the model is included in `car6_crash_v223.k`.


## Running the model from Python

Running the model from Python and performing the Global Sensitivity analysis is 
done with `LSDYNAmodel.py` module. 


## Post-treatment: plot all the Morris analyses

To plot the various Morris analyses, use `GSAutils.py` module. 
It reads and plots data stored manually in `GSA/car_v223_right-impact_v30/morris_nXX_output.md`,
where `XX` stands for the number of Morris trajectories used for each analysis.


## Polynomial Chaos metamodel and Sobol indices

Module `PCEsobol.py` builds a PCE metamodel and computes the Sobol indices. 
It relies on [OpenTURNS](https://openturns.github.io/openturns/latest/index.html). 

The input and output of the car crash model are generated with `LSDYNAmodel.py`,
they are stored in the `LHS` folder.

- `LHS-nnn_X.csv`, where `nnn` stands for the number of samples, contains the input samples;
- and `LHS-nnn_Y_xxxx.csv`, where `xxxx` stands for the name of the output, contains the output values.
