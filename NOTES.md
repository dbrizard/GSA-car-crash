# Side notes related to the project

## 20/09/2024 Test functions
### Available test functions
[Virtual Library of Simulation Experiments: Test Functions and Datasets](https://www.sfu.ca/~ssurjano/index.html)
[UQTestFuns](https://github.com/damar-wicaksono/uqtestfuns) requires Python 3.7, currently running 3.6.8...

SALib 1.3.11 has only: 

* SALib.test_functions.Ishigami
* SALib.test_functions.Sobol_G

OpenTURNS has:

* usecases.ackley_function.AckleyModel()
* usecases.branin_function.BraninModel()
* usecases.cantilever_beam.CantileverBeam()
* usecases.chaboche_model.ChabocheModel([...])
* usecases.deflection_tube.DeflectionTube()
* usecases.flood_model.FloodModel([trueKs, ...])
* usecases.ishigami_function.IshigamiModel()
* usecases.logistic_model.LogisticModel([t0, ...])
* usecases.stressed_beam.AxialStressedBeam()
* usecases.viscous_free_fall.ViscousFreeFall()
* usecases.fireSatellite_function.FireSatelliteModel()
* usecases.wingweight_function.WingWeightModel()
* usecases.oscillator.Oscillator()
* usecases.stiffened_panel.StiffenedPanel()


### Why not do things in R??
We have the `sensitivity` package...


### GSA toolboxes in Python

* [SALib](https://salib.readthedocs.io/en/latest/)
* [OpenTURNS](https://openturns.github.io/openturns/latest/theory/reliability_sensitivity/reliability_sensitivity.html#sensitivity-analysis)
* [SAFE toolbax](https://safetoolbox.github.io/) 
  + requires python 3.7
  + documentation is rather scarce, based on notebooks

### CONCLUSION
Should do first tests with SALib.

## 24/09/2024
Following biblio from `UQTestFuns` JOSS paper, installed with pip
- `Successfully installed pyDOE2-1.3.0 **smt-1.3.0**`

`UQTestFuns` installed in `python38` venv
Did not install `pyapprox` as it requires `torch` and is very heavy in memory (too much download !!)

## 25/09/2024
Number of interesting functions in `UQTestFuns` is limited if I only want uniform distribution inputs.
16 GSA related functions. 

## 29/04/2025
OT was updated from 1.19 to 1.24 for metamodel validation (following SAMO 2025 conversation).

**WARNING** I think bootstrap of Sobol indices was broken with the update.

* [ ] try PCE metamodel with 1000 QMC samples (should be better than LHS!)
  + use LowDiscrepancyExperiment (suite de Sobol)
* [ ] try Morris analysis with Openturns (otmorris)
