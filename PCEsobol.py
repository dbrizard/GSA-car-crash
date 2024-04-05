#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 11:03:36 2024

@author: dbrizard
"""

import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import pylab as plt
import openturns as ot
import openturns.viewer as viewer
ot.Log.Show(ot.Log.NONE)

import pythia as pt


def sobol_function(x, a=None, **kwargs) -> np.ndarray:
    """Sobol function.

    Parameters
    ----------
    x : np.ndarray
        Input values.
    a : np.ndarray | None, optional
        Coefficients, by default None

    Returns
    -------
    :
        Sobol function values.

    Raises
    ------
    ValueError
        Wrong dimension for input `x`.
    ValueError
        Wrong shape of Coefficients `a`.
    """
    if not 0 < x.ndim < 3:
        raise ValueError("Wrong ndim: {}".format(x.ndim))
    if x.ndim == 1:
        x.shape = 1, -1
    if a is None:
        a = np.zeros(x.shape[1])
    elif not a.shape == (x.shape[1],):
        raise ValueError("Wrong shape: {}".format(a.shape))
    return np.prod((abs(4.0 * x - 2.0) + a) / (1.0 + a), axis=1).reshape(-1, 1)


def sobol_sc(a: np.ndarray, dim: int = 1, **kwargs):
    """Sobol function Sobol indices.

    Parameters
    ----------
    a : np.ndarray
        Coefficients.
    dim : int, optional
        Parameter dimension, by default 1.

    Returns
    -------
    :
        Sobol indices of Sobol function.
    """
    sobol = {}
    beta = (1.0 + a) ** (-2) / 3
    var = np.prod(1.0 + beta) - 1.0
    sobol_tuples = pt.index.IndexSet(pt.index.tensor_set([1, 1, 1])).sobol_tuples
    for sdx in sobol_tuples:
        sobol[sdx] = 1.0 / var
        for k in sdx:
            sobol[sdx] *= beta[k - 1]
    if dim > 1:
        return np.array([sobol for _ in range(dim)])
    else:
        return sobol





def target_function(x: np.ndarray) -> np.ndarray:
    """Target function.

    Parameters
    ----------
    x : np.ndarray
    """
    return sobol_function(x, a=a)


if __name__=='__main__':
    plt.close('all')
    
    #%% Openturns on LS-DYNA car simulation data
    if False:
        ot.ResourceMap.SetAsUnsignedInteger("FittingTest-LillieforsMaximumSamplingSize", 100)
        
        # Import X and Y
        problem = {'names': ['tbumper', 'trailb', 'trailf', 'tgrill', 'thood', 'ybumper', 'yrailf', 'yrailb', 'ybody'], 
                   'units': ['mm', 'mm', 'mm', 'mm', 'mm', 'MPa', 'MPa', 'MPa', 'MPa'], 
                   'num_vars': 9, 
                   'bounds': [[2, 4], [1, 3], [3, 7], [0.5, 1.5], [0.5, 1.5], [300, 500], [300, 500], [300, 500], [300, 500]]}

        X = np.loadtxt('LHS4Eric_X.csv', delimiter=',')
        inputSample = ot.Sample(X)
        inputSample.setDescription(problem['names'])
        
        Y = {'dmax':np.loadtxt('LHS4Eric_Y_dmax.csv', delimiter=',')}
        outputSample = ot.Sample(Y['dmax'][:,np.newaxis])
        outputSample.setDescription(['d_max'])
        
        
        # PCE
        ot.ResourceMap.SetAsUnsignedInteger("FunctionalChaosAlgorithm-MaximumTotalDegree", 15)
        algo = ot.FunctionalChaosAlgorithm(inputSample, outputSample)
        algo.run()
        result = algo.getResult()
        metamodel = result.getMetaModel()
        
        chaosSI = ot.FunctionalChaosSobolIndices(result)
        print(chaosSI.summary())
        
    
    #%% Pythia example
    if True:
        print("Tutorial 03 - 2D approximation with PC")

        # target function definition
        a = np.array([1, 2, 3])
        
        # analytical sobol coefficients
        sobol_dict = sobol_sc(a=a, dim=len(a))[0]
        sobol_coefficients = np.array(list(sobol_dict.values())).reshape(-1, 1)
        
        # setup pc surrogate
        params = [
            pt.parameter.Parameter(name=f"x_{j+1}", domain=[0, 1], distribution="uniform")
            for j in range(a.size)
        ]
        
        max_dim = 11
        # limit total polynomial degree of expansion terms to 10
        indices = pt.index.simplex_set(len(params), max_dim - 1)
        index_set = pt.index.IndexSet(indices)
        print("multiindex information:")
        print(f"    number of indices: {index_set.shape[0]}")
        print(f"    dimension: {index_set.shape[1]}")
        print(f"    number of sobol indices: {len(index_set.sobol_tuples)}")
        
        N = 10_000
        print(f"generate training data ({N})")
        s = pt.sampler.WLSTensorSampler(params, [max_dim - 1] * len(params))
        x_train = s.sample(N)
        w_train = s.weight(x_train)
        y_train = target_function(x_train)
        
        print("compute pc expansion")
        surrogate = pt.chaos.PolynomialChaos(params, index_set, x_train, w_train, y_train)
        
        # test PC approximation
        N = 1000
        print(f"generate test data ({N})")
        s_test = pt.sampler.ParameterSampler(params)
        x_test = s_test.sample(N)
        y_test = target_function(x_test)
        y_approx = surrogate.eval(x_test)
        
        error_L2 = np.sqrt(np.sum((y_test - y_approx) ** 2) / N)
        error_L2_rel = error_L2 / np.sqrt(np.sum((y_test) ** 2) / N)
        error_max = np.max(np.abs(y_test - y_approx))
        error_max_rel = np.max(np.abs(y_test - y_approx) / np.abs(y_test))
        
        print(f"test error L2 (abs/rel): {error_L2:4.2e} / {error_L2_rel:4.2e}")
        print(f"test error max (abs/rel): {error_max:4.2e} / {error_max_rel:4.2e}")
        
        # compare Sobol indices
        print("Comparison of Sobol indices")
        print(f" {'sobol_tuple':<12} {'exact':<8}  {'approx':<8}  {'abs error':<9}")
        print("-" * 44)
        for j, sdx in enumerate(sobol_dict.keys()):
            print(
                f" {str(sdx):<11} ",  # Sobol index subscripts
                f"{sobol_coefficients[j, 0]:<4.2e} ",
                f"{surrogate.sobol[j, 0]:<4.2e} ",
                f"{np.abs(sobol_coefficients[j, 0] - surrogate.sobol[j, 0]):<4.2e}",
            )
        
            
    
    #%% Openturns example
    if False:
        """https://openturns.github.io/openturns/latest/auto_meta_modeling/polynomial_chaos_metamodel/plot_functional_chaos.html#sphx-glr-auto-meta-modeling-polynomial-chaos-metamodel-plot-functional-chaos-py
        
        How to provide input and output Sample?
        
        """
        
        ot.RandomGenerator.SetSeed(0)
        dimension = 2
        input_names = ["x1", "x2"]
        formulas = ["cos(x1 + x2)", "(x2 + 1) * exp(x1)"]
        model = ot.SymbolicFunction(input_names, formulas)
        
        distribution = ot.Normal(dimension)
        samplesize = 80
        inputSample = distribution.getSample(samplesize)
        outputSample = model(inputSample)
        
        ot.ResourceMap.SetAsUnsignedInteger("FittingTest-LillieforsMaximumSamplingSize", 100)  # ??
        
        algo = ot.FunctionalChaosAlgorithm(inputSample, outputSample)
        algo.run()
        result = algo.getResult()
        metamodel = result.getMetaModel()



        x1index = 0
        x1value = 0.5
        x2min = -3.0
        x2max = 3.0
        outputIndex = 1
        metamodelParametric = ot.ParametricFunction(metamodel, [x1index], [x1value])
        graph = metamodelParametric.getMarginal(outputIndex).draw(x2min, x2max)
        graph.setLegends(["Metamodel"])
        modelParametric = ot.ParametricFunction(model, [x1index], [x1value])
        curve = modelParametric.getMarginal(outputIndex).draw(x2min, x2max).getDrawable(0)
        curve.setColor("red")
        curve.setLegend("Model")
        graph.add(curve)
        graph.setLegendPosition("bottomright")
        graph.setXTitle("X2")
        graph.setTitle("Metamodel Validation, output #%d" % (outputIndex))
        view = viewer.View(graph)
        
        n_valid = 100
        inputTest = distribution.getSample(n_valid)
        outputTest = model(inputTest)
        
        # Plot the corresponding validation graphics.
        val = ot.MetaModelValidation(inputTest, outputTest, metamodel)
        Q2 = val.computePredictivityFactor()
        graph = val.drawValidation()
        graph.setTitle("Metamodel validation Q2=" + str(Q2))
        view = viewer.View(graph)


        
        
        #---Compute and print Sobol’ indices---
        chaosSI = ot.FunctionalChaosSobolIndices(result)
        print(chaosSI.summary())
        
        
        #---Testing the sensitivity to the degree---
        ot.ResourceMap.GetAsUnsignedInteger("FunctionalChaosAlgorithm-MaximumTotalDegree")
        
        degrees = range(5, 12)
        q2 = ot.Sample(len(degrees), 2)
        for maximumDegree in degrees:
            ot.ResourceMap.SetAsUnsignedInteger(
                "FunctionalChaosAlgorithm-MaximumTotalDegree", maximumDegree
            )
            print("Maximum total degree =", maximumDegree)
            algo = ot.FunctionalChaosAlgorithm(inputSample, outputSample)
            algo.run()
            result = algo.getResult()
            metamodel = result.getMetaModel()
            for outputIndex in range(2):
                val = ot.MetaModelValidation(
                    inputTest, outputTest[:, outputIndex], metamodel.getMarginal(outputIndex)
                )
                q2[maximumDegree - degrees[0], outputIndex] = val.computePredictivityFactor()[0]
        
        graph = ot.Graph("Predictivity", "Total degree", "Q2", True)
        cloud = ot.Cloud([[d] for d in degrees], q2[:, 0])
        cloud.setLegend("Output #0")
        cloud.setPointStyle("bullet")
        graph.add(cloud)
        cloud = ot.Cloud([[d] for d in degrees], q2[:, 1])
        cloud.setLegend("Output #1")
        cloud.setColor("red")
        cloud.setPointStyle("bullet")
        graph.add(cloud)
        graph.setLegendPosition("topright")
        view = viewer.View(graph)
        plt.show()
        
        ot.ResourceMap.Reload() #Reset default settings

