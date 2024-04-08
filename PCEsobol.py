#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Trying different modules for PCE approximation.

- [x] Openturns runnings
- [ ] Pythia-uq 
  - 2.0.0 installed (no documentation online)
  - latest is 4.0.3 (runs with python 3.9)
- [ ] pygpc 0.2.7.5
  - I cannot provide simulation output myself
- [ ] Chaospy 4.2.7 installed. 
  - Usage does not seem straightforward
- [ ] Uncertainpy
- [ ] UQpy 4.1.5 
  - https://uqpyproject.readthedocs.io/en/latest/auto_examples/surrogates/pce/index.html



Created on Fri Apr  5 11:03:36 2024

@author: dbrizard
"""

import os
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import pylab as plt
import openturns as ot
import openturns.viewer as viewer
ot.Log.Show(ot.Log.NONE)


# from UQpy.surrogates.polynomial_chaos.PolynomialChaosExpansion import PolynomialChaosExpansion
# from UQpy.sensitivity.PceSensitivity import PceSensitivity


class EricPCESobol():
    """Import and plot Sobol indices computed by Eric
    
    """
    
    def __init__(self):
        """Load and store data
        
        """
        folder = '/home/dbrizard/Calcul/25_car/GSA/car_v223_right-impact_v30/Eric/EricPCEsobol'
        self.S1 = {}
        self.S1_mean = {}
        self.ST = {}
        self.ST_mean = {}
        self.output = ('dmax', 'fmax', 'IE', 'vfin')
        self.input = ['tbumper', 'trailb', 'trailf', 'tgrill', 'thood', 'ybumper', 'yrailf', 'yrailb', 'ybody']
        
        for oo in self.output:
            self.S1[oo] = np.genfromtxt(os.path.join(folder,'SI-Y-%s.txt'%oo), delimiter='  ')  # loadtxt not working
            self.ST[oo] = np.genfromtxt(os.path.join(folder,'SI-tot-Y-%s.txt'%oo), delimiter='  ')  # loadtxt not working
            self.S1_mean[oo] = self.S1[oo].mean(axis=0)
            self.ST_mean[oo] = self.ST[oo].mean(axis=0)
        
    def plotS1ST(self, figname=None, xmargin=0.2, ylim=True):
        """
        
        :param float xmargin: 
        :param bool ylim: set ylim to [0,1]
        """
        for oo in self.output:
            if figname=='output':
                plt.figure('S1ST-%s'%oo)
            else:
                plt.figure(figname)
            plt.plot(self.S1[oo].T, '+k')
            plt.plot(self.ST[oo].T, '+r')
            plt.plot(self.S1_mean[oo], '+k', label='S1_mean', ms=15)
            plt.plot(self.ST_mean[oo], '+r', label='ST_mean', ms=15)
            plt.xlim(xmin=-xmargin, xmax=len(self.input)-1+xmargin)
            plt.xticks(ticks=range(len(self.input)), labels=self.input, rotation=45)
            plt.title(oo)
            plt.legend()
            if ylim:
                plt.ylim([0,1])
    


if __name__=='__main__':
    plt.close('all')
    
    #%% Plot Eric's results
    if True:
        Eric = EricPCESobol()
        Eric.plotS1ST(figname='output')


    #%% Openturns on LS-DYNA car simulation data
    if True:
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
                
        
        # TODO: compare with Eric
        # TODO: metamodel quality
            
    
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

