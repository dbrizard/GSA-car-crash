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
        
    def plotS1ST(self, figname='', xmargin=0.3, ylim=True):
        """Plot S1 and ST, for each output, on the same graph
        
        :param str figname: prefix for the name of the figures 
        :param float xmargin: x axis margins
        :param bool ylim: set ylim to [0,1]
        """
        for oo in self.output:
            plt.figure('%s-%s'%(figname, oo))
            plt.plot(self.ST[oo].T, '+r')
            plt.plot(self.S1[oo].T, '+k')
            plt.plot(self.ST_mean[oo], '+r', label='ST_mean', ms=15)
            plt.plot(self.S1_mean[oo], '+k', label='S1_mean', ms=15)
            plt.xlim(xmin=-xmargin, xmax=len(self.input)-1+xmargin)
            plt.xticks(ticks=range(len(self.input)), labels=self.input, rotation=45)
            plt.title(oo)
            plt.legend()
            if ylim:
                plt.ylim([0,1])
    

class DenisPCESobol():
    """My take on PCE Sobol indices with Openturns
    
    """
    
    def __init__(self,strategy='cleaning', q=0.4):
        """Set the problem, inputs, outputs and compute PCE metamodel
        
        
        :param str strategy: adaptive strategy ('fixed' or 'cleaning')
        :param float q: q-quasi norm parameter. If not precised, q = 0.4. (see HyperbolicAnisotropicEnumerateFunction)
        """
        mm, mpa = 'mm', 'MPa'
        problem = {'names': ['tbumper', 'trailb', 'trailf', 'tgrill', 'thood',
                             'ybumper', 'yrailf', 'yrailb', 'ybody'],
                   'units': [mm, mm, mm, mm, mm,
                             mpa, mpa, mpa, mpa],
                   'num_vars': 9,
                   'bounds': [[2, 4], [1,3], [3,7], [0.5,1.5], [0.5, 1.5],
                              [300, 500], [300, 500], [300, 500], [300, 500]],
                   }
        
        self.problem = problem
        self.input = problem['names']
        self.output = ('dmax', 'fmax', 'IE', 'vfin')
        
        ot.ResourceMap.SetAsUnsignedInteger("FittingTest-LillieforsMaximumSamplingSize", 100)

        # Import X and Y
        self.X = np.loadtxt('LHS/LHS-120_X.csv', delimiter=',')
        inputSample = ot.Sample(self.X)
        inputSample.setDescription(self.input)

        distlist = [ot.Uniform(aa, bb) for (aa, bb) in problem['bounds']]
        # distribution = ot.ComposedDistribution([ot.Uniform()]*len(self.input))
        distribution = ot.ComposedDistribution(distlist)


        # https://openturns.github.io/openturns/latest/user_manual/response_surface/_generated/openturns.FixedStrategy.html#openturns.FixedStrategy
        polyColl = [0.0]*len(self.input)
        for i in range(distribution.getDimension()):
            polyColl[i] = ot.StandardDistributionPolynomialFactory(distribution.getMarginal(i))
        # enumerateFunction = ot.LinearEnumerateFunction(len(self.input))
        enumerateFunction = ot.HyperbolicAnisotropicEnumerateFunction(len(self.input), q)
        productBasis = ot.OrthogonalProductPolynomialFactory(polyColl, enumerateFunction)

        if strategy=='fixed':
            # Number of terms of the basis.
            degree = 2
            indexMax = enumerateFunction.getStrataCumulatedCardinal(degree)
            # XXX attention à l'overfitting !!
            print('indexMax', indexMax)
            adaptiveStrategy = ot.FixedStrategy(productBasis, indexMax)  # https://openturns.github.io/openturns/latest/user_manual/response_surface/_generated/openturns.FixedStrategy.html#openturns.FixedStrategy
        elif strategy=='cleaning':
            # Maximum index that can be used by the EnumerateFunction to 
            # determine the last term of the basis.
            maximumDimension = 200
            # Parameter that characterizes the cleaning strategy. 
            # It represents the number of efficient coefficients of the basis. 
            # Its default value is set to 20.
            maximumSize = 20
            # Parameter used as a threshold for selecting the efficient coefficients 
            # of the basis. The real threshold represents the multiplication of 
            # the significanceFactor with the maximum magnitude of the current 
            # determined coefficients. Its default value is equal to 1e^{-4}.
            significanceFactor = 1e-4
            adaptiveStrategy = ot.CleaningStrategy(productBasis, maximumDimension, maximumSize, significanceFactor)

        self.Y = {}
        self.chaosSI = {}
        self.metamodel = {}
        self.S1 = {}
        self.ST = {}
        for oo in self.output:
            print('='*20, oo, '='*20, '\n')
            self.Y[oo] = np.loadtxt('LHS/LHS-120_Y_%s.csv'%oo, delimiter=',')
            
            outputSample = ot.Sample(self.Y[oo][:,np.newaxis])
            outputSample.setDescription([oo])
            
            # PCE
            # ot.ResourceMap.SetAsUnsignedInteger("FunctionalChaosAlgorithm-MaximumTotalDegree", 15)
            algo = ot.FunctionalChaosAlgorithm(inputSample, outputSample, distribution, adaptiveStrategy)
            algo.run()
            result = algo.getResult()
            self.metamodel[oo] = result.getMetaModel()
            
            self.chaosSI[oo] = ot.FunctionalChaosSobolIndices(result)
            print(self.chaosSI[oo].summary())
        
            self.S1[oo] = [self.chaosSI[oo].getSobolIndex(ii) for ii in range(len(self.input))]
            self.ST[oo] = [self.chaosSI[oo].getSobolTotalIndex(ii) for ii in range(len(self.input))]
    
    
        # XXX Metamodel validation
        # val = ot.MetaModelValidation(X_test, Y_test, metamodel)
        # https://openturns.github.io/openturns/latest/auto_meta_modeling/polynomial_chaos_metamodel/plot_chaos_sobol_confidence.html#sphx-glr-auto-meta-modeling-polynomial-chaos-metamodel-plot-chaos-sobol-confidence-py
    
    def plotS1ST(self, figname='', xmargin=0.3, xoffset=0.2, ylim=True):
        """Plot S1 and ST, for each output, on the same graph
        
        :param str figname: prefix for the name of the figures 
        :param float xmargin: x axis margins
        :param bool ylim: set ylim to [0,1]
        """
        for oo in self.output:
            plt.figure('%s-%s'%(figname, oo))
            x = np.arange(0, len(self.input)) + xoffset
            plt.plot(x, self.ST[oo], '+r', label='ST_Openturns', ms=10)
            plt.plot(x, self.S1[oo], '+k', label='S1_Openturns', ms=10)
            plt.xlim(xmin=-xmargin, xmax=len(self.input)-1+xmargin)
            plt.xticks(ticks=range(len(self.input)), labels=self.input , rotation=45)
            plt.title(oo)
            plt.legend()
            if ylim:
                plt.ylim([0,1]) 
                
        

if __name__=='__main__':
    plt.close('all')
    
    #%% Plot Eric's results
    if True:
        Eric = EricPCESobol()
        Eric.plotS1ST(figname='S1ST')


    #%% Openturns on LS-DYNA car simulation data
    if True:
        Denis = DenisPCESobol()
        Denis.plotS1ST(figname='S1ST')
                
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

