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
from matplotlib import ticker
# from matplotlib import pylab as plt
import openturns as ot  # ot.__version__ => '1.19.post1' before, now 1.24 installed via yum.
import openturns.experimental as otexp
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
    

class OpenTurnsPCESobol():
    """My take on PCE Sobol indices with Openturns
    
    """
    
    def __init__(self, basepath='LHS/9param/LHS-', ns=120, prob=9):
        """Set the problem, inputs and outputs
        
        :param str basepath: base path for the input and output files
        :param int ns: number of samples in the LHS DOE
        """
        mm, mpa = 'mm', 'MPa'
        car_output = ('dmax', 'fmax', 'IE', 'vfin')
        car_out_unit = ('mm', 'N', '', 'mm/s', 'Nmm')
        if prob==9:
            # problem with 9 uncertain parameters
            problem = {'names': ['tbumper', 'trailb', 'trailf', 'tgrill', 'thood',
                                 'ybumper', 'yrailf', 'yrailb', 'ybody'],
                       'units': [mm]*5 + [mpa]*4,
                       'num_vars': 9,
                       'bounds': [[2, 4], [1,3], [3,7], [0.5,1.5], [0.5, 1.5],
                                  [300, 500], [300, 500], [300, 500], [300, 500]],
                       }
            output = car_output
            out_units = car_out_unit
        elif prob==4:
            # toy car crash problem with 4 uncertain parameters
            problem = {'names': ['tbumper', 'trailb', 'trailf', 'yrailf'],
                       'units': [mm]*3 + [mpa],
                       'num_vars':4,
                       'bounds': [[2,4], [1,3], [3,7], [300,500]]}
            output = car_output
            out_units = car_out_unit
        elif prob==110:
            # Wood barriere problem (Nyobe PhD)
            lbound = [30555.56, 0.001, 19, 9, 4.4e-10, 5, 194, 96, 194, 243]
            ubound = [32694.44, 0.075, 21.5, 20, 7.1e-10, 20, 208, 106, 208, 259]
            
            problem = {'names':['VIV', 'AJM', 'ANI', 'TEE', 'DEN', 'TEM', 
                                'POX', 'POY', 'LIY', 'LIZ'],
                       'units':['mm/s', 't', '°', '%', 't/mm3', '°C']+['mm']*4,
                       'dyna_names': ['speed', 'addmass', 'beta', 'mois', 'rho', 'temp',
                                      'postx', 'posty', 'beamy', 'beamz'],  
                       'num_vars':10,
                       'bounds':[(lb, ub) for lb,ub in zip(lbound, ubound)]}
            output = ('ASI', 'THIV', 'Dm', 'Wm', 
                      'post1', 'post2', 'post3', 'post4')
            out_units = ['-', 'km/h'] + ['mm']*6
            
        self.problem = problem
        self.input = problem['names']
        self.output = output
        self.output_units = out_units
        
        # Import X and Y
        self.Y = {}
        if prob in (9, 4):
            self.X = np.loadtxt('%s%i_X.csv'%(basepath, ns), delimiter=',')
            for oo in self.output:
                self.Y[oo] = np.loadtxt('%s%i_Y_%s.csv'%(basepath, ns, oo), delimiter=',')
        elif prob==110:
            self.X = np.loadtxt('%s%i_X.csv'%(basepath, ns), delimiter=';', skiprows=1)
            Y = np.loadtxt('%s%i_Y.csv'%(basepath, ns), delimiter=';', skiprows=1)
            for ii, oo in enumerate(self.output):
                self.Y[oo] = Y[:,ii]

        # OpenTURNS variables and definitions 
        ot.ResourceMap.SetAsUnsignedInteger("FittingTest-LillieforsMaximumSamplingSize", 100)
        ot.ResourceMap.SetAsBool("FunctionalChaosValidation-ModelSelection", True)
        self.inputSample = ot.Sample(self.X)
        self.inputSample.setDescription(self.input)
        self.outputSample = {}
        for oo in self.output:
            self.outputSample[oo] = ot.Sample(self.Y[oo][:,np.newaxis])
            self.outputSample[oo].setDescription([oo])
            
        distlist = [ot.Uniform(aa, bb) for (aa, bb) in self.problem['bounds']]
        # distribution = ot.ComposedDistribution([ot.Uniform()]*len(self.input))
        self.distribution = ot.ComposedDistribution(distlist)
    
    
    def computeChaosSensitivity(self, strategy='cleaning', q=0.4):
        """Compute PCE metamodel and get Sobol indices for all outputs

        :param str strategy: adaptive strategy ('fixed' or 'cleaning')
        :param float q: q-quasi norm parameter. If not precised, q = 0.4. (see HyperbolicAnisotropicEnumerateFunction)
        """
        self.chaosSI = {}
        self.metamodel = {}
        self.validation = {}
        self.S1 = {}
        self.ST = {}
        for oo in self.output:
            print('='*20, oo, '='*20, '\n')
            self.S1[oo], self.ST[oo] = self._computeChaosSensitivity(self.inputSample, self.outputSample[oo], strategy=strategy, q=q, verbose=False)
            
    
    def _computeChaosSensitivity(self, inputSample, outputSample, 
                                 strategy='cleaning', q=0.4, verbose=False,
                                 validation=True):
        """Compute PCE metamodel and return Sobol indices for a specific output
        
        :param Sample inputSample: 
        :param Sample outputSample: 
        :param str strategy: adaptive strategy ('fixed' or 'cleaning')
        :param float q: q-quasi norm parameter. If not precised, q = 0.4. (see HyperbolicAnisotropicEnumerateFunction)
        """

        # https://openturns.github.io/openturns/latest/user_manual/response_surface/_generated/openturns.FixedStrategy.html#openturns.FixedStrategy
        polyColl = [0.0]*len(self.input)
        for i in range(self.distribution.getDimension()):
            polyColl[i] = ot.StandardDistributionPolynomialFactory(self.distribution.getMarginal(i))
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


        oo = outputSample.getDescription()[0]
        # PCE
        # ot.ResourceMap.SetAsUnsignedInteger("FunctionalChaosAlgorithm-MaximumTotalDegree", 15)
        algo = ot.FunctionalChaosAlgorithm(
            inputSample, outputSample, self.distribution, adaptiveStrategy)
        algo.run()
        result = algo.getResult()
        self.metamodel[oo] = result.getMetaModel()
        
        # SOBOL INDICES
        self.chaosSI[oo] = ot.FunctionalChaosSobolIndices(result)
        if verbose:
            print(self.chaosSI[oo].summary())
    
        S1 = [self.chaosSI[oo].getSobolIndex(ii) for ii in range(len(self.input))]
        ST = [self.chaosSI[oo].getSobolTotalIndex(ii) for ii in range(len(self.input))]
        self.S1[oo] = S1
        self.ST[oo] = ST
    
        # XXX Metamodel validation
        # val = ot.MetaModelValidation(X_test, Y_test, metamodel)
        # https://openturns.github.io/openturns/latest/auto_meta_modeling/polynomial_chaos_metamodel/plot_chaos_sobol_confidence.html#sphx-glr-auto-meta-modeling-polynomial-chaos-metamodel-plot-chaos-sobol-confidence-py
        
        # VALIDATION
        if validation:
            splitterLOO = ot.LeaveOneOutSplitter(len(inputSample))
            validation = otexp.FunctionalChaosValidation(result, splitterLOO)
            r2Score = validation.computeR2Score()
            print('R2 = ', r2Score[0])
            self.validation[oo] = validation
            
        return S1, ST


    def computeBootstrapChaosSobolIndices(self, bootstrap_size, pick=False, verbose=True):
        """
        Computes a bootstrap sample of first and total order indices from polynomial chaos.
        
        https://openturns.github.io/openturns/latest/auto_meta_modeling/polynomial_chaos_metamodel/plot_chaos_sobol_confidence.html#sphx-glr-auto-meta-modeling-polynomial-chaos-metamodel-plot-chaos-sobol-confidence-py

    
        :param interval bootstrap_size: The bootstrap sample size
        :param bool pick: discard indices when 0 or 1
        :param bool verbose: tell human where computer is
        """
        X = self.inputSample
        dim_input = X.getDimension()

        FO = {}
        TO = {}
        FOI = {}
        TOI = {}
        for oo in self.output:
            if verbose:
                print('*'*20, oo, '*'*20, '\n')
            Y = self.outputSample[oo]
            fo_sample = ot.Sample(0, dim_input)
            to_sample = ot.Sample(0, dim_input)
            unit_eps = ot.Interval([1e-9] * dim_input, [1 - 1e-9] * dim_input)
            for ii in range(bootstrap_size):
                if verbose and ii%20==0:
                    print('bootstrap %i/%i'%(ii, bootstrap_size))
                X_boot, Y_boot = multiBootstrap(X, Y)
                first_order, total_order = self._computeChaosSensitivity(X_boot, Y_boot, validation=False)
                if pick:
                    # do not add to sample if any first_order or total_order is 0.0
                    if unit_eps.contains(first_order) and unit_eps.contains(total_order):
                        fo_sample.add(first_order)
                        to_sample.add(total_order)
                else:
                    fo_sample.add(first_order)
                    to_sample.add(total_order)                    

            # compute confidence intervals
            fo_interval, to_interval = computeSobolIndicesConfidenceInterval(fo_sample, to_sample)
            # Store
            FO[oo] = fo_sample
            TO[oo] = to_sample
            FOI[oo] = fo_interval
            TOI[oo] = to_interval
        
        self.bootstrap = {'FO':FO, 'TO':TO, 'FOI':FOI, 'TOI':TOI}


    def plotS1STbootstrap(self, figname='', method='mine', 
                          xmargin=0.2, xoffset=0, xST1=0.1, labelsuffix=''):
        """
        
        :param str figname: prefix for figure name
        :param str method: plotting method ('OT' for OpenTURNS, otherwise mine)
        :param float xmargin: left and right margin around x axis
        :param float xoffset: additionnal x offset (for overlay purpose)
        :param float xST1: x offset between ST and S1
        """
        for oo in self.output:
            fo_sample = self.bootstrap['FO'][oo]
            to_sample = self.bootstrap['TO'][oo]
            fo_interval = self.bootstrap['FOI'][oo]
            to_interval = self.bootstrap['TOI'][oo]
            
            if method=='OT':
                graph = ot.SobolIndicesAlgorithm.DrawSobolIndices(
                    self.inputSample.getDescription(),
                    fo_sample.computeMean(),
                    to_sample.computeMean(),
                    fo_interval,
                    to_interval,
                )
                graph.setTitle(f"Sobol' indices: {oo}")
                fig = plt.figure('%s-%s'%(figname, oo))
                _ = viewer.View(graph, figure=fig)
            else:
                x = np.arange(len(self.input)) + 1
                plt.figure('%s-%s'%(figname, oo))
                # ST
                ST = to_sample.computeMean()
                STerr = np.vstack((ST-to_interval.getLowerBound(), to_interval.getUpperBound()-ST))
                plt.errorbar(x+xoffset, ST, yerr=STerr, label='ST'+labelsuffix,
                             marker='o', elinewidth=2, linestyle='')
                # S1
                S1 = fo_sample.computeMean()
                S1err = np.vstack((S1-fo_interval.getLowerBound(), fo_interval.getUpperBound()-S1))
                plt.errorbar(x+xST1+xoffset, S1, yerr=S1err, label='S1'+labelsuffix,
                             marker='x', elinewidth=2, linestyle='')
                # graph tuning                
                plt.xlim(xmin=x.min()-xmargin, xmax=x.max()+xmargin)
                plt.xticks(x, self.input, rotation=30)
                plt.ylim((0,1))
                plt.legend()
                plt.title(f"Sobol' indices: {oo}")
                # plt.box(False)
        
                
        
    def plotS1ST(self, figname='', color=None, label='', 
                 xmargin=0.3, xoffset=0.2, ylim=True):
        """Plot S1 and ST, for each output, on the same graph
        
        :param str figname: prefix for the name of the figures
        :param str color: color for the markers
        :param float xmargin: x axis margins
        :param float xoffset: horizontal offset of the points (in [0, 1] interval)
        :param bool ylim: set ylim to [0,1]
        """
        for oo in self.output:
            plt.figure('%s-%s'%(figname, oo))
            x = np.arange(0, len(self.input)) + xoffset
            plt.plot(x, self.ST[oo], '+', label='ST_%s'%label, ms=14, color=color)
            plt.plot(x, self.S1[oo], 'x', label='S1_%s'%label, ms=10, color=color)
            plt.xlim(xmin=-xmargin, xmax=len(self.input)-1+xmargin)
            plt.xticks(ticks=range(len(self.input)), labels=self.input , rotation=45)
            plt.title(oo)
            plt.legend()
            if ylim:
                plt.ylim([0,1]) 
    
    
    def plotRanking(self, figname=''):
        """
        
        :param str figname: base name for the figure (S1/ST and output name are added)
        """
        for oo in self.output:
            iS1 = np.argsort(self.S1[oo]) + 1
            iST = np.argsort(self.ST[oo]) + 1
            S1 = np.sort(self.S1[oo])
            ST = np.sort(self.ST[oo])

            fn = '%s_%s_'%(figname, oo)
            plotSobolRanking(iS1[::-1,np.newaxis], figname=fn+'S1',
                             yticks=S1[::-1], xlabel='Sobol S1')
            plotSobolRanking(iST[::-1,np.newaxis], figname=fn+'ST',
                             yticks=ST[::-1], xlabel='Sobol ST')


    def getR2(self,):
        """Get R2 validation values for each output as a dict
        
        """
        R2 = {}
        for oo in self.output:
            R2[oo] = self.validation[oo].computeR2Score()[0]
        return R2
    
    def plotValidation(self,):
        """
        
        """
        R2 = self.getR2()
        for oo, uu in zip(self.output, self.output_units):
            graph = self.validation[oo].drawValidation()
            graph.setTitle('%s /%s [$R^2=%.3f$]'%(oo, uu, R2[oo]))
            # graph.setXTitle('test')  # XXX AttributeError ???
            viewer.View(graph, figure_kw={'num':'validation-%s'%oo})


def plotSobolRanking(matrix, figname=None, xlabel=None, yticks=None):
    """
    
    """
    nparam = matrix.shape[0]
    nrepet = matrix.shape[1]
    
    fig, ax = plt.subplots(num=figname, figsize=(nrepet*1.5, nparam/3))
    ax.matshow(matrix, cmap='viridis_r')
    
    # plt.title(title)
    plt.xlabel(xlabel, ha='left')
    plt.ylabel('ranking')
    if yticks is not None:
        plt.xticks([])
        # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
        yticks = ['%.3f'%tt for tt in yticks]
        plt.yticks(range(nparam), yticks)
        ax.yaxis.tick_right()
    else:
        plt.xticks([])
        plt.yticks([])
    plt.box(False)
    
    for ii in range(nrepet):
        for jj in range(nparam):
            ax.text(ii, jj, str(matrix[jj,ii]), va='center', ha='center')


#%% The following functions were taken from:
# https://openturns.github.io/openturns/latest/auto_meta_modeling/polynomial_chaos_metamodel/plot_chaos_sobol_confidence.html#sphx-glr-auto-meta-modeling-polynomial-chaos-metamodel-plot-chaos-sobol-confidence-py

def multiBootstrap(*data):
    """
    Bootstrap multiple samples at once.

    Parameters
    ----------
    data : sequence of Sample
        Multiple samples to bootstrap.

    Returns
    -------
    data_boot : sequence of Sample
        The bootstrap samples.
    """
    assert len(data) > 0, "empty list"
    size = data[0].getSize()
    selection = ot.BootstrapExperiment.GenerateSelection(size, size)
    return [Z[selection] for Z in data]


def computeSobolIndicesConfidenceInterval(fo_sample, to_sample, alpha=0.95):
    """
    From a sample of first or total order indices,
    compute a bilateral confidence interval of level alpha.

    Estimates the distribution of the first and total order Sobol' indices
    from a bootstrap estimation.
    Then computes a bilateral confidence interval for each marginal.

    Parameters
    ----------
    fo_sample: ot.Sample(n, dim_input)
        The first order indices
    to_sample: ot.Sample(n, dim_input)
        The total order indices
    alpha : float
        The confidence level

    Returns
    -------
    fo_interval : ot.Interval
        The confidence interval of first order Sobol' indices
    to_interval : ot.Interval
        The confidence interval of total order Sobol' indices
    """
    dim_input = fo_sample.getDimension()
    fo_lb = [0] * dim_input
    fo_ub = [0] * dim_input
    to_lb = [0] * dim_input
    to_ub = [0] * dim_input
    for i in range(dim_input):
        fo_i = fo_sample[:, i]
        to_i = to_sample[:, i]
        beta = (1.0 - alpha) / 2
        fo_lb[i] = fo_i.computeQuantile(beta)[0]
        fo_ub[i] = fo_i.computeQuantile(1.0 - beta)[0]
        to_lb[i] = to_i.computeQuantile(beta)[0]
        to_ub[i] = to_i.computeQuantile(1.0 - beta)[0]

    # Create intervals
    fo_interval = ot.Interval(fo_lb, fo_ub)
    to_interval = ot.Interval(to_lb, to_ub)
    return fo_interval, to_interval



if __name__=='__main__':
    plt.close('all')
    
    #%% Plot Eric's results
    if False:
        Eric = EricPCESobol()
        Eric.plotS1ST(figname='S1ST')


    #%% Openturns on LS-DYNA car simulation data: 9 uncertain parameters
    if False:
        bs = 500
        if True:
            OTS50 = OpenTurnsPCESobol(ns=50)
            OTS50.computeChaosSensitivity()
            OTS50.plotS1ST(figname='S1ST', color='C0', label='LHS-50')
            # OTS50.plotRanking(figname='sobol50')
            OTS50.computeBootstrapChaosSobolIndices(bs)  # influence de N ???
            OTS50.plotS1STbootstrap(figname='STS1-50-bs%i'%bs)
            
        if True:
            OTS120 = OpenTurnsPCESobol(ns=120)
            OTS120.computeChaosSensitivity()
            OTS120.plotS1ST(figname='S1ST', color='C1', label='LHS-120')
            # OTS120.plotRanking(figname='sobol120')
                
            OTS120.computeBootstrapChaosSobolIndices(bs)  # influence de N ???
            OTS120.plotS1STbootstrap(figname='STS1-120-bs%i'%bs)
    
        
        if True:
            OTS330 = OpenTurnsPCESobol(ns=330)
            OTS330.computeChaosSensitivity()
            OTS330.plotS1ST(figname='S1ST', color='C2', label='LHS-330')
            # OTS330.plotRanking(figname='sobol330')
            # OTS330.validation['vfin'].drawValidation()
             
            OTS330.computeBootstrapChaosSobolIndices(bs)  # influence de N ???
            OTS330.plotS1STbootstrap(figname='STS1-330-bs%i'%bs)
            # TODO: metamodel quality
        
        if True:
            OTS50.plotS1STbootstrap( figname='STS1-50-120-330_bs%i'%bs, labelsuffix='-50', xST1=0.07)
            OTS120.plotS1STbootstrap(figname='STS1-50-120-330_bs%i'%bs, labelsuffix='-120', xST1=0.07, xoffset=0.2)
            OTS330.plotS1STbootstrap(figname='STS1-50-120-330_bs%i'%bs, labelsuffix='-330', xST1=0.07, xoffset=0.4)
            
        
    #%% Openturns on LS-DYNA car simulation data: 4 uncertain parameters
    if False:
        bs = 500
        if False:
            S100 = OpenTurnsPCESobol(basepath='LHS/4param/LHS-', ns=100, prob=4)
            S100.computeChaosSensitivity()
            S100.plotS1ST(figname='S1ST', label='LHS-100')
            S100.computeBootstrapChaosSobolIndices(bs)
            S100.plotS1STbootstrap(figname='S1ST-100_bs%i'%bs, labelsuffix='-100')
        
        if True:
            S200 = OpenTurnsPCESobol(basepath='LHS/4param/LHS-', ns=200, prob=4)
            S200.computeChaosSensitivity()
            S200.plotS1ST(figname='S1ST', label='LHS-200')
            S200.computeBootstrapChaosSobolIndices(bs)
            S200.plotS1STbootstrap(figname='S1ST-200_bs%i'%bs, labelsuffix='-200')
            S200.plotRanking(figname='sobol200')
        
        if False:
            S100.plotS1STbootstrap(figname='S1ST-100-200_bs%i'%bs, labelsuffix='-100', xST1=0.07)
            S200.plotS1STbootstrap(figname='S1ST-100-200_bs%i'%bs, labelsuffix='-200', xST1=0.07, xoffset=0.2, xmargin=0.3)

    #%% OpenTURNS Nyobe VRS
    if True:
        VRS60 = OpenTurnsPCESobol(basepath='VRS/Sobol-', ns=60, prob=110)
        VRS60.computeChaosSensitivity()
        VRS60.plotS1ST(figname='S1ST')
        VRS60.plotValidation()
        
        # bs= 500
        # VRS60.computeBootstrapChaosSobolIndices(bs)
        # VRS60.plotS1STbootstrap(figname='S1ST_bs%i'%bs)
        
        

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

