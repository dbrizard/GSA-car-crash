#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:53:00 2024

@author: dbrizard
"""

import numpy as np
import matplotlib.pyplot as plt

from SALib.test_functions import Ishigami, Sobol_G
from SALib.sample import morris
from SALib.analyze import morris as momo

import uqtestfuns as uqtf  # requires python38 venv !!

from GSAutils import MorrisResults, GatherMorris

class MorrisOutput(MorrisResults):
    """The same class, except that it takes a SALib result as input rather
    than a text file
    
    """
    def __init__(self, resultDict, outname=''):
        """
        
        :param dict resultDict: ResultDict from SALib
        """
        self._readSAlibDict(resultDict)
        self.outname= outname
        self.sortbyMuStar()
        
    
    def _readSAlibDict(self, rd):
        """Read Morris analysis output dict and store it the same way it is done
        in :meth:`self._readSAlibOutput`
        
        
        :param dict rd: ResultDict from SALib with Morris analysis output
        """
        SI = (rd['mu_star'], rd['mu'], rd['mu_star_conf'], rd['sigma'])
        # mu_star, mu, m_star_conf, sigma
        self.SI = np.vstack(SI).T
        self.param = rd['names']


if __name__=="__main__":
    plt.close('all')


    #%% REPEAT MORRIS ANALYSIS
    if True:
        # Choose benchmark functino here
        funcname = 'BratleyB'        
        nTraj = 20  # number of trajectories for the Morris analysis
        nRep = 100  # number of repetitions of the Morris analyses
        
        # Define the model inputs and function
        if funcname=='Ishigami':
            func = uqtf.Ishigami()
            problem = {'num_vars': 3, 'names': ['x1', 'x2', 'x3'],
                       'bounds': [[-np.pi, np.pi]]*3, 'groups':None}  # required for Morris method
        else:
            if funcname=='Sobol_G*':
                dim = 8
                func = Sobol_G.evaluate  #a = np.array([0, 1, 4.5, 9, 99, 99, 99, 99])
            elif funcname=='Sobol_G':
                dim = 6
                func = uqtf.SobolG(spatial_dimension=dim)
            # elif funcname=='Friedman6D':
            #     dim = 6
            #     func = uqtf.Friedman6D  # apprears in the doc but not available
            elif funcname=='BratleyB':
                dim = 8
                func = uqtf.Bratley1992b(spatial_dimension=dim)
            problem = {'num_vars':dim, 'names':['x%i'%(ii+1) for ii in range(dim)],
                       'bounds':[[0, 1]]*dim, 'groups':None}

                
                        
        # Repeat the analysis
        MOlist = []
        for ii in range(nRep):    
            # Run Morris analysis
            Xm = morris.sample(problem, nTraj, num_levels=4)  # , grid_jump=2
            Ym = func(Xm)
            Sim = momo.analyze(problem, Xm, Ym, conf_level=0.95,
                               print_to_console=True, num_levels=4)  # , grid_jump=2)
            # Store results
            Mo = MorrisOutput(Sim)
            MOlist.append(Mo)
            
        
        GM = GatherMorris(MOlist, 'rep', [ii for ii in range(nRep)], funcname)
        GM.plot2D(conf=False)
        GM.plot(conf=False)
        GM.plotRanking()
