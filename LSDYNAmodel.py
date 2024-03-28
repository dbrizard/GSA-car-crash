#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 11:14:04 2023

@author: dbrizard
"""

import re
import os
import shutil
import warnings
import subprocess as sbp
from timeit import default_timer as timer

import numpy as np
import matplotlib.pyplot as plt

from dynareadout import Binout


from SALib.sample import saltelli, morris, fast_sampler, latin, ff
from SALib.analyze import fast, rbd_fast, sobol, delta
from SALib.analyze import morris as morrisa
from SALib.analyze import ff as ffa


precision = 'single'
# Set LS-DYNA executables
# Do not forget to set the LSTC_FILE in the executable folder to get the license 
lsddir = '/home/dbrizard/Logiciels/lsdyna'
if precision=='single':
    lsdslvr = 'ls-dyna_smp_s_r1010_x64_redhat5_ifort160'
elif precision=='double':
    lsdslvr = 'ls-dyna_smp_d_r1010_x64_redhat5_ifort160'
SOLVER = os.path.join(lsddir, lsdslvr)
LSPP = "/home/dbrizard/Logiciels/lsprepost4.8_common/./lspp48"  # ls prepost executable

    
class LSDYNAmodel:
    """A class to handle a LS-DYNA model for parametric study and GSA.
    
    Run simulations, change parameters, fetch results. 
    
    """
    
    def __init__(self, kfile, param):
        """Define parameters

        :param str kfile: main .k file full path
        :param dict param: {'param':value} dictionnary
        
        :param list pname: full name of parameters
        :param list pvar: symbol for parameters
        :param list punit: units of parameters
        :param list pvalue: nominal values of parameters
        """
        # self.params = {'name':pname, 'var':pvar, 'unit':punit, 'value':pvalue}
        self.param = param
        abspath = os.path.abspath(kfile)
        self.kfile = {'basename':os.path.basename(abspath),
                      'dirname':os.path.dirname(abspath),
                      'abspath':abspath}
        
        self._checkMainK()  # check structure of main file and get include file
        # Copy given main file because it will be modified
        temp = os.path.splitext(abspath)
        paramfile = temp[0]+'_param'+temp[1]
        shutil.copy2(abspath, paramfile)
        self.kfile['paramfile'] = paramfile  # work on param file not to erase initial values in original file
        
        self._setPartNames()
        self.GSA = {}  # will be used to store GSA related things
    
    
    def _checkMainK(self):
        """Check main .k file (keywords) and get *INCLUDE file name
        
        TODO: get parameters names
        
        """
        with open(self.kfile['abspath'], 'r') as f:
            contents = f.read()
            keywords = re.findall('\*[A-Z_0-9]+', contents)
            include = re.findall('[A-Za-z0-9-_/\.]+\.k', contents)  # not sure this grabs all the possible file names
        
        # Check keywords
        KWlist = ['*KEYWORD', '*PARAMETER', '*INCLUDE', '*END']
        if keywords==KWlist:
            ok = True
        else:
            ok = False
            warnings.warn('Check keywords in %s'%self.kfile['abspath'])
            print('Required:', KWlist)
            print('Contains:', keywords)

        # Store include
        if len(include)==1:
            self.kfile['include'] = include[0]
        else:
            print('Problem with *INCLUDE detection', include)
    
    
    def _setPartNames(self):
        """Set list of part names. 
        
        :param list plist: list of part names, in order of increasing ID
        """
        self.partList = ('part1', 'part20', 'shell', 'ground')
    
    
    def overrideParam(self, param, write=True):
        """Rewrite main .k file with overridden parameters
        
        WARNING: parameter names are not checked
        
        :param dict param: {'param':value} dictionnary
        :param bool write: if False, no file written (for debug)
        """
        maink = '*KEYWORD\n*PARAMETER\n'
        # how many parameters are allowed on one line? 4??
        maink += ''.join(['r%s,%g,'%(kk,param[kk]) for kk in param])[:-1]
        maink += '\n*INCLUDE\n%s\n*END\n'%self.kfile['include']
        
        kmain = []
        kmain.append('*KEYWORD\n')
        kmain.append('*PARAMETER\n')
        plist = list(param.items())
        for ii in range(0, len(plist), 4):  #only 4 parameters are allowed per line
            chunk = plist[ii:ii+4]
            kmain.append(''.join(['r%s,%g,'%(pp, vv) for (pp,vv) in chunk])[:-1]+'\n')
        kmain.append('*INCLUDE\n')
        kmain.append(self.kfile['include']+'\n')
        kmain.append('*END\n')
        
        if write:
            with open(self.kfile['paramfile'], 'w') as f:
                f.writelines(kmain)  # coud also use writelines(L) with list L...
        else:
            print(kmain)
    
    
    def run(self, compute=False, ncpu=4, memory='20m', nodump=True, verbose=False):
        """Run LS-DYNA model in its folder
        
        """
        #---Build command---
        preview = '%s i=%s ncpu=%i memory=%s S=intfor'%(SOLVER, self.kfile['paramfile'], ncpu, memory)
        if nodump:
            preview += ' d=nodump'
        if verbose:
            print('PREVIEW: %s'%preview)
        
        basedir = os.getcwd()  # remember where we were
        os.chdir(self.kfile['dirname'])
        
        #---RUN LS-DYNA---
        if compute:
            with open('lsdynalog.txt','w') as f:
                sbp.call(preview, shell=True, stdout=f)         
        os.chdir(basedir)  # go back to basedir
        
        
    
    def fetchGLSTAT(self, comp=None):
        """Fetch energies in GLSTAT section of binout file
        
        :param list comp: list of components to get
        """
        binout = Binout(os.path.join(self.kfile['dirname'], 'binout'))
        components = {'TE':'total_energy', 'KE':'kinetic_energy', 
                      'HE':'hourglass_energy', 'IE':'internal_energy', 't':'time',
                      'GXV':'global_x_velocity'}
        
        if comp is None:
            comp = components
        
        GL = {}
        for ee in comp:
            try:
                GL[ee] = np.array(binout.read('glstat/%s'%components[ee]))
            except RuntimeError:
                print('"%s" not found in GLSTAT'%components[ee])
        
        self.GLSTAT = GL
        self.glstat_comp = components
        
        
    def fetchMATSUM(self, comp=None):
        """Fetch data in MATSUM section of binout file
        
        :param list comp: list of components to get
        """
        binout = Binout(os.path.join(self.kfile['dirname'], 'binout'))
        components = {'IE':'internal_energy', 'KE':'kinetic_energy', 
                      'HE':'hourglass_energy', 't':'time',
                      'xrbv':'x_rbvelocity'}
        
        if comp is None:
            comp = components
        
        MS = {}
        for ee in comp:
            try:
                MS[ee] = np.array(binout.read('matsum/%s'%components[ee]))
            except RuntimeError:
                print('"%s" not found in MATSUM'%components[ee])
        
        self.MATSUM = MS  
        self.matsum_comp = components
        
    
    def plotGLSTAT(self):
        """
        
        """
        plt.figure()
        for ee in self.GLSTAT:
            if not ee=='t':
                plt.plot(self.GLSTAT['t'], self.GLSTAT[ee], label=ee)
        plt.legend()
        plt.xlabel('time')
        plt.title('GLSTAT')
        
    
    def plotMATSUM(self, comp):
        """
        
        :param str comp: component of the MATSUM database to plot
        """
        plt.figure()
        plt.plot(self.MATSUM['t'], self.MATSUM[comp])
        plt.legend(self.partList)
        plt.xlabel('time')
        plt.ylabel(self.matsum_comp[comp])
        plt.title('MATSUM')
        
    
    def fetchSimulationResults(self):
        """Fetch simulation results (eg. in binout file)
        
        """
        #---Fetch results---
        time = np.linspace(0, 3*np.pi)
        signal = np.sin(time)
        maxs = max(signal)
        
        #---Store time series outputs---
        TS = VariableTimeSeries()
        TS.addVariable('signal', 's', time, signal, 's', 'm')
        self.TimeSeries = TS

        #---Store scalar outputs---
        VS = VariableScalar()
        VS.addVariable('max(S)', 'max(S)', maxs, 'm')
        self.ScalarOutput = VS
        


class VariableScalar:
    """A class to handle a list of scalar variables
    
    """
    
    def __init__(self):
        """
        
        """
        self.name = []
        self.var  = []
        self.value= []
        self.unit = []    
    
    
    def addVariable(self, name, var, value, unit):
        """Add a scalar variable to the list
        
        :param str name: full name of the variable
        :param str var: variable symbol
        :param float value: value of the variable
        :param str unit: unit of the variable
        """
        self.name.append(name)
        self.var.append(var)
        self.value.append(value)
        self.unit.append(unit)
    
    
    def __str__(self):
        """
        
        """
        sttr = '%i scalar variables list:\n'%len(self.name)
        for name, var, value, unit in zip(self.name, self.var, self.value, self.unit):
            sttr += '* %s: %s = %g [%s] \n'%(name, var, value, unit)
        return sttr
        


class VariableTimeSeries:
    """A class to handle time series variables
    """
    
    def __init__(self):
        """
        
        """
        self.name = []
        self.var = []
        self.time = []
        self.signal = []
        self.time_unit = []
        self.signal_unit = []
        
    
    def addVariable(self, name, var, time, signal, time_unit, signal_unit):
        """
        
        :param str name: full name of the variable
        :param str var: variable symbol
        :param array time: time vector
        :param array signal: time series vector
        :param str time_unit: unit for time
        :apram str signal_unit: unit for signal
        """
        self.name.append(name)
        self.var.append(var)
        self.time.append(time)
        self.signal.append(signal)
        self.time_unit.append(time_unit)
        self.signal_unit.append(signal_unit)
    
    
    def _getIndex(self, var):
        """
        
        :param str var: name of the variable to get
        """
        return self.var.index(var)
    
    
    def getVariable(self, var):
        """Return a dictionnary representing one time series
        
        :param str var: name of the variable to get
        """
        ind = self._getIndex(var)
        vdict = {'name':self.name[ind],
                 'var':self.var[ind],
                 'time':self.time[ind],
                 'signal':self.signal[ind],
                 'time_unit':self.time_unit[ind],
                 'signal_unit':self.signal_unit[ind]}
        return vdict
    
    
    def plot(self, var=None, figname=None, xylabels=True):
        """
        
        """
        ind =self._getIndex(var=var)
        plt.figure(figname)
        plt.plot(self.time[ind], self.signal[ind], '-')
        if xylabels:
            plt.xlabel('time [%s]'%self.time_unit[ind])
            plt.ylabel('%s [%s]'%(self.name[ind], self.signal_unit[ind]))
    
    
    def plotAll(self, figname=None):
        """
        
        """
        for vv in self.var:
            self.plot(var=vv)
        
        
    def __str__(self):
        """
        
        """
        sttr = '%i time series list:\n'%len(self.name)
        for nn, vv, tt, ss, tu, su in zip(self.name, self.var, self.time, self.signal, 
                                          self.time_unit, self.signal_unit):
            sttr += '* %s: %s (%i points) [%s], wrt time [%s]'%(nn, vv, len(tt), su, tu)
        return sttr
    
    
    def __repr__(self):
        return "%s(%r)"%(self.__class__, self.__dict__)  # not sure this is useful
            


class CAR6model(LSDYNAmodel):
    """Handle car6.k model

    """
    def _setPartNames(self):
        """
        
        """
        self.partList = ('Body', 'Bumper', 'Hood', 'Grill', 'Longeron (front)',
                         'Longeron (rear)', 'Roof', 'Post', 'Ground')
    
    
    def fetchSimulationResults(self):
        """
        
        """
        binout = Binout(os.path.join(self.kfile['dirname'], 'binout'))
        
        time = np.array(binout.read('nodout/time'))        
        xdispl = np.array(binout.read('nodout/x_displacement'))
        self.NODOUT = {'t':time, 'x_displacement':xdispl}
        
        rwtime = np.array(binout.read('rwforc/forces/time'))  
        xf = np.array(binout.read('rwforc/forces/x_force'))
        self.RWFORC = {'t':rwtime, 'x_force':xf}
        
        # # Write RW forces in csv file with cfile run
        # basedir = os.getcwd()  # remember where we were
        # os.chdir(self.kfile['dirname'])
        # lsppcfile = '%s -nographics lspost_getRWforces.cfile'%LSPP
        # with open('lspp.log','w') as f:
        #     sbp.call(lsppcfile, shell=True, stdout=f)
        # os.chdir(basedir)  # go back to working directory
        # # Now import csv file
        # # np.loadtxt() not working        
        # forces = np.genfromtxt(os.path.join(self.kfile['dirname'], 'rwforces.csv'), 
        #                        delimiter=',', skip_header=2, usecols=(0,1,2,3,4))
        # # Time,normal_force @ 2,x_force @ 2,y_force @ 2,z_force @ 2,
        
        # XXX works but takes too much time...
        
        #---Store time series 
        TS = VariableTimeSeries()
        TS.addVariable('X displacement node #167', 'xdisp167', time, xdispl[:,0], 's', 'mm')
        self.TimeSeries = TS
    
    
    def getGSAoutput(self):
        """Get the values of the scalar outputs for the GSA analysis        
        
        """
        out = {'fmax':max(self.RWFORC['x_force'][:,1]),  # force max sur poteau
               # 'dmax':self.NODOUT['x_displacement'][-1, 1] - self.NODOUT['x_displacement'][-1, 3],  # enfoncement maxi
               'dmax':self.NODOUT['x_displacement'][-1, 4],
               # 'vfin':self.GLSTAT['GXV'][-1][:],  # vitesse finale
               'vfin':self.MATSUM['xrbv'][-1, 0], # vitesse finale (part 1 = body)
               'IE':self.GLSTAT['IE'][-1][:]  # internal energy
               }
        # ajouter énergie interne
        return out
        
        
    
    def runGSA(self, N=10, meth='morris'):
        """
        
        
        :param int N: number of trajectories ('morris')
        :param str meth: GSA method ('morris')
        """
        # Define parameters and bounds
        if self.kfile['basename']=='main_v222.k':
            problem = {'names': ['tbumper', 'troof', 'trailb', 'trailf', 'tgrill', 'thood'],
                       'num_vars': 6,
                       'bounds': [[2, 4], [1, 3], [1,3], [3,7], [0.5,1.5], [0.5, 1.5]],
                       # 'groups': ['G1', 'G2', 'G1'], # Sobol and Morris only. See Advanced Examples.
                       # 'dists': ['unif', 'lognorm', 'triang']
                       }
        elif self.kfile['basename']=='main_v223.k':
            problem = {'names': ['tbumper', 'trailb', 'trailf', 'tgrill', 'thood',
                                 'ybumper', 'yrailf', 'yrailb', 'ybody'],
                       'num_vars': 9,
                       'bounds': [[2, 4], [1,3], [3,7], [0.5,1.5], [0.5, 1.5],
                                  [300, 500], [300, 500], [300, 500], [300, 500]],
                       # 'groups': ['G1', 'G2', 'G1'], # Sobol and Morris only. See Advanced Examples.
                       # 'dists': ['unif', 'lognorm', 'triang']
                       }            
                
        # Generate samples
        print("="*50)
        if meth=='morris':
            nlevels = 4
            X = morris.sample(problem, N, num_levels=nlevels)
            print("%i Morris trajectories, %i variables: N*(D+1) = %i simulations"%(N, problem['num_vars'], len(X)))
        elif meth=='sobol':
            second_order = False
            X = saltelli.sample(problem, N, calc_second_order=second_order)
            print("Sobol, N*(D+2) = %i simulations"%len(X))
        print("="*50)
            
        FEAT = []
        startt = timer()
        # Run GSA
        for ii, xx in enumerate(X):
            print("Simulation # %i/%i"%(ii+1, len(X)))
            pdict = dict(zip(problem['names'], xx))
            self.overrideParam(pdict)
            self.run(compute=True)
            # get outputs
            self.fetchGLSTAT(comp=['IE', 't'])
            self.fetchMATSUM(comp=['t', 'xrbv'])
            self.fetchSimulationResults()
            FEAT.append(self.getGSAoutput())
        
        endt = timer()
        print("Took %g s to run %i simulations"%(endt - startt, len(X)))
        
        # REARRANGE FEATURES STORING
        feat = {}
        for kk in FEAT[0].keys():
            feat[kk] = np.squeeze([ffeat[kk] for ffeat in FEAT])  # squeezee instead of array to reduce dimension

        # Store results
        SI = {}
        self.GSA[meth] = {'Y':feat, 'X':X, 'Si':SI}

        # Analyse outputs
        for kk in feat:
            print('='*80)
            print("**%s**"%kk)
            if meth=='morris':
                si = morrisa.analyze(problem, X, feat[kk], conf_level=0.95,
                                     print_to_console=True, num_levels=nlevels)
            elif meth=='sobol':
                si = sobol.analyze(problem, feat[kk], calc_second_order=second_order, 
                                   num_resamples=100, conf_level=0.95, print_to_console=True)
            SI[kk] = si
        
        
    def plotGSAmorris(self, star=True):
        """2D plot of the sensitivity indices of Morris analysis (mu vs sigma)
        
        """
        morris = self.GSA['morris']
        # XXX add star option
        for kk in morris['Si']:
            plt.figure('morris_%s'%kk)
            plt.title(kk)
            # fetch variables to plot
            X = morris['Si'][kk]['mu_star']
            Xrr = morris['Si'][kk]['mu_star_conf']
            Y = morris['Si'][kk]['sigma']
            ann = morris['Si'][kk]['names']
            # plot !
            plt.errorbar(X, Y, yerr=None, xerr=Xrr, fmt='+', ms=10, lw=2, ecolor='0.7')
            plt.axline((0,0), slope=1, color='0.8')
            # plt.plot(X, Y, '+', ms=10)
            for xx, yy, nn in zip(X, Y, ann):
                plt.annotate(nn, (xx, yy))
                
            plt.xlabel('$\\mu^*$')
            plt.ylabel('$\\sigma$')
            # plt.box(False) # not such a good option here
            
        
    def plotGSA(self, meth):
        """Plot sensitivity indices of GSA. Plotting method from SALib.
        
        """
        for kk in self.GSA[meth]['Si']:
            ax = self.GSA[meth]['Si'][kk].plot()
            if meth in ('sobol', 'ff'):
                # subplot
                for ii, axx in enumerate(ax):
                    axx.set_title(kk)
                    axx.axhline(y=0, color='0.8', zorder=0)

            else:
                # no subplot
                ax.set_title(kk)
                ax.axhline(y=0, color='0.8', zorder=0)        
    

if __name__=="__main__":
    plt.close('all')
    
    #%% TEST LSDYNAmodel CLASS
    param = {'tbumper':3, 'troof':2, 'trailb':2, 'trailf':5, 'tgrill':1, 'thood':1}
    if False:
        
        Car = LSDYNAmodel('./lsopt_car6/main_v21.k', param)
        # Car = LSDYNAmodel('./lsopt_car6/car6_crash_v2.k', ['test'], ['t'], ['-'], [0])
        Car._checkMainK()
        Car.overrideParam(param)
        Car.run()
    
    #%% TEST VariableScalar CLASS
    if False:
        V1 = VariableScalar()
        V1.addVariable('Mass', 'm', 56.4, 'kg')
        print(V1)
        
        TS = VariableTimeSeries()
        TS.addVariable('signal', 'S', np.linspace(0,8), np.sin(np.linspace(0,8)), 's', 'm')
        print(TS)
        TS.plot('S')
    
    
    #%% TEST CAR6model CLASS
    if False:
        CAR = CAR6model('./lsopt_car6_v2/main_v222.k', param)
        CAR.run(compute=True)
        CAR.fetchMATSUM()
        CAR.fetchGLSTAT()
        CAR.plotMATSUM('HE')
        CAR.plotGLSTAT()
        CAR.fetchSimulationResults()
        
        # CAR.runGSA(N=5, meth='morris')
        # CAR.plotGSAmorris()
        
        # CAR.runGSA(N=128, meth='sobol')
        # CAR.plotGSA('sobol')

        # compare with Sobol
        # CAR.runGSA(N=10, meth='sobol')
        # CAR.plotGSA(meth='sobol')

    # with open('CAR_sobol_n128.pickle', 'wb') as f:
    # # Pickle the 'data' dictionary using the highest protocol available.
    #     pickle.dump(CAR, f, pickle.HIGHEST_PROTOCOL)
    
    if True:
        CAR = CAR6model('./lsopt_car6_v3/main_v223.k', param)
        CAR.runGSA(N=10, meth='morris')
        CAR.plotGSAmorris()