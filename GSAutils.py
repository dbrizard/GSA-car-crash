#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A few utilities to handle GSA results from SAlib.


Created on Mon Dec 18 15:19:38 2023

@author: dbrizard
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import figutils as fu  # only to save figures


class SobolIndices:
    """
    A class to handle Sobol Indices
    """
    
    def __init__(self, fname, skiprows=3, max_rows=None):
        """
        
        """
        self._readSAlibOuput(fname, skiprows=skiprows, max_rows=max_rows)
        
    def _readSAlibOuput(self, fname, skiprows=3, max_rows=None):
        """
        
        :param str fname: name of file to read
        """
        param = np.loadtxt(fname, delimiter=' ', usecols=[0], dtype=str,
                           skiprows=skiprows, max_rows=max_rows)
        SI = np.loadtxt(fname, delimiter=' ', usecols=[1,2,3,4],
                        skiprows=skiprows, max_rows=max_rows)
        # S1, S1_conf, ST, ST_conf
        self.SI = SI
        self.param = param.tolist()
    
    
    def plotSTS1(self, figname=None, xmargin=0.2, xoffset=0.05):
        """
        
        :param str figname: name for the figure
        :param float xmargin: additional margin for xlim
        :param float xoffset: x offset between S_1 and S_T
        """
        x = np.arange(len(self.param))
        plt.figure()
        plt.axhline(y=0, color='0.8', zorder=0)
        plt.axhline(y=1, color='0.8', zorder=0)
        plt.errorbar(x, self.SI[:,2], yerr=self.SI[:, 3], label='$S_T$', 
                     capsize=2, marker='.', linestyle='')
        plt.errorbar(x+xoffset, self.SI[:,0], yerr=self.SI[:, 1], label='$S_1$',
                     capsize=3, marker='.', linestyle='')
        
        plt.xticks(ticks=x, labels=self.param)
        plt.legend()
        
        plt.xlim(xmin=x.min()-xmargin, xmax=x.max()+xmargin)
        
        
    def _plot(self, si='ST', xoffset=0, label=None):
        """
        
        """
        x = np.arange(len(self.param))
        if si=="ST":
            ind = 2
            if label is None:
                label = '$S_T$'
        elif si=="S1":
            ind = 0
            if label is None:
                label = '$S_1$'
        plt.errorbar(x+xoffset, self.SI[:,ind], yerr=self.SI[:,ind+1], label=label,
                     capsize=2, marker='.', linestyle='')
        

class MorrisResults:
    """
    A class to handle Morris analysis results
    """
    
    def __init__(self, fname, skiprows=3, max_rows=None, outname=''):
        """
        
        :param str fname: name of file to read
        :param int skiprows: given to :func:`np.loadtxt`
        :param int max_rows: given to :func:`np.loadtxt`
        :param str outname: name of the output
        """
        self._readSAlibOutput(fname, skiprows=skiprows, max_rows=max_rows)
        self.outname = outname
        self.sortbyMuStar()
        
        
    def _readSAlibOutput(self, fname, skiprows=3, max_rows=None):
        """
        
        :param str fname: name of file to read
        """
        param = np.loadtxt(fname, usecols=[0], dtype=str,
                           skiprows=skiprows, max_rows=max_rows)
        SI = np.loadtxt(fname, usecols=[1,2,3,4],
                        skiprows=skiprows, max_rows=max_rows)
        # mu_star, mu, m_star_conf, sigma
        self.SI = SI
        self.param = param.tolist()
        
    
    def sortbyMuStar(self,):
        """Sort input parameters with decreasing mu_star
        
        """
        ind = self.SI[:,0].argsort()
        self.muStarSort = {'param':[self.param[ii] for ii in ind], 'ind':ind}
        


    def plotMorris(self, figname=None, conf=True):
        """Plot Morris graph (*sigma* vs *mu_star*)
        
        :param str figname: name for the figure
        :param bool conf: plot confidence interval *mu_star_conf*
        """       
        plt.figure(figname)
        plt.title(self.outname)
        
        X = self.SI[:,0]  # mu_star
        Y = self.SI[:,3]  # sigma
        if conf:
            Xrr = self.SI[:,2]  # mu_star_conf
        else:
            Xrr = 0

        plt.errorbar(X, Y, yerr=None, xerr=Xrr, fmt='+', ms=10, lw=2, ecolor='0.7')
        plt.axline((0,0), slope=1, color='0.8')
        # plt.plot(X, Y, '+', ms=10)
        for xx, yy, nn in zip(X, Y, self.param):
            plt.annotate(nn, (xx, yy), xytext=(2,2), textcoords='offset pixels')
            
        plt.xlabel('$\\mu^*$')
        plt.ylabel('$\\sigma$')
        
        
    def plotMorris_color(self, figname=None, conf=False, nofig=False,
                         xlim=None, ylim=None):
        """Plot Morris graph (*sigma* vs *mu_star*) with colored markers
        
        :param str figname: name for the figure
        :param bool conf: plot confidence interval *mu_star_conf*
        :param bool nofig: do not create a new figure (for subfig purpose, see :meth:`GatherMorris.subplot2D`)
        """
        if not nofig:
            plt.figure(figname)
        plt.title(self.outname)
        
        X = self.SI[:,0]
        Y = self.SI[:,3]
        if conf:
            Xrr = self.SI[:,2]
        else:
            Xrr = [0]*len(X)

        for ii, (xx, yy, rr, nn) in enumerate(zip(X, Y, Xrr, self.param)):
            plt.errorbar(xx, yy, yerr=None, xerr=rr, fmt='+', ms=10, lw=2, ecolor='0.7', mew=2)
            plt.annotate(nn, (xx, yy), xytext=(2,2), textcoords='offset pixels', color='C%i'%ii)
        plt.axline((0,0), slope=1, color='0.8')
            
        plt.xlabel('$\\mu^*$')
        plt.ylabel('$\\sigma$')     
        
        plt.xlim(xlim)
        plt.ylim(ylim)



class GatherSobol:
    """
    A class to handle a list of :class:`SobolIndices`
    """
    
    def __init__(self, SIlist, pname, values):
        """
        
        :param list SIlist: list of :class:`SobolIndices` objects
        :param str pname: name of the parameter
        :param list values: values of the parameter
        """
        self.SIlist = SIlist
        self.param = {'name':pname, 'values':values}
        
        self._aggregate()
    
    
    def _aggregate(self, ):
        """Fetch data and store in a convenient way
        
        """
        self.S1 = np.array([temp.SI[:,0] for temp in self.SIlist]).T
        self.S1_conf = np.array([temp.SI[:,1] for temp in self.SIlist]).T
        self.ST = np.array([temp.SI[:,2] for temp in self.SIlist]).T
        self.ST_conf = np.array([temp.SI[:,3] for temp in self.SIlist]).T
        
    
    def plot(self, figname=None, xmargin=0.2, conf=True):
        """Plot S1 and ST wrt parameter
        
        :param str figname: name for the figure
        :param float xmargin: additional margin for xlim
        :param bool conf: plot confidence intervals
        """
        X = self.param['values']
        plt.figure(figname)

        ZZ = zip([211, 212],
                 ['$S_T$', '$S_1$'],
                 [self.ST, self.S1],
                 [self.ST_conf, self.S1_conf],
                 [self.SIlist[0].param, self.SIlist[0].param])
        
        for sbplt, ylab, SS, SSCC, pp in ZZ:
            plt.subplot(sbplt)
            for ss, sscc, ppp in zip(SS, SSCC, pp):
                if conf:
                    plt.errorbar(X, ss, yerr=sscc, label=ppp,
                                 capsize=3, marker='.')
                else:
                    plt.plot(X, ss, marker='.', label=ppp)

            plt.xlabel(self.param['name'])
            plt.ylabel(ylab)
            plt.xlim(xmin=np.min(X)-xmargin, xmax=np.max(X)+xmargin)
            
            plt.axhline(y=0, color='0.8')
            plt.axhline(y=1, color='0.8')
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        

        # for S1, S1_conf, pp in zip(self.S1, self.S1_conf, self.SIlist[0].param):
        #     if conf:
        #         plt.errorbar(X, S1, yerr=S1_conf, label=pp,
        #                      capsize=3, marker='.')
        #     else:
        #         plt.plot(X, S1, marker='.', label=pp)

    
    def plotSTS1(self, si='ST', figname=None, conf=True):
        """
        
        """
        n = len(self.SIlist)
        xoffset = 1/n/2
        
        plt.figure()
        plt.axhline(y=0, color='0.8', zorder=0)
        plt.axhline(y=1, color='0.8', zorder=0)
        for ii, (SI, vv) in enumerate(zip(self.SIlist, self.param['values'])):
            SI._plot(si=si, xoffset=ii*xoffset, label=vv)
        
        x = np.arange(len(self.SIlist[0].param))
        plt.xticks(ticks=x, labels=self.SIlist[0].param)
        plt.legend(title=self.param['name'])
        
        plt.xlim(xmin=x.min()-0.2)
        plt.ylabel(si)


class GatherMorris:
    """
    A class to handle a list of :meth:`MorrisResults`
    """
    
    def __init__(self, MOlist, pname, values, outname):
        """
        
        :param list SIlist: list of :class:`MorrisResults` objects
        :param str pname: name of the parameter
        :param list values: values of the parameter  
        :param str outname: name of the output
        """
        self.MOlist = MOlist
        self.param = {'name':pname, 'values':values}
        self.outname = outname
        
        self._aggregate()
    
    
    def _aggregate(self, ):
        """Fetch data and store in a convenient way
        
        """
        self.mu_star = np.array([temp.SI[:,0] for temp in self.MOlist]).T
        self.mu = np.array([temp.SI[:,1] for temp in self.MOlist]).T
        self.mu_star_conf = np.array([temp.SI[:,2] for temp in self.MOlist]).T
        self.sigma = np.array([temp.SI[:,3] for temp in self.MOlist]).T  # XXX .T seems useless...
        
        self.paramSort = np.array([temp.muStarSort['param'] for temp in self.MOlist]).T
        self.indSort = np.array([temp.muStarSort['ind'] for temp in self.MOlist]).T
        
        
    def plot(self, figname=None, xmargin=0.2, conf=True):
        """Plot mu_star and sigma wrt param value (subplot)
        
        :param str figname: name for the figure
        :param float xmargin: additional margin for xlim
        :param bool conf: plot confidence intervals
        """
        X = self.param['values']
        plt.figure(figname)
        
        plt.subplot(211)
        plt.title(self.outname)
        for mu_star, mu_star_conf, pp in zip(self.mu_star, self.mu_star_conf, self.MOlist[0].param):
            if conf:
                plt.errorbar(X, mu_star, yerr=mu_star_conf, label=pp,
                             capsize=3, marker='.')
            else:
                plt.plot(X, mu_star, '.-', label=pp)
            
            plt.legend()
            plt.xlabel(self.param['name'])
            plt.ylabel('$\\mu^*$')
            plt.xlim(xmin=np.min(X)-xmargin, xmax=np.max(X)+xmargin)

        plt.subplot(212)
        plt.plot(X, self.sigma.T, '.-')
        plt.legend(self.MOlist[0].param)
        plt.xlim(xmin=np.min(X)-xmargin, xmax=np.max(X)+xmargin)
        plt.xlabel(self.param['name'])
        plt.ylabel('$\\sigma$')
        # for sigma, pp in zip(self.sigma, self.MOlist[0].param):
        #     plt.plot(X, sigma)


    def plot2D(self, figname=None, conf=True, pointID=True, xlim0=True, xlim=None, ylim=None):
        """Gather sigma vs mu_star plots for all the values of param
        
        :param str figname: name for the figure
        :param bool conf: plot confidence intervals
        :param bool pointID: 
        :param bool xlim0: 
        :param float xlim:
        :param float ylim: 
        """
        plt.figure(figname)
        plt.title(self.outname)
        plt.axline((0,0), slope=1, color='0.8')
        
        ZZ = zip(self.mu_star, self.mu_star_conf, self.sigma, self.MOlist[0].param)
        
        for jj, (mu_star, mu_star_conf, sigma, pp) in enumerate(ZZ):
            if conf:
                plt.errorbar(mu_star, sigma, yerr=None, xerr=mu_star_conf, label=pp,
                             fmt='.:', lw=1, ecolor='0.8')
            else:
                plt.plot(mu_star, sigma, '.:', label=pp)
            if pointID:
                for ii, (xx, yy) in enumerate(zip(mu_star, sigma)):
                    plt.text(xx, yy, ii, color='C%i'%jj, )

        plt.legend()
        plt.xlabel('$\\mu^*$')
        plt.ylabel('$\\sigma$')
        if xlim0:
            plt.xlim(xmin=0)
        plt.xlim(xmax=xlim)
        plt.ylim(ymax=ylim)
    
    def subplot2D(self, figname=None, figsize=(19.2, 4.8), conf=False, margin=1.05,
                  xlim=None, ylim=None):
        """Plot all the sigma vs mu_star plots with subplot
        
        :param str figname: name for the figure
        :param tupple figsize: (width, height) of the firure
        :param bool conf: plot confidence intervals
        :param float margin: margin for x and y lim
        """
        plt.figure(figname, figsize=figsize)
        for ii, MO in enumerate(self.MOlist):
            # determine xlim and ylim
            if xlim is None:
                xlim = margin*self.mu_star.max()
            if ylim is None:
                ylim = margin*self.sigma.max()
            # plot Morris graphs
            plt.subplot(1,len(self.MOlist),ii+1)
            MO.plotMorris_color(conf=conf, nofig=True, xlim=(0, xlim), ylim=(0, ylim))
            if ii==0:
                plt.title(self.outname)
        
    
    def plotRanking(self, figname=None, title=None, ticks=False):
        """Plot input ranking for all the Morris repetitions
        
        https://stackoverflow.com/questions/40887753/display-matrix-values-and-colormap
        
        :param str figname: name for the figure
        :param str title: title for the plot
        :param bool ticks: 
        """
        matrix = np.flipud(self.indSort) + 1

        nparam = matrix.shape[0]
        nrepet = matrix.shape[1]
        
        fig, ax = plt.subplots(num=figname, figsize=(nrepet/3, nparam/3))
        ax.matshow(matrix, cmap='viridis_r')
        
        plt.title(title)
        plt.xlabel('repetition')
        plt.ylabel('ranking')
        if ticks:
            plt.xticks(range(nrepet), np.arange(nrepet)+1)
            plt.yticks(range(nparam), np.arange(nparam)+1)
        else:
            plt.xticks([])
            plt.yticks([])
        plt.box(False)
        
        for ii in range(nrepet):
            for jj in range(nparam):
                ax.text(ii, jj, str(matrix[jj,ii]), va='center', ha='center')
        


if __name__=="__main__":
    plt.close('all')
    
    if False:
        # %% TEST SOBOLINDICES
        SI = SobolIndices('gsautils_sobol.csv')
        SI.plotSTS1()
        
        # %% TEST MORRISRESULTS
        MO1 = MorrisResults('gsautils_morris.csv')
        MO1.plotMorris(conf=False)
        
        MO2 = MorrisResults('gsautils_morris2.csv')
        MO2.plotMorris()
        
        # %% TEST GATHERSOBOL
        GS = GatherSobol([SI, SI], 'dummy param', [4, 5])
        GS.plot()
        GS.plot(conf=False)
        
        # %% TEST GATHERMORRIS
        GM = GatherMorris([MO1, MO2], 'dum. p.', [4, 5])
        GM.plot()
        GM.plot2D()
        GM.plot2D(conf=False)
    
    # %% TRY IN REAL SITUATION!
    if False:
        file = '/home/dbrizard/Calcul/24_car/GSA/car_v222_right-impact/sobol-output.log.md'
        
        out = ['fmax', 'dmax', 'vfin']
        offset = [0, 9, 18]
        for outt, of in zip(out, offset):
            # **fmax**
            so1 = SobolIndices(file, skiprows=4+of, max_rows=6)
            so2 = SobolIndices(file, skiprows=34+of, max_rows=6)
            so3 = SobolIndices(file, skiprows=64+of, max_rows=6)
            so4 = SobolIndices(file, skiprows=94+of, max_rows=6)
            so5 = SobolIndices(file, skiprows=124+of, max_rows=6)
            
            SO = GatherSobol([so1, so2, so3, so4, so5], 'N', 
                             [64, 128, 130, 256, 512])
            SO.plot(figname=outt, xmargin=10)
            SO.plotSTS1(si='ST')
            SO.plotSTS1(si='S1')

    
    # %% Morris n5, n10 and n20 repetitions
    if True:
        ntraj = 20
        model = 'v223'
        if model=='v222':
            # OLD MODEL
            offset = [0, 9, 18, 27]  # offset to fecht results for each output variable
            nlignes = 39  # number of lignes for 1 set of morris indices for all outputs
            nparam = 6
            if ntraj==10:
                file = '/home/dbrizard/Calcul/25_car/GSA_old/car_v222_right-impact_v30/morris-n10-output.md'
                nrep = 6
            elif ntraj==5:
                file = '/home/dbrizard/Calcul/25_car/GSA_old/car_v222_right-impact_v30/morris-n5-output.md'
                nrep = 6
        elif model=='v223':
            # CURRENT MODEL
            offset = [0, 12, 24, 36]  # offset to fecht results for each output variable
            nlignes = 51  # number of lignes for 1 set of morris indices for all outputs
            nparam = 9  # number of uncertain parameters
            folder = '/home/dbrizard/Calcul/25_car/GSA/car_v223_right-impact_v30/'
            file = os.path.join(folder, 'morris_n%i_output.md'%ntraj)
            NREP = {5:6, 10:6, 20:6}  # number of repetitions of the morris analysis wrt number of trajectories
            nrep = NREP[ntraj]
            ffolder = 'GSA/car_v223_right-impact_v30/gathermorris_n%i/'%ntraj  # folder for saving fig output
            xlim = {'dmax':None, 'fmax':700000, 'IE':9e6, 'vfin':1700}
            ylim = {'dmax':None, 'fmax':370000, 'IE':7.5e6, 'vfin':900}
        elif model=='v223_4':
            # CURRENT MODEL but 4 param instead of 9
            offset = [0, 7, 14, 21]  # offset to fecht results for each output variable
            nlignes = 31  # number of lignes for 1 set of morris indices for all outputs
            nparam = 4  # number of uncertain parameters
            folder = '/home/dbrizard/Calcul/25_car/GSA/car_v223_right-impact_v30_prob4/'
            file = os.path.join(folder, 'morris_n%i_output.md'%ntraj)
            NREP = {5:6, 10:6, 20:6}  # number of repetitions of the morris analysis wrt number of trajectories
            nrep = NREP[ntraj]
            ffolder = 'GSA/car_v223_right-impact_v30_prob4/gathermorris_n%i/'%ntraj  # folder for saving fig output            
            xlim = {'dmax':None, 'fmax':700000, 'IE':8.5e6, 'vfin':1800}
            ylim = {'dmax':None, 'fmax':510000, 'IE':7e6, 'vfin':1100}
                

        out = ['fmax', 'dmax', 'vfin', 'IE']
        GM = []
        for outt, of in zip(out, offset):
            MOlist = []
            for ii in range(nrep):
                mo1 = MorrisResults(file, skiprows=ii*nlignes+4+of, max_rows=nparam)
                MOlist.append(mo1)
            
            MO = GatherMorris(MOlist, 'rep', [ii for ii in range(len(MOlist))], outt)
            MO.plot2D(figname='morris_%i_%s'%(ntraj, outt), conf=False, xlim=xlim[outt], ylim=ylim[outt])
            MO.subplot2D(figname='smorris_%i_%s'%(ntraj, outt), xlim=xlim[outt], ylim=ylim[outt])
            MO.plotRanking(figname='rank_%i_%s'%(ntraj, outt), title='r=%i'%ntraj)
            GM.append(MO)
        
        if True:
            fu.savefigs(ffolder, overw=True)
