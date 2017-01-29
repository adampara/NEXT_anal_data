

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 20:12:25 2017
Analyze PMT information as a function of time

@author: para
"""
import argparse
from ROOT import TFile, gDirectory, TH1F
from plot_histogram import plot_histogram
from histograms_from_file import ldir
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np

def PMT_vs_time_inp():
    """
    parse input pvalues
    """
    
#    parse input arguments for analysis 

    parser = argparse.ArgumentParser()
    parser.add_argument("-an_fil", type=str, help="files to analyze",
                        default='run_3112',dest='an_fil')

    
    ar = parser.parse_args()
    
    
    print '---------------------------------------------------------------------'
    print '  files to analyze               ',ar.an_fil
    print '---------------------------------------------------------------------'
     
    return ar.an_fil
    
    
def h_stat(h):
    """
    return the mean value and its error for a given histogram
    """
    n_PMT = 12
    mean = h.GetMean()
    RMS = h.GetRMS()
    ent = h.GetEntries()/n_PMT
    err = RMS/sqrt(ent)
    
    return mean, err
    
import sys
import glob

an_fil = PMT_vs_time_inp()
print an_fil
filelist = glob.glob(an_fil+'_*hst')

cont = ['t_pe_post', 'q_pe_post', 'n_pe_post',  'ampl_post',
        't_pe_front', 'q_pe_front', 'n_pe_front',  'ampl_front',
        'q_sig' ]

title = ['mean arrival time of pe, late', 'charge of pe, late',
         'number of pe, late', 'amplitude of pe, late',
         'mean arrival time of pe, early', 'charge of pe, early',
         'number of pe, early', 'amplitude of pe, early',
         'charge, S2 signal']
title_dict = dict(zip(cont,title))
# y_lower = [ 600, 0, 0, 0, 0, 0, 0, 0, 0, ]
set_mean = {}
set_err = {}
time = []
for ff in filelist:
    print ff
    fnum = float(ff.split('_')[2].split('.')[0])
    print fnum
    time.append(fnum)
    
    MyFile = TFile(ff)  
    gDirectory.pwd()
    dir = gDirectory

    lhist = ldir(dir,MyFile)
    for hh in lhist:
        name = hh.GetName()
        if name in set_mean:
            pass
        else:
            set_mean[name] = []
        if name in set_err:
            pass
        else:
            set_err[name] = []
                
        h = hh.ReadObj()
        (mean,err) = h_stat(h)
        print mean,err
        set_mean[name].append(mean)
        set_err[name].append(err)
        #plot_histogram(h)
    
    MyFile.Close()

print time
ifig = 0
for name in set_mean:
    ifig += 1
    print name
    print set_mean[name]
    plt.figure(ifig)
    # plt.plot(time, np.array(set_mean[name])/set_mean[name][-1], marker='*')
    #plt.errorbar(time, np.array(set_mean[name])/set_mean[name][-1],
    #             yerr=np.array(set_err[name])/set_mean[name][-1], marker='*')
    ymean = np.mean(np.array(set_mean[name][-5:]))
    for t,y,dy in zip(time,np.array(set_mean[name]),np.array(set_err[name])):
        plt.errorbar(t,y/ymean,yerr=dy/ymean,marker='*', color='red')
    plt.title(title_dict[name])
    plt.ylim([0.7,1.1])
    plt.xlim([time[0]-5,time[-1]+5])
    plt.grid(True)
    plt.xlabel('Time in the run (0 - 100)')
    plt.ylabel('arb. units, rnd of run = 1.')
    
ifig += 1
plt.figure(ifig)
name = 'n_pe_post'
plt.errorbar(time, np.array(set_mean[name]),
               yerr=np.array(set_err[name]), marker='*', linestyle='None')
name = 'n_pe_front'
plt.errorbar(time, np.array(set_mean[name]),
               yerr=np.array(set_err[name]), marker='*', linestyle='None')
plt.xlim([time[0]-5,time[-1]+5])
plt.grid(True)
plt.xlabel('Time in the run (0 - 100)')
plt.ylabel('average number of photoelectrons')

ifig += 1
plt.figure(ifig)
name = 'q_pe_post'
plt.errorbar(time, np.array(set_mean[name]),
               yerr=np.array(set_err[name]), marker='*', linestyle='None')
name = 'q_pe_front'
plt.errorbar(time, np.array(set_mean[name]),
               yerr=np.array(set_err[name]), marker='*', linestyle='None')
plt.xlim([time[0]-5,time[-1]+5])
plt.ylim([0,25])
plt.grid(True)
plt.xlabel('Time in the run (0 - 100)')
plt.ylabel('average charge of photoelectrons, arb. units')

ifig += 1
plt.figure(ifig)
name = 'ampl_post'
plt.errorbar(time, np.array(set_mean[name]),
               yerr=np.array(set_err[name]), marker='*', linestyle='None')
name = 'ampl_front'
plt.errorbar(time, np.array(set_mean[name]),
               yerr=np.array(set_err[name]), marker='*', linestyle='None')
plt.xlim([time[0]-5,time[-1]+5])
plt.ylim([0,5])
plt.grid(True)
plt.xlabel('Time in the run (0 - 100)')
plt.ylabel('average amplitude of photoelectrons, arb. units')

plt.show()