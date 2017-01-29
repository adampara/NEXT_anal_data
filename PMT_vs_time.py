#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 20:12:25 2017
Analyze PMT information as a function of time

@author: para
"""
import argparse
from ROOT import TFile, gDirectory
from plot_histogram import plot_histogram
from histograms_from_file import ldir
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
    
import sys
import glob
from plot_histogram import plot_histogram

err
an_fil = PMT_vs_time_inp()
print an_fil
filelist = glob.glob(an_fil+'_*hst')
for ff in filelist:
    print ff
    fnum = ff.split('_')[2].split('.')[0]
    print fnum
    MyFile = TFile(ff)  
    gDirectory.pwd()
    dir = gDirectory
    if debug:    

    lhist = ldir(dir,MyFile)
    for hh in lhist:
         h = hh.Readobj()
         plot_histogram(h)
    
    MyFile.Close()
    