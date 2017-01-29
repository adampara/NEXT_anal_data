#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 12:13:42 2017
analyze SiPM response as a function of radius
@author: para
"""
import math
from histograms_from_file import ldir
from ROOT import TProfile, TFile,gDirectory, TH2F
from plot_histogram import plot_histogram

debug = True

f_SiPM = 'SiPM_map.txt'
fs = open(f_SiPM)
SiPM_map = fs.read().splitlines()
print len(SiPM_map)

radius = {}
xpos = {}
ypos = {}

for line in SiPM_map:
    print len(line), line
    if len(line)>0:
        tok = line.split('|')
        if len(tok)>0:
            print len(tok[1]), tok[1], tok[4], tok[5]
            radius[tok[1].strip()] = math.sqrt(int(tok[4])**2 + int(tok[5])**2)
            xpos[tok[1].strip()] = float(tok[4])
            ypos[tok[1].strip()] = float(tok[5])
#print radius
  
an_fil = 'run_3112'     
hist_file = an_fil+'SiPM.hst'

MyFile = TFile(hist_file)

gDirectory.pwd()
dir = gDirectory
if debug:    
    print 'file directory  ',dir

hist = ldir(dir,MyFile)

 
if debug:   
    # print hists
    print 'list of histograms'
    #print hist
    
prof = TProfile('rad','signal vs radius', 250, 0., 250. )
plan = TH2F('plane',' x vs y ',50,-250.,250.,50,-250,250)
for hh in hist:
    #print hh
    h = hh.ReadObj()
    #print h
    hname = h.GetName()
    tok = hname.split('_')
    print h.GetName()
    print hname, tok
    rad = -999.
    if tok[0] == 'sipm':
        try:
            rad = radius[tok[1]]
            plan.Fill(xpos[tok[1]],ypos[tok[1]],h.GetMean()) 
            prof.Fill(rad,h.GetMean())

        except KeyError:
            print 'missing SiPM',tok[1]
            pass
        print hname, rad
    
    if rad > 175. and rad < 200.:
        print 'radius ',rad
        #plot_histogram(h)
plot_histogram(prof)
plot_histogram(plan,option='LEGO')