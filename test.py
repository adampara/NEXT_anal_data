#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 17:46:35 2017

@author: para
"""
import os
import glob

filedir = '.'
filelist = os.listdir(filedir)
print filelist
filt = glob.glob('*hst')
print filt