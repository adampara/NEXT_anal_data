#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 19:31:51 2017

@author: para
"""
import sys
import os


def Spyder():

    return any('SPYDER' in name for name in os.environ)
    
if sys.flags.interactive:
    print 'interactive'
else:
    print ' not interactive '
    

if not 'PYTHONSTARTUP' in os.environ:
    print 'not'
else:
    print 'yes'
    
print os.environ

for name in os.environ:
    print name

if any('SPYDER' in name for name in os.environ):
    print 'aa'
else:
    print 'bb'
 #   print 'any'
#else:        
 #   print 'else'

print 'interactive? ', Spyder()