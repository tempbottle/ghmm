#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
'''
file: testghmm.py
author: xjump.me#at#gmail#dot#com
'''

from ghmm import *
sigma = IntegerRange(1,7)
A = [[0.9,0.1],[0.3,0.7]]
efair = [3.0/13, 3.0/13, 2.0/13, 2.0/13, 2.0/13, 1.0/13]
eloaded = efair
efair = [1.0/6]*6
B = [efair, eloaded]
pi = [0.5]*2
m = HMMFromMatrices(sigma, DiscreteDistribution(sigma), A, B, pi)
print m
