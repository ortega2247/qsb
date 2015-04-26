# %load qsb_model.py
"""
Created on Wed Apr 22 01:56:25 2015

@author: lolab
"""


import csv
import matplotlib.pyplot as plt 
import numpy as np
from scipy.optimize import curve_fit
from pysb import *
from pysb.macros import *
import pylab as pl
from pysb.integrate import odesolve

# some functions to make life easy
site_name = 'b'
def catalyze_b(enz, sub, product, klist):
    """Alias for pysb.macros.catalyze with default binding site 'b'.
    """
    return catalyze(enz, site_name, sub, site_name, product, klist)

def bind_table_b(table):
    """Alias for pysb.macros.bind_table with default binding sites 'bf'.
    """
    return bind_table(table, site_name, site_name)


KF = 1e-6
KR = 1e-3
KC = 1    

# Default forward, reverse, and catalytic rates
KFc = 8.6897973656754578e-07 #1e-6
KRc = 0.00026884384823502833 #1e-3
KCc = 7.9623024965305493 #1

KFcc = 8.6897973656754578e-07 #1e-6
KRcc = 0.00026884384823502833 #1e-3
KCcc = 7.9623024965305493 #1

KFb = 4.149733352245889e-09 #1e-7
KRb = 0.00050110975371602537 #1e-3
KCb = 0.39197357243506348 #1

KFba = 2.2299766742556005e-08 #1e-6
KRba = 0.0018670344174003893 #1e-3
KCba = 0.41085725049645588 #1

KFm = 1.2128334849034545e-07 #2e-6
KRm = 0.011242081848190973 #10e-3 
KCm = 17.003767066674314 #10

#C3 claeaves PARP
KFs = 4.0379463916345119e-05 #1e-6
KRs = 0.0024327371579203338 #10e-3 
KCs = 0.012083736473176411 #1





# Bcl2 Inhibition Rates
# OLD bcl2_bid_rates = [1e-6, 1e-3] # 1.0e-6/v_mito
bcl2_bid_rates = [2.3147601898026888e-06, 0.0017312607994352624] # 1.0e-6/v_mito
bcl2_bax_rates = [9.9990412921086819e-05, 2.0740711292172149e-05]
XIAP_C3_rates =  [9.9990412921086819e-05, 2.0740711292172149e-05]
XIAP_Smac_rates = [2.3147601898026888e-06, 0.0017312607994352624]

# instantiate a model
Model()

# declare monomers
Monomer('C8', ['b'])
Monomer('Bid',    ['b', 'S'], {'S':['n', 't']})
#Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':[ 'd', 'a']})
Monomer('Bax',   ['b', 'S'], {'S':['d', 'a']})
Monomer('M',   ['b', 'S'], {'S':['c', 'p']})
Monomer('Smac', ['b', 'S'], {'S':['m', 'r']})
Monomer('PARP', ['b', 'S'],  {'S':['n', 'c']})
Monomer('Bcl2', ['b'])
Monomer('C3', ['b', 'S'], {'S':['d','a']})
Monomer('XIAP', ['b'])

synthesize(C8(b=None), 0.032866241994187327)
catalyze_b(C8, C3(S='d'), C3(S='a'), [KFcc, KRcc, KCcc])
catalyze_b(C3(S='a'), PARP(S='n'), PARP(S='c'), [KFs, KRs, KCs])
catalyze_b(C8, Bid(S='n'), Bid(S='t'), [KFc, KRc, KCc])
# Activate Bax
catalyze_b(Bid(S='t'), Bax(S='d'), Bax(S='a'), [KFb, KRb, KCb])


# Activate Bax/Bak
catalyze_b(Bax(S='a'), M(S='c'), M(S='p'), [KFba, KRba, KCba])
catalyze_b(M(S='p'), Smac(S='m'), Smac(S='r'), [KFm, KRm, KCm])
#catalyze_b(Smac(S='r'), PARP(S='n'), PARP(S='c'), [KFs, KRs, KCs])

bind_table_b([[                  Bcl2,           XIAP],
              [Bid(S='t'),   bcl2_bid_rates,     None        ],
              [Bax(S='a'),   bcl2_bax_rates,     None        ],
              [C3,           None,          XIAP_C3_rates    ],
              [Smac(S='r'), None ,         XIAP_Smac_rates  ]])

# Bid, Bax, BclxL "transport" to the membrane
#equilibrate(Bid(b=None, S='t'),   Bid(b=None, S='m'), [1e-1, 1e-3])
#equilibrate(Bax(b=None, S='m'),   Bax(b=None, S='a'), [1e-1, 1e-3])
#equilibrate(BclxL(b=None, S='c'), BclxL(b=None, S='m'), [1e-1, 1e-3])


# initial conditions
Parameter('C8_0',    0)
Parameter('Bid_0',   4e4)
Parameter('Bax_0',  1e5)
Parameter('M_0',  5e5)
Parameter('Smac_0', 1e5)
Parameter('PARP_0',  1e6)
Parameter('Bcl2_0',  2e4)
Parameter('C3_0',  2e4)
Parameter('XIAP_0',  2e4)

Initial(C8(b=None), C8_0)
Initial(Bid(b=None, S='n'), Bid_0)
Initial(Bax(b=None, S='d'), Bax_0)
Initial(M(b=None, S='c'), M_0)
Initial(Smac(b=None, S='m'), Smac_0)
Initial(PARP(b=None, S='n'), PARP_0)
Initial(Bcl2(b=None), Bcl2_0)
Initial(C3(b=None, S='d'), C3_0)
Initial(XIAP(b=None), XIAP_0)

# Observables
Observable('mBid', Bid(b=None, S='t'))
Observable('aSmac', Smac(b=None, S='r'))
Observable('cPARP', PARP(b=None, S='c'))


