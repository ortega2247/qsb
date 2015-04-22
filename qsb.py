# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/lolab/.spyder2/.temp.py
"""

import csv
import matplotlib.pyplot as plt 
import numpy as np
from scipy.optimize import curve_fit
from pysb import *
from pysb.macros import *
import pylab as pl
from pysb.integrate import odesolve

exp_data = open('/home/lolab/Downloads/EC-RP_IMS-RP_IC-RP_data_for_models.csv')
csv_exp = csv.reader(exp_data)

# IC-RP: initiator caspase reporter protein 
# IMS-RP: Mitochondrial inter-membrane space reporter protein
# EC-RP: Effector caspase reporter protein

time = []
IC_RP = []
IMS_RP = []
EC_RP = []

firstline = True
for row in csv_exp:
    if firstline:
        firstline = False
        continue
    time.append(float(row[0]))
    IC_RP.append(float(row[2]))
    IMS_RP.append(float(row[5]))
    EC_RP.append(float(row[8]))

time = np.array(time)
IC_RP = np.array(IC_RP)
IMS_RP = np.array(IMS_RP)
EC_RP = np.array(EC_RP)

plt.plot(time,IC_RP)
#plt.plot(time,IMS_RP)
#plt.plot(time,EC_RP)



#We have to define the functions that we want to fit

def linea(x,a,b):
    return a*x+b


def expon(x,a,b,c):
    return a*np.exp(b*x)+c
    
    
def poli(x, a, b, c, d):
    return a*x**3 + b*x**2 +c*x + d
    
poptl, pcovl = curve_fit(linea,time,IC_RP)
poptp, pcovp = curve_fit(poli,time,IC_RP)
poptex, pcovex = curve_fit(expon,time,IC_RP, p0=(0.1, 1e-4,10))
 
linea_fit = linea(time,*poptl)
poli_fit = poli(time,*poptp)
expo_fit = expon(time,*poptex)

plt.plot(time,linea_fit)
plt.plot(time,poli_fit)
plt.plot(time,expo_fit)
plt.ylim(0,1)

plt.show()

############Modeling enzymatic reaction with PYSB



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

# Default forward, reverse, and catalytic rates
KF = 1e-6
KR = 1e-3
KC = 1

# Bid activation rates
bid_rates = [        1e-7, 1e-3, 1] #

# Bcl2 Inhibition Rates
bcl2_rates = [1.428571e-05, 1e-3] # 1.0e-6/v_mito

# instantiate a model
Model()

# declare monomers
Monomer('C8', ['b'])
Monomer('Bid',    ['b', 'S'], {'S':['n', 't']})
Monomer('Bax',   ['b', 'S'], {'S':['d', 'a']})
Monomer('M',   ['b', 'S'], {'S':['c', 'p']})
Monomer('Smac', ['b', 'S'], {'S':['m', 'r']})
Monomer('Parp', ['b', 'S'],  {'S':['n', 'c']})


catalyze_b(C8, Bid(S='n'), Bid(S='t'), [KF, KR, KC])
# Activate Bax
catalyze_b(Bid(S='t'), Bax(S='d'), Bax(S='a'), [KF, KR, KC])

# Activate Bax/Bak
catalyze_b(Bax(S='a'), M(S='c'), M(S='p'), bid_rates)
catalyze_b(M(S='p'), Smac(S='m'), Smac(S='r'), bid_rates)
catalyze_b(Smac(S='r'), Parp(S='n'), Parp(S='c'), bid_rates)


# Bid, Bax, BclxL "transport" to the membrane
#equilibrate(Bid(b=None, S='t'),   Bid(b=None, S='m'), [1e-1, 1e-3])
#equilibrate(Bax(b=None, S='m'),   Bax(b=None, S='a'), [1e-1, 1e-3])
#equilibrate(BclxL(b=None, S='c'), BclxL(b=None, S='m'), [1e-1, 1e-3])


#bind_table_b([[                  Bcl2,  BclxL(S='m'),       Mcl1],
#              [Bid(S='m'), bcl2_rates,  bcl2_rates,   bcl2_rates],
#              [Bax(S='a'), bcl2_rates,  bcl2_rates,         None],
#              [Bak(S='a'),       None,  bcl2_rates,   bcl2_rates]])

# initial conditions
Parameter('C8_0',    1e4)
Parameter('Bid_0',   1e4)
Parameter('Bax_0',  .8e5)
Parameter('M_0',  .2e5)
Parameter('Smac_0', 1e3)
Parameter('Parp_0',  1e3)
==
Initial(C8(b=None), C8_0)
Initial(Bid(b=None, S='n'), Bid_0)
Initial(Bax(b=None, S='d'), Bax_0)
Initial(M(b=None, S='c'), M_0)
Initial(Smac(b=None, S='m'), Smac_0)
Initial(Parp(b=None, S='n'), Parp_0)

# Observables
Observable('obstBid', Bid(b=None, S='t'))
Observable('obsBax', Bax(b=None, S='a'))
Observable('obsSmac', Smac(b=None, S='r'))
Observable('obsParp', Parp(b=None, S='c'))


t = pl.linspace(0, 20000)
yout = odesolve(model, t)

pl.ion()
pl.figure()
pl.plot(t, yout['obstBid'], label="tBid")
pl.plot(t, yout['obsC8'], label="C8")
pl.plot(t, yout['obsBax'], label='activated bax')
pl.plot(t, yout['obsSmac'], label='Smac')
pl.plot(t, yout['obsParp'], label='cParp')
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules/cell")
pl.show()