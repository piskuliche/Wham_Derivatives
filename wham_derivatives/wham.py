#!/usr/bin/env python
import numpy as np
from numpy import inf,nan
import argparse,pickle,os
import matplotlib.pyplot as plt

"""
This is a simple Weighted Histogram Analysis Code that can be used to calculate the Potential of Mean Force.
Please read the readme file for more specific notes about its operation.

Copyright 2021, Zeke Piskulich, University of Kansas
For questions, concerns, email piskulichz@gmail.com
"""

# Constants
kb=0.0019872041
T=298.15
kbT=kb*T

def choose_bias(opt):
    if opt == 1:
        def calc_bias(x,k):
            # This function calculates the bias potential value 
            return k*x**2.
        return calc_bias
    elif opt == 2:
        def calc_bias(x,k):
            # This function calculates the bias potential value
            return 0.5*k*x**2.
        return calc_bias
    else:
        exit("improper bias")

def calc_p(F):
    # This function calculates the unweighted probability distribution
    numerator=np.sum(ni,axis=0)
    denominator=np.sum(np.multiply(cnt,np.exp(np.divide(-np.subtract(U.T,F-np.min(F)),kbT))),axis=1)
    return numerator/denominator

def calc_f(P):
    # This calculates the weight of each window (and then takes the log to get the free energy)
    F=-kbT*np.log(np.sum(np.multiply(P,np.exp(-U/kbT)),axis=1))
    F[F==inf] = 0.0
    return F

def read_ener(etypes,window):
    energy={}
    if shouldrerun == 0:
        for ener in etypes:
            energy[ener]=np.genfromtxt(str(window)+"/"+ener+"_init.out",usecols=0,unpack=True)
        pickle.dump(open("wham_pckl/energy_+%d.pckl"%window,'wb'))
    else:
        energy = pickle.load(open("wham_pckl/energy_%d.pckl"%window,'rb'))
    return energy

def calc_dp(P,F,dpbias,dF):
    denominator=np.sum(np.multiply(cnt,np.exp(np.divide(-np.subtract(U.T,F-np.min(F)),kbT))),axis=1)
    term1=np.sum(np.multiply(cnt,dpbias.T),axis=1)
    # Somewhere in term2 is the thing that is too big.
    term2=P*np.sum(np.multiply(np.multiply(cnt,np.exp(np.divide(-np.subtract(U.T,F-np.min(F)),kbT))),(dF/kbT+F-np.min(F)-U.T)),axis=1)
    #print(term1,term2)
    return (term1 - term2)/denominator

def calc_dF(P,F,dP):
    term1=-kbT*F
    term2num=kbT*np.sum(np.multiply(np.subtract(dP,np.multiply(U,P)),np.exp(-U/kbT)),axis=1)
    term2den=np.sum(np.exp(-U/kbT)*P,axis=1)
    dF = term1+term2num/term2den
    print(term1,term2num/term2den)
    dF[dF==inf] == 0.0
    return dF


def wham_iteration(F):
    P=calc_p(F)
    F=calc_f(P)
    # Calculates and writes out the error
    Ferr = np.sum(np.abs(np.subtract(F,F_old)))
    print("Error: %s " % Ferr)
    # Checks if convered, if it is, writes files and exits calculation.
    isconverged=False
    if (Ferr < tolerance):
        isconverged=True
        np.savetxt("wham_pmf.out", np.c_[xvals,-kbT*np.log(P)-np.min(-kbT*np.log(P)),P])
        np.savetxt("wham_F.out", np.c_[F-np.min(F)])
    return F, isconverged

def dwham_iteration(F,dF,key):
    P=calc_p(F)
    dP = calc_dp(P,F,dhval[key],dF)
    dF = calc_dF(P,F,dP)
    dFerr = np.sum(np.abs(np.subtract(dF,dF_old)))
    print("Error: %s " % dFerr)
    print("dF: %s" % np.sum(dF))
    isconverged=False
    if (dFerr < tolerance):
        isconverged=True
        np.savetxt(key+"_wham_pmf.out", np.c_[xvals,-kbT*(-kbT*np.log(P)-np.min(-kbT*np.log(P))+dP/P),dP])
        np.savetxt(key+"_wham_dF.out", np.c_[dF])
    return dF, isconverged



# START

#Command Line Arguments
parser=argparse.ArgumentParser()
parser.add_argument('-Nw', default=30,type=int, help='Number of Windows')
parser.add_argument('-rlow', default=1.5,type=float, help='Low cutoff in angstroms')
parser.add_argument('-rhigh', default=8,type=float, help='High cutoff in angstroms')
parser.add_argument('-nbin', default=100, type=int, help='Histogram bins')
parser.add_argument('-k', default=11.0, type=float, help='Force constant')
parser.add_argument('-plot', default=0, type=int, help='If 1, plots histograms in matplotlib')
parser.add_argument('-subfile', default="lif.distance", type=str, help="File name for colvar in subdirectory")
parser.add_argument('-subcol',default=1, type=int, help="Column of subfile for colvar")
parser.add_argument('-unit', default=1, type=int, help="[0] Angstrom [1] bohr (output always angstrom, kcal/mol)")
parser.add_argument('-skip', default=10000, type=int, help="How many elements to skip")
parser.add_argument('-opt', default=1, type=int, help="Functional form")
parser.add_argument('-deriv',default=0, type=int, help="[0] Turn off derivative calculation (default) [1] Turn on derivative calculation")
parser.add_argument('-enerfile', default="flucts.inp", type=str, help="File name with energy information, default flucts.inp")
parser.add_argument('-rerun', default=0, type=int, help="[0] read in files, [1] read pckl rerun file")
parser.add_argument('-T', default=298.15, type=float, help="T in units of Kelvin")
args = parser.parse_args()

Nwindows=args.Nw
rlow=args.rlow
rhi=args.rhigh
nbins=args.nbin
k=args.k
unit=args.unit
subfile=args.subfile
shouldplot=args.plot
skip=args.skip
subcol=args.subcol
opt=args.opt
shouldderiv=args.deriv
enerfile=args.enerfile
shouldrerun=args.rerun
T=args.T

# Constants
kb=0.0019872041
kbT=kb*T

# Checks if wham_pckl directory exists/makes the directory
if not os.path.exists('wham_pckl'):
    os.makedirs('wham_pckl')

# Check energy file
etypes=[]
if shouldderiv == 1:
    etypes=np.genfromtxt(enerfile,usecols=0,dtype=str)
    print("Derivative calculation turned on")
    print("Selected energy types:",*etypes)

calc_bias=choose_bias(opt)

# Set Units
BohrToAng=0.529177
convdist=1.0
convk=1.0
if unit==1:
    print("Converting r,k to units of angstroms, kcal/mol/ang^2 from bohr, hartree/bohr^2")
    convdist=BohrToAng
    convk=627.509/(BohrToAng**2.)

# Does unit conversions
k=k*convk
print("New low: %s New hi: %s" % (rlow, rhi))
print("new k: %s" % k)

# Sets up Arrays
U, ni, cnt, hval, center = [],[],[],[],[]
dhval={}
for key in etypes:
    dhval[key]=[]

xvals=np.linspace(rlow,rhi,num=nbins)
# Reads in bias potential locations
xc=np.genfromtxt("wham_metadata.info",usecols=1,unpack=True)
# Sets up the Histograms
svdata=[] # Saves the data for easy rerunning
for window in range(Nwindows):
    if window%10 == 0: print(window)
    data, en, den, dhist = [],{},{},{}
    # Reads Data In and Skips (if necessary) 
    # Also writes pickle files for ease. 
    if shouldrerun == 0:
        data=np.genfromtxt(str(window)+"/"+subfile, usecols=subcol,unpack=True)
        data=data*convdist # Converts data
        pickle.dump(data,open("wham_pckl/window_%d.pckl"%window,'wb'))
        data=data[skip:-1]
        if shouldderiv==1:
            for key in etypes:
                en[key]=np.genfromtxt(str(window)+"/"+key+"_init.out", usecols=0)
            pickle.dump(en,open("wham_pckl/ener_%d.pckl"%window,'wb'))
    else:
        data=pickle.load(open("wham_pckl/window_%d.pckl"%window,'rb'))
        data=data[skip:-1]
        if shouldderiv==1: 
            en=pickle.load(open("wham_pckl/ener_%d.pckl"%window,'rb'))
    hist,bins=np.histogram(data,bins=nbins,range=(rlow,rhi),density=False)
    # Derivative of biased distribution
    if shouldderiv ==1:
        for key in en:
            en[key]=en[key][skip:]
            den[key]=en[key]-np.average(en[key])
            dhist[key],dbins=np.histogram(data,bins=nbins,range=(rlow,rhi),density=False,weights=-den[key])
            dhval[key].append(dhist[key]/np.sum(hist)) # Array for easy plotting
    # Plotting Stuff
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    xvals=center
    if shouldplot==1: plt.bar(center, hist/np.sum(hist), align='center', width=width)
    # Bookkeeping
    ni.append(hist) # Array that stores counts for each x, window
    hval.append(hist/np.sum(hist)) # Array for outputting easy plotting
    cnt.append(np.sum(hist)) # Total number of counts for each window
    dx=np.subtract(xc[window],center) # difference between the bias location and each x
    U.append(calc_bias(dx,k)) # Array of the bias potential (2d array of shape dx, nwindow)


U=np.array(U)
for key in etypes:
    dhval[key]=np.array(dhval[key])

np.savetxt("hist_vals.dat", np.c_[center,np.array(hval).T])
for key in etypes:
    np.savetxt(key+"_hist_vals.dat", np.c_[center,np.array(dhval[key]).T])

# Plots the histograms using matplotlib
if shouldplot==1: plt.show()

# Basic Settings
isconverged = False
F_old=np.zeros(Nwindows)
tolerance = 1e-6
maxiter=100000
iteration=0
#does the primary wham iterations 
while ( isconverged == False):
    F_old,isconverged=wham_iteration(F_old)
    iteration+=1
    if iteration > maxiter:
        exit("Too many iterations")
isconverged = False
dF_old=np.ones(Nwindows)
iteration=0
key="TotEng"
if shouldderiv==1:
    # Does secondary wham iterations
    while ( isconverged == False):
        dF_old,isconverged=dwham_iteration(F_old,dF_old,key)
        iteration+=1
        if iteration > maxiter:
            exit("Too many iterations in derivative")
     



