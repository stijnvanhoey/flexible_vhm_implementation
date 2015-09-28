#-------------------------------------------------------------------------------
# Name:        VHM - python version 1
# Purpose:
#
# Author:      VHOEYS
#
# Created:     12/03/2011
# Copyright:   (c) VHOEYS 2011
# Licence:     BSD 3-Clause License
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import numpy as np
import random as rd
import matplotlib.pyplot as plt
from pylab import *

from scipy.integrate import odeint
from scipy.interpolate import interp1d

#%%
# HELP-FUNCTIONS
################################################################################

def normalize_pdf(inputarr):
    '''sum of all elements becomes 1
    '''
    normed = inputarr/inputarr.sum()
    return normed

def fracthandling(f_stor, f_ro, f_in, f_b, method='relative'):
    '''
    Handling of fractions if > 1.0

    The way this is done is unclear in the original submission of Willems '14,
    as such, we implemented different options to mimick the behaviour of
    known model outputs. Relative is probably the one used in the original
    implementation, but actually there isno specific reasoning behind it
    from a model-conceptual thinking to prefer this over the others (exept of
    using another model conceptualisation)
    '''
    f_stor = max(f_stor, 0.0)
    f_ro = max(f_ro, 0.0)
    f_in = max(f_in, 0.0)
    f_b = max(f_b, 0.0)

    if method == 'relative':
        fractions = np.array([f_stor, f_ro, f_in, f_b])
        input_cfr = np.zeros(len(fractions))
        fractions = np.maximum(input_cfr, fractions)
        fract_norm = normalize_pdf(fractions)
        return fract_norm

    if method == 'sequential1': #stor>ro>in>base
        if f_stor > 1.0:
            f_stor = np.float64(1.0)
            f_ro = 0.0
            f_in = 0.0
            f_b = 0.0
        if f_stor + f_ro > 1.0:
            f_ro = np.float64(1.0) - f_stor
            f_in = 0.0
            f_b = 0.0
        if f_stor + f_ro + f_in > 1.0:
            f_in = np.float64(1.0) - f_stor - f_ro
            f_b = 0.0
        f_b = 1.0 - f_stor - f_ro - f_in
        return array([f_stor, f_ro, f_in, f_b])

    if method=='sequential2':#stor>in>ro>base
        if f_stor>1.0:
            f_stor=np.float64(1.0)
            f_ro=0.0
            f_in=0.0
            f_b=0.0
        if f_stor+f_in>1.0:
            f_in=np.float64(1.0)-f_stor
            f_ro=0.0
            f_b=0.0
        if f_stor+f_ro+f_in>1.0:
            f_in=np.float64(1.0)-f_stor-f_ro
            f_b=0.0
        f_b=1.0-f_stor-f_ro-f_in
        return array([f_stor,f_ro,f_in,f_b])

    if method=='sequential3':#stor>base>ro>in
        if f_stor>1.0:
            f_stor=np.float64(1.0)
            f_ro=0.0
            f_in=0.0
            f_b=0.0
        if f_stor+f_b>1.0:
            f_b=np.float64(1.0)-f_stor
            f_ro=0.0
            f_in=0.0
        if f_stor+f_ro+f_b>1.0:
            f_ro=np.float64(1.0)-f_stor-f_b
            f_in=0.0
        f_in=1.0-f_stor-f_ro-f_b
        return array([f_stor,f_ro,f_in,f_b])

    if method=='sequential4':#stor>base>in>ro
        if f_stor>1.0:
            f_stor=np.float64(1.0)
            f_ro=0.0
            f_in=0.0
            f_b=0.0
        if f_stor+f_b>1.0:
            f_b=np.float64(1.0)-f_stor
            f_ro=0.0
            f_in=0.0
        if f_stor+f_in+f_b>1.0:
            f_in=np.float64(1.0)-f_stor-f_b
            f_ro=0.0
        f_ro=1.0-f_stor-f_in-f_b
        return array([f_stor,f_ro,f_in,f_b])

    if method=='rediff':
        fractions = np.array([f_stor, f_ro, f_in, f_b])
        diff = (fractions.sum()-1.)/2.
        f_ro = f_ro - diff
        f_in = f_in - diff
        return array([f_stor,f_ro,f_in,f_b])

def anteced_rain(data, degree, dropVals=False, visualise=False): #preprocessing of data
    '''
    performs sum of antecedet rainfall for degree timesteps (data-steps)
    * note that if dropVals is False, output length will be identical
    to input length, but with copies of data at the flanking regions
    '''
    smoothed=[]
    for i in range(degree,len(data)):
        point = data[i-degree:i]
        smoothed.append(sum(point)/degree*1.0)
    if visualise:
        pylab.plot(anteced_rain(data, 3), "-",
                   label="Antecedent rainfall with with window 3")
        pylab.plot(data, "k.-", label="Original data", alpha=.6)
        pylab.title("Antecedent rain")
        pylab.grid(alpha=.5)
        pylab.legend()
        pylab.savefig('Antecedentrain', dpi=600)

    if dropVals:
        return smoothed

    while len(smoothed) < len(data):
        smoothed.insert(0, sum(data[:len(data)-len(smoothed)-1])/degree*1.0)
    return smoothed

def linres(n_res ,q_init, co, k):  #nog VHM gewijs met cov te doen van vorige tijdstap
    if n_res==1:
        q_init[0]=q_init[0]*np.exp(-1/k) + co*(1 - np.exp(-1/k))
        return q_init
    else:
        q_init[n_res-1]=q_init[n_res-1]* np.exp(-1/k)+linres(n_res-1,q_init,co,k)[n_res-2]*(1 - np.exp(-1/k))
        return q_init

def linresv(n_res, q_init, co, v, k):  #nog VHM gewijs met cov te doen van vorige tijdstap
    if n_res==1:
        q_init[0]=q_init[0]*np.exp(-1/k) + co*(1 - np.exp(-1/k))*v
        return q_init
    else:
        q_init[n_res-1]=q_init[n_res-1]* np.exp(-1/k)+linres(n_res-1,q_init,co,k)[n_res-2]*(1 - np.exp(-1/k))*v
        return q_init

#%%

def VHM_flexible(pars, constants, init_conditions,
                    structure_options, rain, pet):
    """
    Flexible and straight forward application of the VHM modelling approach
    as proposed by Willems P. (2014) and mase flexible with a set of modelling
    options

    Parameters
    ------------
    pars : list
        18 parameters for the model listed
    constants : list
        area and timestep in a sinlge list
    init_conditions : list
        initial conditions for soil storage and different routing elements
    structure_options : list
        * fracthand 'relative' or 'sequentialx' with x [1-4]
        * storhand 'linear' or 'nonlinear'
        * interflowhand True or False
        * infexcesshand True or False
        * nres_g/nres_i/nres_o string of 3 options, each [1-2], eg 211, 121,...
    rain : ndarray
        array with input rain observations
    pet : ndarray
        array with input pet observations

    Returns
    ---------
    outflows : ndarray
        columns q_out, q_overland, q_interflow, q_baseflow
        total flow, overland, interflow and baseflow
    fractions : ndarray
        clumns fs, fi, fb, fu, i.e. fractions of overland flow, interflow,
        base flow and infiltration
    moisture : ndarray
        mooisture content in time
    """
    #Define the parameters
    ##################################
    #Storage model parameters
    umax = np.float64(pars[0])         #Maximum storage for soil moiture sompartiment
    uevap = np.float64(pars[1])        #soilmoisture content for  evapotranspiration
    c1s = np.float64(pars[2])
    c2s = np.float64(pars[3])
    c3s = np.float64(pars[4])
    #OF model parameters
    c1o = np.float64(pars[5])
    c2o = np.float64(pars[6])
    c3o = np.float64(pars[7])          #eventueel aanpassen, zodat enkel c1 of c3 in te geven, vermits zelfde parameter...
    c4o = np.float64(pars[8])
    #IF model parameters
    c1i = np.float64(pars[9])
    c2i = np.float64(pars[10])
    c3i = np.float64(pars[11])         #eventueel aanpassen, zodat enkel c1 of c3 in te geven, vermits zelfde parameter...
    c4i = np.float64(pars[12])
    ns_ro = int(pars[13])              #antecedent day = # voorgaande dagen waarvan bui no in rekening gebracht
    ns_in = int(pars[14])
    #Flow Routing parameter
    Kg = np.float64(pars[15])
    Ki = np.float64(pars[16])
    Ko = np.float64(pars[17])

    #Define the constantsants
    ##################################
    area = np.float64(constants[0])     #catchment area
    timestep = np.float64(constants[1])

    totn = rain.size     #Timesteps (eigenlijk niet nodig als gepast inputs!)

    #Define the options
    ##################################
    fracthand = structure_options[0]
    storhand = structure_options[1]
    interflowhand = structure_options[2]
    infexcesshand = structure_options[3]
    nres_o = int(structure_options[4][0])
    nres_i = int(structure_options[4][1])
    nres_g = int(structure_options[4][2])

    #Define the initial conditions
    ##################################
    u = np.float64(init_conditions[0])     #Soil moiosture storage
    qg = np.float64(init_conditions[1])*np.ones(nres_g)    #baseflow
    cg = np.float64(init_conditions[2])   #baseflow current timestep
    qi = np.float64(init_conditions[5])*np.ones(nres_i)    #interflow
    ci = np.float64(init_conditions[6])    #interflow current timestep
    qo = np.float64(init_conditions[3])*np.ones(nres_o)    #overland flow
    co = np.float64(init_conditions[4])    #overland flow current timestep
    fb = 0.

    #qr = qo + qi + qg
    v = np.float64(area * 1000.0 / (60.0 * 60.0))    #mm to m3/s

    moisture = np.zeros(totn, np.float64)
    outflows = np.zeros((totn, 4), np.float64)
    fractions = np.zeros((totn, 4), np.float64)

    #Define array for moving window -antecedent rainfall
    ##################################
    sb_ro = np.zeros(ns_ro+2)
    sb_in = np.zeros(ns_in+2)

    #Start dynamic run - explicit euler solution
    #############################################
    for t in range(totn):
        neersl = np.float64(rain[t])
        evap = np.float64(pet[t])

        #Soil Moisturs storage submodel
        if storhand == 'linear':                          #for relation between soil moisture storage and fraction interflow/overland
            fu = c1s - c2s*u/umax
        else:
            fu = c1s - np.exp(c2s*(u/umax)**c3s)    #correct?!?hetgeen in exp zit (c2*u/umax**c3 moet dus eigenlijk negatief zijn...niet dus

        #Evapotranspiration submodel  (Als input procentueel van totaal, dan aanpassen)
        if u > uevap:
            fe = evap
        else:
            fe = evap*u/(uevap)                     #lineair

        #ANTECEDENT rainfall OVERLAND to represent soil surface wetness
        sb_ro[ns_ro + 1] = neersl
        s_ro = np.float64(0.0)
        for ts in range(2,ns_ro + 2):
            sb_ro[ts-1] = sb_ro[ts]
            s_ro = s_ro + sb_ro[ts-1] / (ns_ro*1.0)

        #ANTECEDENT rainfall INTERFLOW to represent soil surface wetness
        sb_in[ns_in + 1] = neersl
        s_in = np.float64(0.0)
        for ts in range(2,ns_in + 2):
            sb_in[ts-1] = sb_in[ts]
            s_in = s_in + sb_in[ts-1] / (ns_in*1.0)

        # structure option infiltration excess handling with antecedent rain
        if not infexcesshand == True:
            s_ro = 1.
            s_in = 1.

        #Fractions overland - pre
        fs1 = np.exp(c1o + c2o*u/umax)      #c1o en c3o zijn te combineren tot 1 parameter!
        if s_ro > 0.0:
            fs2 = np.exp(c3o + c4o * np.log(s_ro))
        else:
            fs2 = np.float64(0.0)
        fs = fs1*fs2     #Fractions overland

        #fracties interflow - pre
        fi1 = np.exp(c1i + c2i*u/umax)      #c1i en c3i zijn te combineren tot 1 parameter!
        if s_in > 0.0:
            fi2 = np.exp(c3i + c4i * np.log(s_in))
        else:
            fi2 = np.float64(0.0)
        fi = fi1*fi2     #Fractions interflow

        # structure handling with the interflow incorporation
        if not interflowhand == True:
            fi = 0.

        fb = 1. - fu - fs - fi
        #print fu, fs, fi, fb, u, s_in, s_ro

        #fractionhandling - not the best part of the model :-(
        update_fract = fracthandling(fu, fs, fi, fb, method=fracthand)  # return array([f_stor,f_ro,f_in,f_b])
        fu=update_fract[0]
        fs=update_fract[1]
        fi=update_fract[2]
        fb=update_fract[3]

        #Soil moisture balance
        u = u + (fu*neersl - fe)*timestep   # Soil storage + neerslag - evapotranspiratie
        co = fs * neersl                    # overland flow
        ci = fi * neersl                    # interflow
        cg = fb * neersl

        # Routing
        # Linear reservoir for the  3 routing types
        qo = linresv(nres_o, qo, co, v, Ko)
        qi = linresv(nres_i, qi, ci, v, Ki)
        qg = linresv(nres_g, qg, cg, v, Kg)

        q = qo[nres_o-1] + qi[nres_i-1] + qg[nres_g-1]

        # extra outputs
        outflows[t, :] = np.array([q, qo[nres_o-1], qi[nres_i-1], qg[nres_g-1]])
        fractions[t, :] = np.array([fs, fi, fb, fu])
        moisture[t] = u

    return outflows, fractions, moisture