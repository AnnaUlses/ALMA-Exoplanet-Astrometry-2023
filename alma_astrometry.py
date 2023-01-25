#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 11:06:06 2022

@author: annagraceulses
"""
#this file calculates data to be used for plotting and interpretation, but if you want to run it you have to download and import all the csv files in
#this repository

import pandas as pd
import numpy as np
from astropy.io.votable import parse
import os

wd = os.chdir('/Users/annagraceulses/Library/Mobile Documents/com~apple~CloudDocs/2022-2023/FYP/Data and Scripts')
photometry = pd.read_csv('brighter3mag_photometry.csv')
filterss = photometry.filterr.unique()[140:176] #alma filters 

alma_phot = photometry.loc[(photometry.filterr == 'WAV880')]
                           
alma_phot = alma_phot.reset_index(drop = True)
alma_phot.to_csv('alma_phot.csv', index = False)

alma_fwhm = pd.read_csv('alma_fwhm.csv')
band_freq = [100,150,185,230,345,460,650,870] #in GHz
alma_fwhm.insert(1, 'freq', band_freq)

#convert vot to dataframe for ease
def votable_to_pandas(votable_file):
    votable = parse(votable_file)
    table = votable.get_first_table().to_table(use_names_over_ids=True)
    return table.to_pandas()

fitted_pars = votable_to_pandas("brighter3mag.xml")
observables = fitted_pars.loc[(fitted_pars.dej2000 < 37) & (fitted_pars.dej2000 > -50)] #stars observable\
    #by ALMA based on dec
observables = observables.copy()
observables.e_a_v = 1 / observables.plx_arcsec
observables.rename(columns = {'e_a_v' : 'distance'}, inplace = True) #add distance to df
observables.to_csv('observable_sample.csv')

hip2_223 = pd.read_csv('hip2_223.csv')
hip2_223['HIP'] = 'HIP_' + hip2_223['HIP'].astype(str)
hip2_notobs = hip2_223.drop(observables.index)
hip2_observable = hip2_223.drop(hip2_notobs.index)
hip2_observable = hip2_observable.reset_index(drop = True)
HIP_id = hip2_observable.HIP #id of objects left in the sample
alma_phot = alma_phot.drop(hip2_notobs.index)
observables = observables.reset_index(drop = True)

decs_n = np.linspace(0,90, 19)
decs_s = np.linspace(0,-90, 19)
observable_sample = []
for n,s in zip(decs_n, decs_s):
    observable = fitted_pars.loc[(fitted_pars.dej2000 < n) & (fitted_pars.dej2000 > s)]
    observable_sample.append(len(observable))

def wave_to_freq(wavelength):
    freq = 3e8 / wavelength / 1000
    return freq #returns in GHz

def freq_to_wave(frequency):
    frequency *= 1e9
    wave = 3e8 / frequency #frequency is in Hz
    wave *= 1000
    return wave #returns in microns

band_wave = []
for w in band_freq:
    wv = freq_to_wave(w)
    band_wave.append(wv)

def extrapolate_flux(band_freq):
    flux_2 = alma_phot.star_jy * (band_freq**2) / (wave_to_freq(880)**2)
    return flux_2 #returns jy

fluxs = pd.DataFrame()
for b, x in zip(band_freq, range(len(band_freq))):
    f = extrapolate_flux(b)
    fluxs.insert(0 + x, 'band'+str(x+3), f)

alma_sensitivity = pd.read_csv('sensitivity_alma.csv') 
fluxs = fluxs.T
fluxs.columns = HIP_id

def snr(signal, noise): #both are in mJy
    return (signal * 1000) / noise 

sens_15 = alma_sensitivity.sens_15 #these are calculated from the ALMA sens. calculator
sens_60 = alma_sensitivity.sens_60

def snr_id(hipid):
    snrs = []
    for c,f in zip(fluxs[hipid], sens_15): #switch to 60 to test for an hour 
        snrs.append(snr(c,f))
    return snrs

snrr = pd.DataFrame()
for h,x in zip(HIP_id, range(len(hip2_observable.HIP))):
    s = snr_id(h)
    snrr = snrr.copy()
    snrr.insert(0+x, ''+str(h), s)

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, 'ALMA band snr plots/')

if not os.path.isdir(results_dir): #for saving the figures 
    os.makedirs(results_dir)    


def astrometric_acc(snr, fwhm): #fwhm = alma_fwhm.fwhm
    sigma = fwhm / (snr/2) 
    return sigma #in arcsec

def astro_acc_id(hip_id):
    accs = []
    for e,l in zip(snrr[hip_id], alma_fwhm.fwhm):
        accs.append(astrometric_acc(e,l))
    return accs

ast_acc = pd.DataFrame()
for a,g in zip(HIP_id, range(len(hip2_observable.HIP))):
    ast = astro_acc_id(a)
    ast_acc = ast_acc.copy()
    ast_acc.insert(0+g, ''+str(a), ast)
    
min_ast = ast_acc.min() #4 corresponds to band 7
max_snr = snrr.idxmax() #also band 7

print('Ideal ALMA observation band for this sample is band 7')

z = np.linspace(0,len(observables.id),len(observables.id))

ast_acc.index = alma_fwhm.Band

script_dir1 = os.path.dirname(__file__)
results_dir1 = os.path.join(script_dir1, 'Misle. Plots/')

if not os.path.isdir(results_dir1): #for saving the figures 
    os.makedirs(results_dir1)

def angular_radius_star(diameter, distance):
    return 206265 * (diameter/distance)

stellar_diameter = (observables.rstar * 696340) * 2 #convert to km 
stellar_dist = (1 / observables.plx_arcsec) * 3.1e13 #distance in km
resolution = pd.DataFrame({'HIP_id': HIP_id, 'ang_size':angular_radius_star(stellar_diameter, stellar_dist)})

unresolved = resolution.loc[(resolution.ang_size < alma_fwhm.fwhm[4])] #unresolved stars in Band 7
resolved = resolution.loc[(resolution.ang_size >= alma_fwhm.fwhm[4])]

#determine distance limit for sample
observables.insert(60, 'ang_size_band7', resolution.ang_size)
resolved_full = observables.loc[(observables.ang_size_band7 >= alma_fwhm.fwhm[4])]
distance_limit = 1 / max(resolved_full.plx_arcsec)
unresolved_full = observables.loc[(observables.ang_size_band7 < alma_fwhm.fwhm[4])]

#take angular stellar a to be 3 * ast_acc to begin with
#can get mass ratio and planet a from that 
#this is the detectable stellar semi-major axis based off the positional uncertainties\
    #associated with each particular star in the sample using ALMA Band 7
ast_acc_au = ast_acc.iloc[4].values * (1 / observables.plx_arcsec) # astrometric accuracy is in AU now 
stellar_a = 3 * (ast_acc_au) #in AU equal to alpha 


band7_ast = pd.DataFrame({'id': HIP_id,'band7':ast_acc.iloc[4].values * 1000, 'radius': observables.rstar, 'lum': observables.lstar, 'T': observables.teff}) #ast acc is in mas here
#sensitivity curves 

a_p = np.linspace(0.1,10, 50) #range of planet semi-major axis in AU

def mass_ratio(sigma3, ap): #3sigma is stellar_a (minimum detectable star semi-major axis)
    return sigma3 / ap

#separate ms from non-ms stars
main_seq = observables.loc[observables['sp_type'].str.contains('V', case = True)]
non_ms = main_seq.loc[main_seq['sp_type'].str.contains('I', case = True)]
main_seq = main_seq.drop(non_ms.index)
main_seq_sp = main_seq.sp_type.to_csv('main seq spec type.csv')
main_seq_mass = pd.read_csv('main seq spec type.csv')
main_seq = main_seq.drop([7, 76, 101, 126])

#convert stellar radius to AU to check if potential planets would be within the star 
radius_au = observables.rstar * 696340 / 1.5e8 # in AU 
#stars with radius larger than 1 AU

non_ms = observables.drop(main_seq.index)

c = [] #coefficient for calculating non-ms stellar mass
for g in observables.logg: 
    if g > 4:
        c.append(1)
    if 3 < g < 4:
        c.append(1.2)
    if 2 < g < 3:
        c.append(1.4)
    if 1.6 < g < 2:
        c.append(2)
    if 0.9 < g < 1.6:
        c.append(3)
    if 0 < g < 0.9:
        c.append(4)
    if g <= 0:
        c.append(5)

def non_ms_mass(c, t):
    mass = c * (t/5770)**2
    return mass #in solar mass units 

def ms_mass (lstar): #calculate mass using mass-luminosity relation gives completely different masses
    return 10 ** (1/3.5 * np.log10(lstar))

msmass = []
for mm in main_seq.lstar:
    ma = ms_mass(mm)
    msmass.append(ma)
main_seq.insert(3, 'mass', msmass, True)

stellar_mass = []
for cc,tt in zip(c, observables.teff): #fix this so it only calculates the mass for the non-ms stars
    masss = non_ms_mass(cc,tt)
    stellar_mass.append(masss)
    
observables.insert(3, 'mass', stellar_mass, True)

def mass_cutoff(mass):
    return (3 **(2/3)) * ((mass) ** (1/3)) #sets a limit on a_p according to observation time 
# 3 here is delta T (observation time)

cutoff = []
for v in stellar_mass:
    cut = mass_cutoff(v)
    cutoff.append(cut)

#make different a_p arrays based on the left and right cutoff 
ap_arrays = [] #make the observable a_p based off the minimum and maximum (R* and deltaT)
for aa in range(len(radius_au)):
    ar = np.linspace(radius_au[aa], cutoff[aa], 149)
    ap_arrays.append(ar)

#%%
def planet_mass(ratio, starmass):
    return ratio * starmass #solar units

mass_ratios = [] #constrain the lines by stellar mass and radius
#re-calculate the Mp/M* that has the lower and upper limits 
for i,j in zip(stellar_a, ap_arrays):
    m = mass_ratio(i,j)
    mass_ratios.append(m)

mp = [] #og planet mass for uncorrelated to relationiship between period and observation time
for k,l in zip(mass_ratios, observables.mass):
    t = planet_mass(k,l) * 1047.57 #put into M_J unit
    mp.append(t)

#dont actually need orbital period values for following calculations, because boundary conditions are put in terms of \
    #a planet and observation time

#%%ADDING FURTHER COMPLEXITY TO SENSITIVITY CURVES 
deltaT = 3
ap_long = [] #planet semi-major axes in long period regime 
for ee in range(len(observables.mass)):
    minn = (((4 * deltaT) / 3) ** (2/3)) * (observables.mass[ee] ** (1/3))
    long = np.linspace(minn, 10, 149) #minimum set by observation time and period,\
        #max set by original ap_array
    ap_long.append(long)

ap_mid = [] #planet semi-major axes in mid-range period regime 
for ff in range(len(observables.mass)):
    max_mid = (((4 * deltaT) /3) ** (2/3)) * (observables.mass[ff] ** (1/3))
    min_mid = (deltaT ** (2/3)) * (observables.mass[ff] ** (1/3))
    mid = np.linspace(min_mid, max_mid, 149)
    ap_mid.append(mid)

ap_short = [] #planet semi-major axes in short period regime 
for gg in range(len(observables.mass)):
    max_short = (deltaT ** (2/3)) * (observables.mass[gg] ** (1/3))
    minn_ap = np.linspace(min(ap_arrays[gg]), max_short, 149)
    ap_short.append(minn_ap)

#sigma is stellar_a, i.e 3* ast_acc_au, or A1s in the Eisner paper
def long_p(a_planet, sigma, m_s):  #mp in the long period regime 
    return ((2 * np.sqrt(2) * sigma) / (1 - np.cos((np.pi* deltaT) * (a_planet**(-3/2) * m_s**(1/2))))) * (m_s/a_planet)
   
def mid_p(a_planet, sigma, m_s): #mp in the mid range period regime 
    return sigma * (m_s / a_planet) #remove the sqrt(2) to make the lines join for p>deltat

def short_p(a_planet, sigma, m_s): #mp in the short range period regime, same as before
    return sigma * (m_s/ a_planet)

#%% calculating the mp arrays 
long_mp = []
for yy,uu,cc in zip(ap_long, stellar_a, observables.mass):
    tt = long_p(yy,uu,cc)  * 1047.57
    long_mp.append(tt)

mid_mp = []
for ii, jj, xx in zip(ap_mid, stellar_a, observables.mass):
    plant = mid_p(ii,jj,xx)  * 1047.57
    mid_mp.append(plant)
    
short_mp = []
for aa,bb,dd in zip(ap_short, stellar_a, observables.mass):
    jeepers = short_p(aa,bb,dd)  * 1047.57
    short_mp.append(jeepers)

#factor difference should be always 4/3 discontinuity
#%%
print('Factor difference from mid_mp to long_mp:' ,long_mp[0][0] / mid_mp[0][-1])
print('Factor difference from short_mp to mid_mp:' ,mid_mp[0][0] / short_mp[0][-1])
print('Factor difference from short_mp to long_mp: ', long_mp[0][0] / short_mp[0][-1])

#%%
min_mass = []
min_ap = []
max_mass = []
max_ap = []
for vv in range(len(mid_mp)):
    minn = min(mid_mp[vv])
    min_mass.append(minn)
    maxx = max(short_mp[vv])
    max_mass.append(maxx)
    ap = max(ap_mid[vv])
    min_ap.append(ap)
    maxx_ap = max(ap_long[vv])
    max_ap.append(maxx_ap)





