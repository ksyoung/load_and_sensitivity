# -*- coding: utf-8 -*-
'''
Cleaned up mapping speed and loading codes.  Original was from Ben Westbrook.
This is the stripped down version to do only what I need, but no more.



This file is functions to calculate load per element, spillover, ???,
etc.

Called by calc_mapping_speed.py

'''


import os, sys
import matplotlib
#import pdb
import pandas
import time
import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt
from pprint import pprint
from copy import copy, deepcopy


### needed constants
k_b = 1.38064852e-23 # J / K
c   = 299792458 # m/sec
one_over_c = 1./c 
h_  = 6.6260699e-34  # J*s
Tcmb = 2.725
###

def planck_func(nu,T):  
  
  val =  2 * h_ * nu * nu * nu * one_over_c * one_over_c / \
         (np.exp( (h_ * nu) / (k_b * T)) - 1) 
  return val

def planck_func_lambda_sq(nu,T):
  # really is power from an element.  
  return 2 * h_ * nu / (np.exp((h_ * nu) / (k_b * T)) - 1)

def planck_grey(nu,T,beta,nu_0):
  print 'fail'
  sys.exit()
  return

def power_emitted_per_element(T,emiss,band, through='lambda_sq'):
  '''
  calculates in band power emitted by a given element.
  inputs:  temperature, emissivity, band (2 element array, tophat assumed)
  outputs: power(W)

  future, add ability for any band.
  will need to use different integration method.

  Assuming throughput = lamda**2 as norm.
  for 100% of light, i.e. both polarizations included. (would be 1/2 for a single pol)
  '''
  if through is 'lambda_sq':
    if len(band) is 2:
      power, error = scint.quad(planck_func_lambda_sq,band[0],band[1],args=(T))
      #power, error = scint.quad(planck_func,band[0],band[1],args=(T))
      power = power*emiss
    else:
      print '\nTophat bands only please!\n'
      power, error = [np.nan,np.nan]
  else:
    print '\nCode isn\'t smart enough for your fancy throughputs. Rewrite needed'
    power, error = [np.nan,np.nan]

  return power, error

def losstan_to_absorption(freq,loss_tan,index,thick):
  # assume 1 freq at band center. usually a fair-ish average.
  # all SI units.
  # abs = emiss
  pass
  return 1-np.e**(-thick*loss_tan*2.0*np.pi*index*freq/c)

def calc_mirror_transmission(emiss150, band, metal='alum'):
  #assume aluminum 
  if alum:
    pass
    # write some code!
  else:
    print 'I don\'t know that metal'
    sys.exit()

  return

#def calc_lens_transmission(emiss (or loss tan), thick, band):
#  return


def calc_beam(D_m,edge_db,wavelength):
  # approximate Karl found.  Need to track down source.  From radio astro peeps.
  FWHM_rad = (1.02 + 0.0135*edge_db) * wavelength/(D_m)   #in radians 
  FWHM =  FWHM_rad*(180.0/np.pi)*60 #convert to arcminutes
  return FWHM

def gauss(x,sigma_sq):
  # feed me sigma squared!!!!!!!!!!!!
  return np.exp(-(x**2) /(2*sigma_sq))

def efield_integrand(theta, sigma_sq):
  # integrand the e-field squared (gaussian assumed) for spillover calc.
  # from Ben's code and http://gmrt.ncra.tifr.res.in/gmrt_hpage/Users/doc/WEBLF/LFRA/node172.html
  return (gauss(theta,sigma_sq) ** 2) * np.sin(theta)

def efield_integrand_illum_top(theta, sigma_sq):
  return gauss(theta,sigma_sq) * np.tan(theta * 0.5)

def get_spill_illum_effs(edge_db, f_num): 
  # integrate e-field over mirror and over all space, divide to get spill.
  # compare to flat to get illumination.
  edge_power = 10**(edge_db*-.1)
  theta_edge = np.arctan(1./(2.*f_num))
  sigma_sq = theta_edge**2 / (- np.log(edge_power))
  mirror_field_sq = scint.quad(efield_integrand,0,theta_edge,args=(sigma_sq))[0]
  #print edge_power,theta_edge, sigma_sq,mirror_field_sq
  spill = mirror_field_sq / scint.quad(efield_integrand,0,np.pi,args=(sigma_sq))[0]
  illum = 2. * (1. / np.tan(theta_edge*0.5)**2) * scint.quad(efield_integrand_illum_top \
          ,0,theta_edge,args=(sigma_sq))[0]**2 / mirror_field_sq

  return spill, illum

def get_spill_illum_effs_from_theta_px(theta_px, f_num): 
  # integrate e-field over mirror and over all space, divide to get spill.
  # compare to flat to get illumination.
  theta_edge = np.arctan(1./(2.*f_num))
  sigma_sq = 0 # ??? need to convert theta_px into sigma of a gaussian.  In electric field.
  mirror_field_sq = scint.quad(efield_integrand,0,theta_edge,args=(sigma_sq))[0]
  #print edge_power,theta_edge, sigma_sq,mirror_field_sq
  spill = mirror_field_sq / scint.quad(efield_integrand,0,np.pi,args=(sigma_sq))[0]
  illum = 2. * (1. / np.tan(theta_edge*0.5)**2) * scint.quad(efield_integrand_illum_top \
          ,0,theta_edge,args=(sigma_sq))[0]**2 / mirror_field_sq
  print 'NOT WORKING YET!!'
  sys.exit()
  return spill, illum

def poisson_integrand_lambda_sq(nu,T):
  # an extra 2 * h_ * nu over power integrand
  return (2. * h_ * nu) * 2 * h_ * nu / (np.exp((h_ * nu) / (k_b * T)) - 1) 

def bunch_integrand_lambda_sq(nu,T):
  # power integrand ** 2.   
  return (2. * h_ * nu / (np.exp((h_ * nu) / (k_b * T)) - 1))**2.

def bunch_integrand_lambda_sq_multiple_surfaces(nu,T_arr,emiss_arr,cum_eff_arr):
  # integrand for all surfaces simultaneously.
  # still assuming lambda_sq
  #  can do as array! handy.  needs np.array inputs.
  # returns (NEP_bunch)**2
  return (2. * h_ * nu * np.sum(emiss_arr*cum_eff_arr/(np.exp((h_ * nu) / (k_b * T_arr)) - 1)))**2.


def NEP_photon_per_element(T,emiss,band,cum_eff, photon_corr=1,through='lambda_sq'):
  '''
  equations /method from Richards bolo review (and maybe toki thesis?)
  another possible good paper. https://www.osapublishing.org/ao/fulltext.cfm?uri=ao-42-25-4989&id=74102  but a lot of quantum details to follow.

  Inputs:
   blackbody (T) + emiss  ## will have to handle mirror specially. h hmm..mm....
   band
   cum_eff
   Assume correlation = 1
  returns: 
   total noise          ##  All in units NEP**2
   1st term (poisson)
   2nd term (bunching)
  '''
  if through is 'lambda_sq':
    if len(band) is 2:
      nep_poisson, p_error = scint.quad(poisson_integrand_lambda_sq,band[0],band[1],args=(T))
      nep_poisson = nep_poisson * emiss * cum_eff

      nep_bunch, b_error = scint.quad(bunch_integrand_lambda_sq,band[0],band[1],args=(T))
      nep_bunch = nep_bunch * (emiss * cum_eff)**2. # 2 factors of emiss and cum_eff.
    else:
      print '\nTophat bands only please!\n'
      return
  else:
    print '\nCode isn\'t smart enough for your fancy throughputs. Rewrite needed'
    return

  return nep_poisson, p_error, nep_bunch, b_error


def bolo_properties(p_opt,t_bath,safety_factor,n,bolo_Rn,bias_point):
  if n == 1.0:
    t_c = 2.134 * t_bath  # For a metal (Au) link
  elif n == 3.0:
    t_c = 1.705 * t_bath  # For a SiN link
  else:
    t_c = np.interp(n,[1.,3.],[2.134,1.705])

  p_sat = safety_factor * p_opt
  p_elec = p_sat - p_opt
  g_bar = p_sat / (t_c - t_bath)
  G_dyn = g_bar * (n+1) * ((1-t_bath/t_c) / ( 1 - (t_bath/t_c)**(n+1) ))
  v_bias = ((bolo_Rn*bias_point)*p_elec)**(0.5)

  gamma = ((n + 1) / (2*n + 3)) * ((1 - (t_bath / t_c) ** (2*n + 3)) / (1 - (t_bath / t_c) ** (n + 1)))

  return p_sat, g_bar, G_dyn, t_c, v_bias, gamma

def calc_bolo_noises(G_dyn,t_c, R_bolo_biased, v_bias, gamma, readout_noise_amps):
  # R_bolo_biased is resistance at bias.  ??? need to check with franky
  # readout_noise_amps is number from franky of readout in A*sqrt(s). need to check value.
  # rest are obvious? i hope...
  phonon = np.sqrt(4 * gamma * k_b * t_c ** 2 * G_dyn)
  johnson = np.sqrt(4 * k_b * t_c / R_bolo_biased) * v_bias
  readout = readout_noise_amps*v_bias/np.sqrt(2)

  return phonon,johnson,readout

def dPdT_integrand_lambda_sq(nu,T):
  return ( (1 / (np.e ** ((h_ * nu) / (k_b * T)) - 1))**2. * (nu / T)**2. *np.e**((h_ * nu) / (k_b * T)) )

def nep_to_net_Kcmb(nep,cum_eff,band):
  # convert given nep to net_Kcmb.  (from Toki's thesis, added 2 is becuase he ignores 2 pols.)
  # need band.
  # lambda squared still assumed.
  # integrate dP/dT, const and cum_eff in front of int.
  print  scint.quad(dPdT_integrand_lambda_sq,band[0],band[1],args=(Tcmb))
  dPdT, error = scint.quad(dPdT_integrand_lambda_sq,band[0],band[1],args=(Tcmb))
  dPdT = dPdT * (2* h_ ** 2 / k_b) * cum_eff  # the coeffs in front of integral
  return nep/(np.sqrt(2)*dPdT)
