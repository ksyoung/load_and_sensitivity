# -*- coding: utf-8 -*-
'''
Cleaned up mapping speed and loading codes.  Original was from Ben Westbrook.
This is the stripped down version to do only what I need, but no more.



This file is functions to calculate load per element, spillover, noise, bolo properties, 
 ..., etc.

Called by calc_mapping_speed.py

'''


import os, sys
import matplotlib
import pdb
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
  # approximation Karl found.  in Goldsmith's Quasioptical systmes p 139 (Shaul's copy.)  From radio astro peeps.
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
  sigma_sq = theta_edge**2 / (- np.log(edge_power))  ## converting to sigma of the guassian beam. E_field = exp[-x**2/(2*sigma**2)]
  mirror_field_sq = scint.quad(efield_integrand,0,theta_edge,args=(sigma_sq))[0]
  #print edge_power,theta_edge, sigma_sq,mirror_field_sq
  spill = mirror_field_sq / scint.quad(efield_integrand,0,np.pi,args=(sigma_sq))[0]
  illum = 2. * (1. / np.tan(theta_edge*0.5)**2) * scint.quad(efield_integrand_illum_top \
          ,0,theta_edge,args=(sigma_sq))[0]**2 / mirror_field_sq

  return spill, illum

def get_theta_px(wavelength, D_px,waist_ratio):
  # relation from Toki's thesis for SAMPs.  Still holds?
  return waist_ratio*wavelength / (np.pi * D_px)

def get_D_px(wavelength,theta_px,waist_ratio):
  # inverse of relation from Toki.
  return waist_ratio*wavelength / (np.pi * theta_px)

def get_spill_illum_effs_from_theta_px(theta_px, f_num): 
  # integrate e-field over mirror and over all space, divide to get spill.
  # compare to flat to get illumination.
  theta_edge = np.arctan(1./(2.*f_num))
  sigma_sq = 0.5 * theta_px**2.  # this is sigam of E-field. assuming theta_px is where E-field is down by 1/e. or intensity is down by 1/e**2.
  mirror_field_sq = scint.quad(efield_integrand,0,theta_edge,args=(sigma_sq))[0]
  #print edge_power,theta_edge, sigma_sq,mirror_field_sq
  spill = mirror_field_sq / scint.quad(efield_integrand,0,np.pi,args=(sigma_sq))[0]
  illum = 2. * (1. / np.tan(theta_edge*0.5)**2) * scint.quad(efield_integrand_illum_top \
          ,0,theta_edge,args=(sigma_sq))[0]**2 / mirror_field_sq
  return spill, illum

def theta_px_from_edge_db(edge_db,f_num):
  theta_edge = np.arctan(1./(2.*f_num))
  edge_power = 10**(edge_db*-.1)
  theta_px = theta_edge * np.sqrt(2. / -np.log(edge_power))
  return theta_px

def edge_db_from_theta_px(theta_px,f_num):
  theta_edge = np.arctan(1./(2.*f_num))
  edge_power = np.exp(-2. * (theta_edge / theta_px)**2.)
  edge_db = -10 * np.log10(edge_power) # I'm using db below as positive.  just laziness.
  return edge_db

def scale_db_between_bands(db_in, nu_in, nu_out):
  # scale db from one band to another assuming constant pixel size
  # and theta \propto lambda.
  # assume fixed f_num. ratio is always const since I use f_num to go to and from theta_px
  f_num  = 1.
  theta_px1 = theta_px_from_edge_db(db_in,f_num)
  theta_px2 = theta_px1 * nu_in / nu_out
  db_out = edge_db_from_theta_px(theta_px2,f_num)
  return db_out

def poisson_integrand_lambda_sq(nu,T):
  # an extra 2 * h_ * nu over power integrand
  return (2. * h_ * nu) * 2 * h_ * nu / (np.exp((h_ * nu) / (k_b * T)) - 1) 

def bunch_integrand_lambda_sq(nu,T):
  # power integrand ** 2.   
  # factor of 2 in front from Roger O'Brient
  return 2*(2. * h_ * nu / (np.exp((h_ * nu) / (k_b * T)) - 1))**2.

def bunch_integrand_lambda_sq_multiple_surfaces(nu,T_arr,emiss_arr,cum_eff_arr):
  # integrand for all surfaces simultaneously.
  # still assuming lambda_sq
  #  can do as array! handy.  needs np.array inputs.
  # returns (NEP_bunch)**2
  # factor of 2 in front from Roger O'Brient
  return 2*(2. * h_ * nu * np.sum(emiss_arr*cum_eff_arr/(np.exp((h_ * nu) / (k_b * T_arr)) - 1)))**2.


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
  '''
  old method of getting t_c, from Toki's thesis. below in Franky's, actually doing the optimization each time.
  if n == 1.0:
    t_c = 2.134 * t_bath  # For a metal (Au) link
  elif n == 3.0:
    t_c = 1.705 * t_bath  # For a SiN link
  else:
    t_c = np.interp(n,[1.,3.],[2.134,1.705]) * t_bath
  '''
  ## new t_c method from Franky.  acutally find minimum of function.
  dx = 0.01
  x = np.arange(1.+dx, 10., dx)
  y = (x**(2*n+3) - 1) / (x**(n+1) - 1)**2
  min_index = np.argmin(y)
  min_x = x[min_index]
  t_c = min_x * t_bath

  p_sat = safety_factor * p_opt
  p_elec = p_sat - p_opt
  g_bar = p_sat / (t_c - t_bath)
  G_dyn = g_bar * (n+1) * ((1-t_bath/t_c) / ( 1 - (t_bath/t_c)**(n+1) ))
  v_bias = ((bolo_Rn*bias_point)*p_elec)**(0.5)

  gamma = ((n + 1.) / (2.*n + 3.)) * ((1. - (t_bath / t_c) ** (2.*n + 3.)) / (1. - (t_bath / t_c) ** (n + 1.)))

  return p_sat, g_bar, G_dyn, t_c, v_bias, gamma

def calc_bolo_noises(G_dyn,t_c, R_bolo_biased, v_bias, gamma, readout_noise_amps):
  # R_bolo_biased is resistance at bias.  
  # readout_noise_amps is number from franky of readout in A*sqrt(s). need to check value. This now replaced by franky's code in final calculation.  Here for reference.
  # rest are obvious? i hope...
  phonon = np.sqrt(4 * gamma * k_b * t_c ** 2 * G_dyn)
  johnson = np.sqrt(4 * k_b * t_c / R_bolo_biased) * v_bias
  readout = readout_noise_amps*v_bias/np.sqrt(2) # should call Franky's readout noise calc function.

  return phonon,johnson,readout

def TDM_time_const(Ldc, C0, G_dyn, Tc, T0, R_op, Rs_Rop, L_sq_nyq, bI=0.,g=1.):
  # from Roger Obrient's code.  calc all applicable taus for TDM.
  # inputs:
  # Ldc - DC loop gain.
  # C0 - heat capacity, J/K
  # G_dyn - dynamic thermal conductance
  # g  - ??????????????????????????  Is 1 in RO's code though.
  # Tc - transition temp.  Kelvin
  # T0 - bath temp   K
  # R_op, bolo R at operation     Ohm
  # Rs_op, shut resistor R at operation   Ohm
  # L_sq_nyq - Squid nyquist inductor, 2 uH in Roger's inputs.
  # bI - describes change in TES resistance vs current, normally 0.  (from Roger's presentation.)
  # output dictionary of taus.
  tau = {}
  
  # untested.
  # all taus in seconds.
  tau['o'] = (C0*(Tc/T0)**g)/G_dyn; #s, time const, no feedback
  tau['e'] = L_sq_nyq/(Rs_Rop*R_op+R_op*(1+bI))
  tau['LG'] = C0*(Tc/T0)**g/(G_dyn*(1-Ldc))
  tau['minus'] = 1./((0.5/tau['e'])+(0.5/tau['LG'])-0.5*np.sqrt((((1./tau['e'])-(1./tau['LG']))**2)-4*R_op*Ldc*(2+bI)/(L_sq_nyq*tau['o'])))
  tau['plus'] = 1./((0.5/tau['e'])+(0.5/tau['LG'])+0.5*np.sqrt((((1./tau['e'])-(1./tau['LG']))**2)-4*R_op*Ldc*(2+bI)/(L_sq_nyq*tau['o'])))
  tau['mux_stab_fac'] = ((1./tau['e']-1./tau['LG'])**2)/(4*R_op*Ldc*(2+bI)/(Ldc*tau['o']))

  return tau

def TDM_readout_noise(NEP_photon, NEP_phonon, P_inc, Tc, T0, G_dyn, alpha, n, R_op, Rs_Rop, safety_factor, C0, L_sq_nyq,bI=0.):
  # taken from Roger Obrient's code. should be verbatim except for name changes to match my conventions.
  # inputs are:
  # NEP_photon
  # NEP_phonon
  # P_inc - load power   Watts
  # Tc - transition temp.  Kelvin
  # T0 - bath temp   K
  # G_dyn - same def. as earlier.   W/K
  # alpha - slope of transition, typically ~100
  # n, standar thermal exponent.  usually 2-3.
  # R_op, bolo R at operation     Ohm
  # Rs_op, shut resistor R at operation   Ohm
  # safety_factor - ratio of Psat to Pload
  # C0
  # R_dyn_SQ - ????????????????
  # Lstray - ????????????
  # In_SQ - ?????????????????
  # bI - describes change in TES resistance vs current, normally 0.  (from Roger's presentation.)
  
  # outputs johnson, shunt, alias noise.  Everything else (photon, phonon) is same for TDM/FDM.

  # untested
  P_bias = P_inc*(safety_factor-1.); # headroom, equivalently assumed bias power.   Watt
  I0=np.sqrt(P_bias/R_op); #Current through bolo @ operation   Amps
  Ldc=alpha*P_bias/(G_dyn*Tc); #Loop gain, DC bias
  # pdb.set_trace()
  # uncorrected
  johnson = 5.*np.sqrt(4.*k_b*Tc*R_op)*I0/Ldc;  #W/rtHz   ????????????????? why a five???
  shunt = I0*((Ldc-1)/Ldc)*np.sqrt(4*k_b*T0*Rs_Rop*R_op);  #W/rtHz   ??????? where is (1 + 2pi f tau) term?

  #get taus
  tau = TDM_time_const(Ldc, C0, G_dyn, Tc, T0, R_op, Rs_Rop, L_sq_nyq, bI=0.)

  S_DC=-(1./(I0*R_op*(2+bI)))*(1-(tau['plus']/tau['LG']))*(1-tau['minus']/tau['LG']); #uA/W, sensitivity
  ##
  R_dyn_SQ =1.
  Lstray = 1.
  Nmux = 128
  In_SQ =1.
  NEP_no_alias = 1.
  ###
  f_nyq=1e3*(R_dyn_SQ/(2.*np.pi*Lstray*2*Nmux)); #kHz
  det_alias_sum=0;
  for ii in range(1,11):  # vector of 1 to 10.
    XX=(NEP_phonon**2 + NEP_photon**2 + 
        johnson**2.*(1 + (4*np.pi*ii*(1e-3*f_nyq)*tau['o'])**2) + 
        shunt**2.*(1 + (4*np.pi*ii*(1e-3*f_nyq)*tau['e'])**2)) / (
        (1+(4*np.pi*ii*(1e-3*f_nyq)*tau['minus'])**2)*(1+(4*np.pi*ii*(1e-3*f_nyq)*tau['plus'])**2))
    det_alias_sum = det_alias_sum + XX;

  NEP_det_alias=np.sqrt(2*det_alias_sum);
  NEP_sq_alias=1e-18*In_SQ*np.sqrt(1+2*Nmux**2*np.sum(1./(np.arange(1,1001))**2.+Nmux**2))/np.abs(S_DC);
  NEP_sq_no_alias=1e-18*In_SQ/np.abs(S_DC);

  NEP_read_tot=np.sqrt(NEP_no_alias**2+NEP_det_alias**2+NEP_sq_alias**2);
  #NEP_MUX_penalty=NEP_tot/NEP_no_alias;

  alias = 1.#NEP_read_tot #placeholder.
  #pdb.set_trace()
  return johnson, shunt, alias

def dPdT_integrand_lambda_sq(nu,T):
  return ( (1 / (np.e ** ((h_ * nu) / (k_b * T)) - 1))**2. * (nu / T)**2. *np.e**((h_ * nu) / (k_b * T)) )

def nep_to_net_Kcmb(nep,cum_eff,band):
  # convert given nep to net_Kcmb.  (from Toki's thesis, added 2 is becuase he ignores 2 pols.)
  # need band.
  # lambda squared still assumed.
  # integrate dP/dT, const and cum_eff in front of int.
  #
  dPdT, error = scint.quad(dPdT_integrand_lambda_sq,band[0],band[1],args=(Tcmb))
  dPdT = dPdT * (2* h_ ** 2 / k_b) * cum_eff  # the coeffs in front of integral

  return nep/(np.sqrt(2)*dPdT)

def nep_to_net_Krj(nep,cum_eff,band):
  # convert nep to K_rj, NET. 
  # assuming lambda squared.
  return nep/(np.sqrt(2)*cum_eff*2*k_b*(band[1]-band[0]))
