# -*- coding: utf-8 -*-
'''
Cleaned up mapping speed and loading codes.  Original was from Ben Westbrook.
This is the stripped down version to do only what I need, but no more.

Files important for this code:
loading_functions.py
???


'''

import os, sys
import matplotlib
import pandas
import time
import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt
from pprint import pprint
from copy import copy, deepcopy

import pdb

import loading_functions as lf
from default_settings import *


## outline

# load element parameters 
element_df = pandas.read_csv(settings.elements_path)


# make data frame for output. matching row number to elements_in
#(possible to add columns later, see:  https://stackoverflow.com/questions/12555323/adding-new-column-to-existing-dataframe-in-python-pandas
element_out_df = pandas.DataFrame(index=element_df.index, columns=[
                          'element', 'emissivity', 'temperature', 'power_emit', 
                          'transmission', 'cum_eff', 'power_absorb', 'nep_poisson', 
                          'nep_bunch' ])

# code that starts everything at the bottom.


# wrap main calculation as function
def power_noise_calculation(element_df, element_out_df,results_df, band, c_freq, verbose = True, location=0):
  # inputs are:
    # data frames (input, elements out, results)
    # band (this may change...)
    # center freq, a string for labeling.  really is band name.
    # verbose=T/F

  # calc illum and spillover efficiencies.  From edge taper.
  spill, illum = lf.get_spill_illum_effs(settings.edge_db,settings.f_number)
  c_lambda = 299792458 / np.mean(band)
  FWHM = lf.calc_beam(settings.aperture_radius*2., settings.edge_db, c_lambda) # telescope beam in arcmin

  for i,elem in enumerate(element_df.element):
    element_out_df.element[i] = elem # put values in export data frame
    element_out_df.temperature[i] = element_df.temperature[i] # put values in export data frame

    # calc element parameters. emissivity.
    if ('lens' in elem) or ('Lens' in elem):
      pass  # add code to get emissivity from loss tan.
    elif ('Mirror' in elem) or ('mirror' in elem):
      pass # do scaling of emiss by root(Hz) as per Shaul.  Assuming tophat band, single emiss.
      emissivity = element_df.emissivity[i]*np.sqrt(np.mean(band)/150e9)
    else:
      emissivity = element_df.emissivity[i]
    element_out_df.emissivity[i] = emissivity

    # calc efficiency (i.e. transmission through that layer) for element. save to new throughput frame
    if ('stop' in elem) or ('Stop' in elem):
      element_out_df.transmission[i] = spill
      cum_eff = np.prod(element_out_df.transmission[0:i])*(1-spill) # since 1-spill portion of lyot is seen by detector side elements.
    elif 'CMB' in elem:
      element_out_df.transmission[i] = 1  # due to how I'm using emissivity = absorbtion
      cum_eff = np.prod(element_out_df.transmission[0:(i+1)])
    else:
      element_out_df.transmission[i] = 1-emissivity-element_df.reflect_loss[i]
      cum_eff = np.prod(element_out_df.transmission[0:(i+1)])
    element_out_df.cum_eff[i] = cum_eff

  # calc emitted power per element.  save to new throughput frame
  # calc NEP in detector per element.  (uses cum_eff)
  # assuming throughput = lambda ** 2 (in loading_functions.py)
    if ('Mirror' in elem) or ('mirror' in elem):
      # SHOULD change to emiss dependence on freq for alum. 
      power_emit,pow_error = lf.power_emitted_per_element(element_df.temperature[i],
                 emissivity, band)
      # noises are NEP**2  !!!
      poisson_nepsq, poi_err, bunch_nepsq, bunch_err = lf.NEP_photon_per_element(
                   element_df.temperature[i], emissivity, band, cum_eff)
    else:
      power_emit,error = lf.power_emitted_per_element(element_df.temperature[i],
                         emissivity, band)
      poisson_nepsq, poi_err, bunch_nepsq, bunch_err = lf.NEP_photon_per_element(
                   element_df.temperature[i], emissivity, band, cum_eff)

    element_out_df.power_emit[i] = power_emit
    element_out_df.power_absorb[i] = power_emit*cum_eff
    element_out_df.nep_poisson[i] = np.sqrt(poisson_nepsq)
    element_out_df.nep_bunch[i] = np.sqrt(bunch_nepsq)

  # save elements data frame
  save_path = os.path.join(settings.base_path,'outputs/%s' %settings.version)
  
  
  if not os.path.exists(save_path):
    os.mkdir(save_path)

  element_out_df.to_csv(os.path.join(save_path, '%s_elements_out.csv' %(c_freq)), index=True)
    
  p_opt = np.sum(element_out_df.power_absorb)

  # calc NEP_bunch for full system.
  nep_all_bunch_sq, error = scint.quad(lf.bunch_integrand_lambda_sq_multiple_surfaces, 
                            band[0], band[1],args=(element_df.temperature, 
                            element_out_df.emissivity, element_out_df.cum_eff))

  # NEP poisson for all
  nep_all_poisson_sq = np.sum(element_out_df.nep_poisson**2.)
  nep_photon = np.sqrt(nep_all_bunch_sq + nep_all_poisson_sq)

  # bolometer properties
  p_sat, G_bar, G_dyn, t_c, v_bias, gamma = lf.bolo_properties(p_opt,settings.t_bath,
                                     settings.safety_factor,settings.n,settings.bolo_Rn,
                                     settings.bias_point)

  # NEPs, bolometer
  nep_phonon,nep_johnson,nep_readout = lf.calc_bolo_noises(G_dyn,t_c, settings.bolo_resistance, 
                                       v_bias, gamma, settings.readout_noise_amps)

  nep_total = np.sum(np.array([nep_photon,nep_phonon,nep_johnson,nep_readout])**2.)**0.5

  # convert NEP to NET
  net_total = lf.nep_to_net_Kcmb(nep_total,cum_eff,band)  # K_cmb / rt(Hz)


  # write out NEPs
  if verbose:
    print element_out_df
    print 'Running :    ', c_freq, settings.version 
    print 'spillover: ', spill
    print 'beam FWHM: ', FWHM
    print 'eff@ bolo: ', cum_eff
    print 'illumin  : ', illum
    print 'total pow: ', p_opt
    print 'all bunch: ', nep_all_bunch_sq**.5
    print 'all poiss: ', nep_all_poisson_sq**.5
    print 'all phot : ', nep_photon
    print 'all phon : ', nep_phonon
    print 'all john : ', nep_johnson
    print 'all read : ', nep_readout
    print 'total NEP: ', nep_total
    print 'total NET: ', net_total


  results_df.spill_eff[location]   = spill
  results_df.illum_eff[location]   = illum
  results_df.FWHM[location]        = FWHM
  results_df.total_pow[location]   = p_opt
  results_df.NEP_total[location]   = nep_total
  results_df.NET_total[location]   = net_total
  results_df.NEP_poisson[location] = nep_all_poisson_sq**.5
  results_df.NEP_photon[location]  = nep_photon 
  results_df.NEP_phonon[location]  = nep_phonon
  results_df.NEP_johnson[location] = nep_johnson
  results_df.NEP_readout[location] = nep_readout 

  return element_out_df, results_df


#pdb.set_trace()


if settings.mult_bands is True:
  # read in bands csv
  bands = pandas.read_csv(settings.bands_path)
  # make output df
  results_df = pandas.DataFrame(index=bands.index, columns=['Band','nu','nu_low',
                      'nu_high','spill_eff','illum_eff','FWHM', 'total_pow','NEP_total',
                      'NET_total', 'NEP_poisson', 'NEP_photon', 'NEP_phonon', 
                      'NEP_johnson','NEP_readout'])
  # copy data over
  results_df.Band = bands.Band
  results_df.nu = bands.nu
  results_df.nu_low = bands.nu_low
  results_df.nu_high = bands.nu_high



  for i,center_nu in enumerate(bands.nu):
    band = np.array([bands.nu_low[i],bands.nu_high[i]])*1e9
    c_freq = '%g_GHz' %center_nu
    element_out_df, results_df = power_noise_calculation(element_df, element_out_df,
                                 results_df, band, c_freq,location=i)

else:
  results_df = pandas.DataFrame(index=[0], columns=['Band','nu','nu_lo',
                      'nu_high','spill_eff','illum_eff', 'total_pow','NEP_total',
                      'NET_total', 'NEP_poisosn', 'NEP_photon', 'NEP_phonon', 
                      'NEP_johnson','NEP_readout'])
  element_out_df, results_df = power_noise_calculation(element_df, element_out_df, 
                               results_df, settings.band, settings.freq)


results_df.to_csv(os.path.join(settings.base_path, 
                 'outputs/%s/All_results_out.csv ' %(settings.version)), 
                 index=True)


pdb.set_trace()
