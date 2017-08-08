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

# calc illum and spillover efficiencies.  From edge taper.
spill, illum = lf.get_spill_illum_effs(settings.edge_db,settings.f_number)
FWHM = lf.calc_beam(settings.aperture_radius*2., settings.edge_db, settings.c_lambda) # telescope beam in arcmin

for i,elem in enumerate(element_df.element):
  element_out_df.element[i] = elem # put values in export data frame
  element_out_df.temperature[i] = element_df.temperature[i] # put values in export data frame

  # calc element parameters. emissivity.
  if ('lens' in elem) or ('Lens' in elem):
    pass  # add code to get emissivity from loss tan.
  if ('Mirror' in elem) or ('mirror' in elem):
    pass # do scaling of emiss by root(Hz) as per Shaul.  Assuming tophat band, single emiss.
    emissivity = element_df.emissivity[i]*np.sqrt(np.mean(settings.band)/150e9)
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
               emissivity, settings.band)
    # noises are NEP**2  !!!
    poisson_nepsq, poi_err, bunch_nepsq, bunch_err = lf.NEP_photon_per_element(
                 element_df.temperature[i], emissivity, settings.band, cum_eff)
  else:
    power_emit,error = lf.power_emitted_per_element(element_df.temperature[i],
                       emissivity, settings.band)
    poisson_nepsq, poi_err, bunch_nepsq, bunch_err = lf.NEP_photon_per_element(
                 element_df.temperature[i], emissivity, settings.band, cum_eff)

  element_out_df.power_emit[i] = power_emit
  element_out_df.power_absorb[i] = power_emit*cum_eff
  element_out_df.nep_poisson[i] = np.sqrt(poisson_nepsq)
  element_out_df.nep_bunch[i] = np.sqrt(bunch_nepsq)

# save elements data frame
element_out_df.to_csv(settings.elements_out_path, index=True)

p_opt = np.sum(element_out_df.power_absorb)

# calc NEP_bunch for full system.
nep_all_bunch_sq, error = scint.quad(lf.bunch_integrand_lambda_sq_multiple_surfaces, 
                          settings.band[0], settings.band[1],args=(element_df.temperature, 
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
net_total = lf.nep_to_net_Kcmb(nep_total,cum_eff,settings.band)  # K_cmb / rt(Hz)


# write out NEPs
  
print element_out_df
print 'Running :    ', settings.freq, settings.version 
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

pdb.set_trace()




# all above should be wrapped so new input edge taper can rerun it all.
