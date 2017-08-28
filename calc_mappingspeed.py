# -*- coding: utf-8 -*-
'''
Cleaned up mapping speed and loading codes.  Original was from Ben Westbrook.
This is the stripped down version to do only what I need, but no more.
Currently does noise per frequency.  Doesn't map vs pixel size also.  Just uses edge taper.

Files important for this code:
a settings file!!
loading_functions.py
FA_noise_functions.py (franky's noise functions)
pywtl (for franky's readout noise)
a csv of element inputs.
a csv of band inputs (if doing multiple bands).


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
from pywtl.core.wtl_ConvertUtils import convert_squid
import pywtl.common.analysis.noise.analysis.NoisePred as NP
import pywtl.common.analysis.noise.analysis.ParameterLib as PL
import pywtl.common.analysis.noise.analysis.DataHandlingLib as DataHand
import FA_noise_functions as FA

import pdb

import loading_functions as lf
from default_settings import *

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
def power_noise_calculation(element_df, element_out_df,results_df, band, c_freq, location=0):
  # inputs are:
    # data frames (input, elements out, results)
    # band (this may change...)
    # center freq, a string for labeling.  really is band name.
    

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
    elif 'CMB' in elem:
      element_out_df.transmission[i] = 1  # due to how I'm using emissivity = absorbtion
    else:
      element_out_df.transmission[i] = 1-emissivity-element_df.reflect_loss[i]

    # calc cumulated efficiency (i.e. power emitted from this layer to detector)
    if i == 0:
        cum_eff = 1 # since there is no element between the first element and detector
    elif ('stop' in elem) or ('Stop' in elem):
        cum_eff = np.prod(element_out_df.transmission[0:i])*(1-spill) # since 1-spill portion of lyot is seen by detector side elements.
    else:
        cum_eff = np.prod(element_out_df.transmission[0:i])
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
                   element_df.temperature[i], emissivity, band, element_out_df.cum_eff[i])
    else:
      power_emit,error = lf.power_emitted_per_element(element_df.temperature[i],
                         emissivity, band)
      poisson_nepsq, poi_err, bunch_nepsq, bunch_err = lf.NEP_photon_per_element(
                   element_df.temperature[i], emissivity, band, element_out_df.cum_eff[i])

    element_out_df.power_emit[i] = power_emit
    element_out_df.power_absorb[i] = power_emit*element_out_df.cum_eff[i]
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
  # write bolo properties to bolo_char, dfmux_settings, and squid_settings
  # allows me to call Franky's functions directly.
  bolo_char = settings.bolo_char
  dfmux_settings = settings.dfmux_settings
  squid_settings = settings.squid_settings # all 3 lines for less typing later.

  bolo_char['nu'] = np.mean(band)
  bolo_char['dnu'] = band[1]-band[0]
  bolo_char['Tc'] = t_c
  bolo_char['Tbolo'] = t_c
  bolo_char['Popt'] = p_opt
  bolo_char['Psat'] = p_sat
  bolo_char['Gbar'] = G_bar
  bolo_char['G'] = G_dyn
  bolo_char['gamma'] = gamma

  dfmux_settings['Vb'] = v_bias
  dfmux_settings['Ac'] = dfmux_settings['Vb'] / settings.conv.DDStoVbias(Carrier_amplitude=1,
                                                                Carrier_gain=dfmux_settings['Gc'],
                                                                firmware_version='16ch')
  bolo_char['Si'] = NP.Si(dfmux_settings['Vb'])  # assumes deep in transition

  R_gain = [2000., 820., 200., 0.]
  dfmux_settings['An'] = dfmux_settings['Ac']/3. * \
     (R_gain[dfmux_settings['Gn']] + 100.) / (R_gain[dfmux_settings['Gc']] + 100.)

  # NEPs, bolometer  (Karl's method.  inaccurate readout)
  nep_phonon,nep_johnson,nep_readout_KY = lf.calc_bolo_noises(G_dyn,t_c, settings.bolo_resistance, 
                                       v_bias, gamma, settings.readout_noise_amps)

  # NEPs, bolometer in transition (Francois's method, accurate readout)
  if settings.noise_type != 'transition':
    raise ValueError("Only 'transition' is a valid nosie type. Please check settings.")

  noise_pred = NP.addTransitionNoise(dfmux_settings, squid_settings, bolo_char, 
                                     settings.R_wire, settings.cryo)
  # parameter and noise text
  tags = {'title': '', 'squid': '', 'wafer': '', 'dfmux': ''}
  full_text = DataHand.CreateParameterText(tags, 0, np.nan, np.nan, dfmux_settings, 
                                           squid_settings, bolo_char, settings.R_wire, 
                                           noise_pred[0]*1e12, 0., 0., '', noise_pred[1])
  # removes summary part
  full_text = full_text.split('Predicted noise')[0]
  full_text = "\n".join([settings.name, full_text])

  # gets noise by type
  noise = FA.extract_noise_from_text(settings, noise_pred[1]) ## in pA / rt(Hz)

  # custom noise summary
  noise_summary, noise_list = FA.summarize_noise(settings, noise)
  if settings.verbose:
    print full_text+noise_summary

  #pdb.set_trace()
  # get out franky noise that I need.  No trailing subscript = trusted code.  KY or FA means
  # Karl's of Franky's calculations that are less trusted and not used in final product.
  nep_readout =  np.sum(np.array([noise['warm'], noise['cold']])**2.)**0.5 * \
                        1.e-12 / NP.Si(settings.dfmux_settings['Vb'], False, 
                                    settings.bolo_char['L'])
  nep_johnson_FA = noise['johnson'] * 1.e-12 / NP.Si(settings.dfmux_settings['Vb'], False, 
                                                  settings.bolo_char['L'])
  nep_photon_FA = noise['photon'] * 1.e-12 / NP.Si(settings.dfmux_settings['Vb'], False, 
                                                settings.bolo_char['L'])
  nep_phonon_FA = noise['phonon'] * 1.e-12 / NP.Si(settings.dfmux_settings['Vb'], False, 
                                                settings.bolo_char['L'])

  # sum up noise, various combos as double checks.
  nep_total_KY = np.sum(np.array([nep_photon,nep_phonon,nep_johnson,nep_readout_KY])**2.)**0.5
  nep_total_FA = noise_list[1]
  nep_total = np.sum(np.array([nep_photon,nep_phonon,nep_johnson,nep_readout])**2.)**0.5

  # convert NEP to NET
  net_total = lf.nep_to_net_Kcmb(nep_total,element_out_df.cum_eff[len(element_df)-1],band)  # K_cmb / rt(Hz)


  # write out NEPs
  if settings.verbose:
    print element_out_df
    print 'Running :    ', c_freq, settings.version 
    print 'spillover: ', spill
    print 'beam FWHM: ', FWHM
    print 'eff@ bolo: ', element_out_df.cum_eff[len(element_df)-1]
    print 'illumin  : ', illum
    print 'total pow: ', p_opt
    print 'all bunch: ', nep_all_bunch_sq**.5
    print 'all poiss: ', nep_all_poisson_sq**.5
    print 'all phot : ', nep_photon
    print 'all phon : ', nep_phonon
    print 'all john : ', nep_johnson
    print 'FA  read : ', nep_readout
    print 'KY  read : ', nep_readout_KY
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

  bolo_char_out_df.Psat[location] = p_sat
  bolo_char_out_df.Gbar[location] = G_bar
  bolo_char_out_df.Gdyn[location] = G_dyn
  bolo_char_out_df.Tc[location] = t_c
  bolo_char_out_df.Vb[location] = v_bias

  #pdb.set_trace()

  return element_out_df, results_df, bolo_char_out_df


#pdb.set_trace()


if settings.mult_bands is True:
  # read in bands csv
  bands = pandas.read_csv(settings.bands_path)
  # make output df
  results_df = pandas.DataFrame(index=bands.index, columns=['Band','nu','nu_low',
                      'nu_high','spill_eff','illum_eff','FWHM', 'total_pow','NEP_total',
                      'NET_total', 'NEP_poisson', 'NEP_photon', 'NEP_phonon', 
                      'NEP_johnson','NEP_readout'])

  bolo_char_out_df = pandas.DataFrame(index=bands.index, columns=['Band', 'Psat','Gbar','Gdyn', 'Tc', 'Vb'])

  # copy data over
  results_df.Band = bands.Band
  results_df.nu = bands.nu
  results_df.nu_low = bands.nu_low
  results_df.nu_high = bands.nu_high

  bolo_char_out_df.Band = bands.Band


  for i,center_nu in enumerate(bands.nu):
    band = np.array([bands.nu_low[i],bands.nu_high[i]])*1e9
    c_freq = '%g_GHz' %center_nu
    element_out_df, results_df, bolo_char_out_df = power_noise_calculation(element_df, element_out_df,
                                 results_df, band, c_freq,location=i)

else:
  results_df = pandas.DataFrame(index=[0], columns=['Band','nu','nu_lo',
                      'nu_high','spill_eff','illum_eff','FWHM', 'total_pow','NEP_total',
                      'NET_total', 'NEP_poisson', 'NEP_photon', 'NEP_phonon', 
                      'NEP_johnson','NEP_readout'])
  bolo_char_out_df = pandas.DataFrame(index=bands.index, columns=['Psat','Gbar','Gdyn', 'Tc', 'Vb'])

  element_out_df, results_df, bolo_char_out_df = power_noise_calculation(element_df, element_out_df, 
                               results_df, settings.band, settings.freq)


  bolo_char_out_df.Band = settings.freq


results_df.to_csv(os.path.join(settings.base_path, 
                 'outputs/%s/All_results_out.csv' %(settings.version)), 
                 index=True)

bolo_char_out_df.to_csv(os.path.join(settings.base_path, 
                 'outputs/%s/Bolo_char_out.csv' %(settings.version)), 
                 index=True)

pdb.set_trace()
