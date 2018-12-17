# -*- coding: utf-8 -*-
'''
Cleaned up mapping speed and loading codes.  Original was from Ben Westbrook.
This is the stripped down version to do only what I need, but no more.
Currently does noise per frequency. 
Does include MCPs now.
Can calculate pixel size from edge taper or vice-versa
can scan in edge taper to find best pixel size.


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
import datetime
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
import cosmolab.probe.noise_prediction.noise_prediction as NP_FA
import cosmolab.probe.noise_library.probe_noise_pred as NP_lib
from joblib import Parallel, delayed

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
def power_noise_calculation(element_df, element_out_df, bolo_char_out_df, results_df, bands, band, c_freq, save_path, location=0):
  # inputs are:
    # data frames (input, elements out, bolos out, results)
    # band, 2 element array. in Hz
    # center freq, a string for labeling.  really is band name.
  
  c_lambda = 299792458 / np.mean(band)
    
  if settings.MCP:   # use pixel sizes and lambda to get spill and edge_taper. Need to determine how bands relate amongst pixels.
    if settings.use_edge_dB:
      # calc spill from dB per pixel
      spill, illum = lf.get_spill_illum_effs(results_df.loc[location,('edge_dB')],settings.f_number)  
      theta_px = lf.theta_px_from_edge_db(results_df.loc[location,('edge_dB')],settings.f_number)
      # results_df.loc[location,('D_px')] = bands.loc[location,('D_px')] # lf.get_D_px(c_lambda,theta_px,settings.diameter_to_waist_ratio) decoupling diameter and dB for antenna array case.
	  results_df.loc[location,('D_px')] = lf.get_D_px(c_lambda,theta_px,settings.diameter_to_waist_ratio)
	  

    elif settings.use_D_px:
      # calc spill from D_px.
      theta_px = lf.get_theta_px(c_lambda, results_df.loc[location,('D_px')],settings.diameter_to_waist_ratio)
      spill, illum = lf.get_spill_illum_effs_from_theta_px(theta_px,settings.f_number)
      results_df.loc[location,('edge_dB')] =  lf.edge_db_from_theta_px(theta_px,settings.f_number)

    else:
      print '\n\n No clear idea of what pixel size or edge taper to use!\n\n Check settings file\n\n'
      sys.exit()

    FWHM = lf.calc_beam(settings.aperture_radius*2., results_df.loc[location,('edge_dB')], c_lambda) # telescope beam in arcmin

  else:      # use edge taper to get spillover.
    # calc illum and spillover efficiencies.  From edge taper.
    spill, illum = lf.get_spill_illum_effs(settings.edge_db,settings.f_number)
    FWHM = lf.calc_beam(settings.aperture_radius*2., settings.edge_db, c_lambda) # telescope beam in arcmin
    # also calc D_px just because I can.  And I will want it later.
    theta_px = lf.theta_px_from_edge_db(settings.edge_db,settings.f_number)
    results_df.loc[location,('D_px')] = lf.get_D_px(c_lambda,theta_px,settings.diameter_to_waist_ratio)
    results_df.edge_db = settings.edge_db # reset this since I'm using edge_db everywhere.

  for i,elem in enumerate(element_df.element):
    element_out_df.element[i] = elem # put values in export data frame
    element_out_df.temperature[i] = element_df.temperature[i] # put values in export data frame

    # calc element parameters. emissivity.
    if ('lens' in elem) or ('Lens' in elem):
      if element_df.loss_tan[i] == 'PP':  # pull out polypropylene case so can have loss tan vary with frequency.
        # hardcoded loss_tan!!  :-(  # 21 values for the 21 bands.  
        # frequencies [20.8,  25,  30,  36,  43,     52,     62,     75,   90,  108,    129,    155,    186,  223,    268,  321,   385,     462,   555,   666,  799]
        loss_tan =    np.array([5e-4,5e-4,5e-4,5e-4,5e-4, 5.8e-4, 6.4e-4, 7.2e-4, 8e-4, 7e-4, 6.6e-4, 6.4e-4, 6.2e-4, 6e-4, 5.8e-4, 6e-4, 9e-4,  13.5e-4, 20e-4, 31e-4, 50e-4]) / 5.  # reduce by 5 for 300 K to 4 K difference
        # hardcoded thickness scaling as well!  Scale thickness by center wavelength of band.  Scale upper 3 bands by 1 number.  scaling is 1 at 150 GHz.  scaling is linear with wavelength.
        thick_sc = [5.,  4.2,  5.,  4.2,  5., 4.2,  1.7,  1.4,   1.7,  1.4,  1.7,  1.4,  .56,  .47, .56,  .47, .56, .47, 0.225, 0.225,  0.225]
        emissivity = 1-np.exp(-2*np.pi*element_df.opt_index[i] * loss_tan[location]
                             * element_df.thickness[i] * thick_sc[location] / c_lambda)   # location is band we are on.
        # pdb.set_trace()
      else:
        # code to get emissivity from loss tan index and thickness.
        emissivity = 1-np.exp(-2*np.pi*element_df.opt_index[i] * element_df.loss_tan[i]
                             * element_df.thickness[i] / c_lambda)
    elif ('Mirror' in elem) or ('mirror' in elem):
      if settings.emiss_rt_GHz:
        # do scaling of emiss by root(Hz) as per Shaul.  Assuming tophat band, single emiss.
        emissivity = element_df.emissivity[i]*np.sqrt(np.mean(band)/150e9)
      else:
        # fixed emissivity at all frequencies ala Planck HFI instrument paper. Ade 2011,https://doi.org/10.1051/0004-6361/201116486 
        emissivity = element_df.emissivity[i]		
    elif ('ATM' in elem) or ('atm' in elem):  
      # load data from atm code, orginally from Ben W.
	  import atmospheric_calculator as ac
	  
	  ######################################   HARD CODED   ########################  atm file.
	  
	  atmosphere_data = ac.atmospheric_calculator('./inputs/LDB_34km_30deg_el.out',
                                                    band[0], band[1])  # freqs in Hz.
	  # print 'Atmosphere properteris are emissivity = %.4f, transmission = %.4f' %(atmosphere_data[0],atmosphere_data[1])
	  # print 'only using emissivity, assuming R = 0. T = 1-E'
	  # f = open('temp.txt', 'a')
	  # f.write('%.2f - %.2f GHz, emiss = %.5f\n' %(band[0]*1e-9, band[1]*1e-9, atmosphere_data[0]))
	  # f.close()
	  
	  emissivity = atmosphere_data[0]
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
    power_emit,error = lf.power_emitted_per_element(element_df.temperature[i],
                         emissivity, band)
    poisson_nepsq, poi_err, bunch_nepsq, bunch_err = lf.NEP_photon_per_element(
                   element_df.temperature[i], emissivity, band, element_out_df.cum_eff[i])

    element_out_df.power_emit[i] = power_emit
    # if 'ATM' in elem: 
	  # print 'ATM power: (emit), (absorb)', power_emit, power_emit*element_out_df.cum_eff[i], band
	  
    element_out_df.power_absorb[i] = power_emit*element_out_df.cum_eff[i]
    element_out_df.nep_poisson[i] = np.sqrt(poisson_nepsq)
    element_out_df.nep_bunch[i] = np.sqrt(bunch_nepsq)
    
  ## last element calced is CMB.  save that term seperately
  NEP_bunch_CMB = element_out_df.nep_bunch[i]

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
  '''
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
  '''
  #################################3 Franky's new optimized noise below.
   
  mux_app = NP_FA.NoisePred() ## load noise class
  mux_app.settings = settings
  ## noise library loaded before all the parameters are calculated
  mux_app.NP = NP_lib.NoisePredLib(settings.system)
  
  ## calculating derived parameters
  
  if settings.verbose == False: sys.stdout = open(os.devnull, "w")
  mux_app.calculate_parameters()
  if settings.verbose == False: sys.stdout = sys.__stdout__
  #HOW TO load FA coad without loading wrong settings file? it looks in current directory and loads right thing. yay.

  ## noise library loaded after all the parameters are calculated
  mux_app.NP = NP_lib.NoisePredLib(settings.system)

  ## creates output
  #self.create_output()

  # # loading parameters
  # dfmux_settings = self.settings.dfmux_settings
  # squid_settings = self.settings.squid_settings
  # bolo_char = self.settings.bolo_char

  valid_noise_type = ["transition"]  # , "darktransition", "overbias", "resistor", "darksquid"]
  if settings.noise_type not in valid_noise_type:
    raise ValueError("Noise type '{0:s}' not valid, enter one of these: {1:s}".format(settings.noise_type,
                                                                                        valid_noise_type))
  # noise calculation
  if settings.noise_type == "transition":
    noise_pred = mux_app.NP.addTransitionNoise(dfmux_settings, squid_settings, bolo_char)

  # parameter and noise text
  tags = {'title': '', 'squid': '', 'wafer': '', 'dfmux': ''}
  full_text = mux_app.CreateParameterText(dfmux_settings, squid_settings, bolo_char, noise_pred[0]*1e12,
                                       noise_pred[1])

  # removes summary part
  full_text = full_text.split('Predicted noise')[0]
  full_text = "\n".join([settings.name, full_text])

  # gets noise by type
  noise = mux_app.extract_noise_from_text(noise_pred[1])

  # custom noise summary
  noise_summary, noise_list = mux_app.summarize_noise(noise)
  if settings.verbose:
    print full_text+noise_summary

  ###########################################
  # TDM readout noise calculated ala Roger O'brient. readout and johnson are both different.  
  # will keep or toss various terms later as needed.

  #nep_johnson_TDM, nep_shunt_TDM, NEP_alias_TDM = lf.TDM_readout_noise(nep_photon, nep_phonon, 
  #                                                p_opt, t_c, settings.t_bath, G_dyn, 
  #                                                settings.alpha, settings.n, settings.R_op_TDM, 
  #                                                settings.Rs_op_TDM, settings.safety_factor, settings.C0,
  #                                                settings.L_sq_nyq,)
  
  #pdb.set_trace()
      
  ##############################################


  # get out franky noise that I need.  No trailing subscript = trusted code.  KY or FA means
  # Karl's or Franky's calculations that are less trusted and not used in final product.
  #
  # or can use TDM readout and johnson as a good alternative. fun times.


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
  net_total = lf.nep_to_net_Kcmb(nep_total,element_out_df.cum_eff[len(element_df)-1],band)  # K_cmb * rt(sec)
  net_rj_total = lf.nep_to_net_Krj(nep_total,element_out_df.cum_eff[len(element_df)-1],band)  # K_rj * rt(sec)

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


  results_df.loc[location,('spill_eff')] = spill
  results_df.loc[location,('illum_eff')]   = illum
  results_df.loc[location,('FWHM')]        = FWHM
  results_df.loc[location,('total_pow')]   = p_opt
  results_df.loc[location,('NEP_total')]   = nep_total
  results_df.loc[location,('NET_total')]   = net_total
  results_df.loc[location,('NET_rj_total')]   = net_rj_total
  results_df.loc[location,('NEP_poisson')] = nep_all_poisson_sq**.5
  results_df.loc[location,('NEP_photon')]  = nep_photon
  results_df.loc[location,('NET_bunch_all')]  = lf.nep_to_net_Kcmb(nep_all_bunch_sq**.5,element_out_df.cum_eff[len(element_df)-1],band) 
  results_df.loc[location,('NET_bunch_CMB')]  = lf.nep_to_net_Kcmb(NEP_bunch_CMB,element_out_df.cum_eff[len(element_df)-1],band)
  results_df.loc[location,('NEP_phonon')]  = nep_phonon
  results_df.loc[location,('NEP_johnson')] = nep_johnson
  results_df.loc[location,('NEP_readout')] = nep_readout 

  bolo_char_out_df.loc[location,('Psat')] = p_sat
  bolo_char_out_df.loc[location,('Gbar')] = G_bar
  bolo_char_out_df.loc[location,('Gdyn')] = G_dyn
  bolo_char_out_df.loc[location,('Tc')] = t_c
  bolo_char_out_df.loc[location,('Vb')] = v_bias
  bolo_char_out_df.loc[location,('gamma')] = gamma

  #pdb.set_trace()

  return element_out_df, results_df, bolo_char_out_df


#pdb.set_trace()

def make_data_frames(index=[0]):
  results_cols=['Band','nu','nu_low','nu_high','D_px','edge_dB','spill_eff','illum_eff','FWHM', 'total_pow','NEP_total',
                'NET_total','NET_rj_total', 'NEP_poisson', 'NEP_photon', 'NEP_phonon','NEP_johnson','NEP_readout',
                'NET_bunch_all','NET_bunch_CMB']

  results_df = pandas.DataFrame(index=index, columns=results_cols)
  bolo_char_out_df = pandas.DataFrame(index=index, columns=['Band', 'Psat','Gbar','Gdyn', 'Tc', 'Vb','gamma'])

  return results_df, bolo_char_out_df

def calc_MCP_sizes(results_df, bolo_char_out_df, bands):
  # copy data over
  results_df.Band = bands.Band
  results_df.nu = bands.nu
  results_df.nu_low = bands.nu_low
  results_df.nu_high = bands.nu_high

  bolo_char_out_df.Band = bands.Band
  if settings.MCP:
    # calculate (or load) pixel properties to use
    if settings.use_edge_dB:      # calc edge dB for each band given fixed pixels.

      pixel_types = []
      for i in bands.loc[:,('pixel')]:
        if i not in pixel_types:
          pixel_types.append(i)

      for pixel in pixel_types:
        ## get the bands per pixel
        px_bands = bands.loc[bands.loc[:,('pixel')]==pixel,('Band')]  # bands in this pixel.

        ## based on length do soemthign --- # 3 cases, tri-chroic, bi-chroic, single color.
        if len(px_bands) == 3:

          '''                    
          # set low band to edge_dB
          results_df.loc[(px_bands.index[0]),('edge_dB')] = settings.edge_db
          # calc mid and upper bands
          #lower
          results_df.loc[(px_bands.index[1]),('edge_dB')] = lf.scale_db_between_bands(settings.edge_db, 
                                                              bands.nu[px_bands.index[0]], bands.nu[px_bands.index[1]])
          #upper
          results_df.loc[(px_bands.index[2]),('edge_dB')] = lf.scale_db_between_bands(settings.edge_db, 
                                                              bands.nu[px_bands.index[0]], bands.nu[px_bands.index[2]])
          '''
          # set mid band to edge_dB
          results_df.loc[(px_bands.index[1]),('edge_dB')] = settings.edge_db
          # calc low and upper bands
          #lower
          results_df.loc[(px_bands.index[0]),('edge_dB')] = lf.scale_db_between_bands(settings.edge_db, 
                                                              bands.nu[px_bands.index[1]], bands.nu[px_bands.index[0]])
          #upper
          results_df.loc[(px_bands.index[2]),('edge_dB')] = lf.scale_db_between_bands(settings.edge_db, 
                                                              bands.nu[px_bands.index[1]], bands.nu[px_bands.index[2]])


        elif len(px_bands) == 2:
          # set edge_dB between bands at geometric mean.
          mid_freq = np.sqrt(bands.nu[px_bands.index[0]] * bands.nu[px_bands.index[1]])

          # lower band
          results_df.loc[(px_bands.index[0]),('edge_dB')] = lf.scale_db_between_bands(settings.edge_db, 
                                                              mid_freq, bands.nu[px_bands.index[0]])
          # upper band
          results_df.loc[(px_bands.index[1]),('edge_dB')] = lf.scale_db_between_bands(settings.edge_db, 
                                                              mid_freq, bands.nu[px_bands.index[1]])
          '''
          # set lower band to edge_dB
          results_df.loc[(px_bands.index[0]),('edge_dB')] = settings.edge_db
          # calc edge_dB of upper band
          results_df.loc[(px_bands.index[1]),('edge_dB')] = lf.scale_db_between_bands(settings.edge_db, 
                                                              bands.nu[px_bands.index[0]], bands.nu[px_bands.index[1]])
          '''

        elif len(px_bands) == 1:
          results_df.loc[(px_bands.index[0]),('edge_dB')] = settings.edge_db  ## indexed row by band # minus 1. diff or 0--20 or 1--21 in df.index and df.bands.
        else:
          print 'Some error with number of bands in pixel %s. \nExiting...' %pixel
          sys.exit()

    if settings.use_D_px:  # load from bands.csv dataframe
      results_df.loc[:,('D_px')] = bands.loc[:,('D_px')]
  else:
    results_df.loc[:,('edge_dB')] = settings.edge_db
    
  return results_df, bolo_char_out_df, bands

def run_bands_case(element_df, element_out_df, save_path):
  # (element_df, element_out_df,results_df, band, c_freq,  ## from 1 freq noise function. 
  if settings.mult_bands is True:
    bands = pandas.read_csv(settings.bands_path)
    results_df, bolo_char_out_df = make_data_frames(index=bands.index)
    results_df, bolo_char_out_df, bands = calc_MCP_sizes(results_df, bolo_char_out_df, bands)

    # run the main code
    for i,center_nu in enumerate(bands.nu):
      band = np.array([bands.nu_low[i],bands.nu_high[i]])*1e9
      c_freq = '%g_GHz' %center_nu
      element_out_df, results_df, bolo_char_out_df = power_noise_calculation(element_df, element_out_df,
                                   bolo_char_out_df, results_df, bands, band, c_freq, save_path, location=i)

  else:
    results_df, bolo_char_out_df = make_data_frames()
    element_out_df, results_df, bolo_char_out_df = power_noise_calculation(element_df, element_out_df, 
                                 bolo_char_out_df, results_df, bands, settings.band, settings.freq, save_path)

    bolo_char_out_df.Band = settings.freq

  return results_df, bolo_char_out_df
'''
def add_NET_and_correlated_noise(results):
  # adds NET_array and pol_weight for correlated noise values
  # assumes a fixed area of diameter = 50*F*lambda.  Only run when optimizing pixel size.
  # assumes 2 polarizations are uncorrelated.
  results.loc[:,'num_bolos'] = pandas.Series(np.nan, index=results.index)
  results.loc[:,'NET_array'] = pandas.Series(np.nan, index=results.index)
  results.loc[:,'NET_rj_array'] = pandas.Series(np.nan, index=results.index)
  results.loc[:,'pol_weight'] = pandas.Series(np.nan, index=results.index)
  results.loc[:,'corr_NET_array'] = pandas.Series(np.nan, index=results.index)
  results.loc[:,'corrCMB_NET_array'] = pandas.Series(np.nan, index=results.index)
  results.loc[:,'corr_pol_weight'] = pandas.Series(np.nan, index=results.index)

  # get px size -- calc N_airy
  c_lambda = 0.299792458 / results.nu  # wavelengths
  D_airy = 1.22*c_lambda*settings.f_number * 2  # equation is for radius of airy, x2 for diameter.
  N_airy = (D_airy / results.D_px)**2. * 0.9069 # both are circles. 0.9 is for hex packing.
  N_airy = np.array([1. if i < 1 else i for i in N_airy]) # allowing fractional pixels, but not values < 1

  # num pixels
  #area = (50*settings.f_number*c_lambda)**2.  # do 50 F lambda diameter
  area = (0.20)**2.  # do 20cm diameter
  N_px = area*0.9069 / results.D_px**2. # really just ratio of diameters squared. pi's cancell, hex packing is factor of .9069.
  results.num_bolos = N_px

  results.NET_array = results.NET_total / N_px**.5
  results.NET_rj_array = results.NET_rj_total / N_px**.5
  results.corr_NET_array = (results.NET_array**2. + results.NET_bunch_all**2./(N_px*2.)*(N_airy-1))**.5
  results.corrCMB_NET_array = (results.NET_array**2. + results.NET_bunch_CMB**2./(N_px*2.)*(N_airy-1))**.5

  # NET_new**2 = NET_ar**2 + bose**2/N(N_airy-1)  
  yrs2sec_deg2arcmin = np.sqrt(settings.sky_area * 3600. / ( 31557600.))
  results.corr_pol_weight = results.corr_NET_array * yrs2sec_deg2arcmin  / np.sqrt(settings.mission_length) * np.sqrt(2.)
  results.pol_weight = results.NET_array * yrs2sec_deg2arcmin  / np.sqrt(settings.mission_length) * np.sqrt(2.)

  return results
'''

# save path for data frames
main_path = os.path.join(settings.base_path,'outputs/%s' %settings.version)

if not os.path.exists(main_path):
  os.mkdir(main_path)

now = datetime.datetime.now().strftime('%Y%m%d_%H%M')
save_path = os.path.join(main_path, '%s' %(now))
if not os.path.exists(save_path):
  os.mkdir(save_path)

if settings.dB_scan:
  # make place to save pages.
  results_dict = {}
  bolo_dict = {}

  ## add joblib here.
  def loop_func(taper):
    settings.edge_db = taper
    # run all.  
    #print 'Running edge taper: ', taper
    results_df, bolo_char_out_df = run_bands_case(element_df, element_out_df, save_path) 
    results_df = add_NET_and_correlated_noise(results_df)
    return results_df, bolo_char_out_df

  results = Parallel(verbose=10, n_jobs=5)(delayed(loop_func)(taper) for taper in settings.dB_array)

  for i,taper in enumerate(settings.dB_array):  
    results_dict[taper] = results[i][0]
    bolo_dict[taper] = results[i][1]

  '''
  for taper in settings.dB_array:  
    settings.edge_db = taper
    # run all.  
    print 'Running edge taper: ', taper
    results_df, bolo_char_out_df = run_bands_case(element_df, element_out_df, save_path) 
    results_df = add_NET_and_correlated_noise(results_df)

    # save as 1 frame of dict.
    results_dict[taper] = results_df
    bolo_dict[taper] = bolo_char_out_df
  '''
  
  ## save final stuff
  writer1 = pandas.ExcelWriter(os.path.join(save_path, 'All_results_out.xlsx'))
  writer2 = pandas.ExcelWriter(os.path.join(save_path,'Bolo_char_out.xlsx'))

  for taper in settings.dB_array:  
    results_dict[taper].to_excel(writer1,'%.2f' %taper)
    bolo_dict[taper].to_excel(writer2,'%.2f' %taper)

  writer1.save()
  writer2.save()

else:
  # run all.  
  results_df, bolo_char_out_df = run_bands_case(element_df, element_out_df, save_path) 

  ## save stuff
  results_df.to_csv(os.path.join(save_path, 'All_results_out.csv'), index=True)
  bolo_char_out_df.to_csv(os.path.join(save_path,'Bolo_char_out.csv'), index=True)


# link to most recent data.
os.system('ln -sfn %s outputs/%s/current_directory' %(save_path,settings.version))

pdb.set_trace()
