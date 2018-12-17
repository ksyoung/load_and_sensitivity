# -*- coding: utf-8 -*-
'''
Code to read in multiple .csv data frames from calc_mapping_speed.py and plot some comparison plots.


Karl Young, August 28 2017

'''

import os, sys
import matplotlib
import pandas
#import time
import numpy as np
#import scipy.integrate as scint
import matplotlib.pyplot as plt
import optparse
#from pprint import pprint
#from copy import copy, deepcopy
from tabulate import tabulate
import pdb

from default_settings import *

matplotlib.rcParams.update({'font.size': 18}) # set good font sizes

# optparse it!
usage = "usage: %prog  -f <folder_to_read_data_from>"
parser = optparse.OptionParser(usage)
parser.add_option('-f', dest='folder', action='store', type='str', default=None, 
                  help='folder which all results_df is loaded from.  Finds All_results_out.csv in this folder.')
(option, args) = parser.parse_args()

# load bands (has px numbers)
bands = pandas.read_csv(settings.bands_path)
# path to newest data.  May need to change this.......

if option.folder is None:
  print 'Loading newest data for %s\n' %settings.version
  load_path = os.path.join(settings.base_path,'outputs/%s/current_directory' %settings.version)
else:
  print 'Loading data for %s from %s\n' %(settings.version, option.folder)
  load_path = os.path.join(settings.base_path,'%s' %(option.folder))
  

# load bolo char (shouldn't need it)
bolo_char = pandas.read_csv(os.path.join(load_path,'Bolo_char_out.csv'),index_col=0)
# load all results.
results = pandas.read_csv(os.path.join(load_path, 'All_results_out.csv'),index_col=0)

## add new columns to results
results.loc[:,'num_bolo'] = pandas.Series(np.nan, index=results.index)
results.loc[:,'NET_array'] = pandas.Series(np.nan, index=results.index)
results.loc[:,'weight'] = pandas.Series(np.nan, index=results.index)
results.loc[:,'pol_weight'] = pandas.Series(np.nan, index=results.index)


## def up some sweet sweet functions
def check_areas(px_count,px_areas,area_diffs):
  # inputs: number px per type, area used per type, availible strehl per type.
  used_areas = px_count*px_areas # recalc areas.
  if any(used_areas > area_diffs):
    print '\t\tPossible warning!!\n Too many pixel in Band(s) %s\n' %', '.join(pixel_types[used_areas > area_diffs])
  
  # check if all inside a given strehl contour are ok.
  sum_used_A = np.flipud(np.cumsum(np.flipud(used_areas)))
  sum_fp_areas = np.flipud(np.cumsum(np.flipud(area_diffs)))

  if any(sum_used_A > sum_fp_areas):
    print '\t\tWARNING!!\n Too many pixels inside countour(s) %s\n' %', '.join(pixel_types[sum_used_A > sum_fp_areas])
    print 'Reduce area(s) by, ', (sum_used_A-sum_fp_areas)[sum_used_A > sum_fp_areas]

  if np.sum(used_areas) < np.sum(area_diffs):
    print 'Total area is ok.\n\tAvailible:  %.3f\n\t     Used:  %.3f' %(np.sum(area_diffs),np.sum(used_areas))
  else:
    print 'Total area FAILS.\n\tAvailible:  %.3f\n\t     Used:  %.3f\n\n' %(np.sum(area_diffs),np.sum(used_areas))

  return

def print_status(px_types, px_count, single_px_area, used_areas, band_areas):
  # display current state of pixels.
  print 'Current pixel counts and areas:\n'
  table = [row for row in zip(px_types,px_count,single_px_area,used_areas,band_areas)]
  print tabulate(table, headers=['pixel type', 'count','area of pixel', 'area used', 'area availible'],tablefmt='orgtbl') 

  return


def add_correlated_noise(results):
  # adds NET_array and pol_weight for correlated noise values
  # needs to be called after pixel numbers are calculated
  # assumes 2 polarizations are uncorrelated.
  results.loc[:,'corr_NET_array'] = pandas.Series(np.nan, index=results.index)
  results.loc[:,'corrCMB_NET_array'] = pandas.Series(np.nan, index=results.index)
  results.loc[:,'corr_pol_weight'] = pandas.Series(np.nan, index=results.index)

  # get px size -- calc N_airy
  c_lambda = 0.299792458 / results.nu  # wavelengths
  D_airy = 1.22*c_lambda*settings.f_number * 2  # equation is for radius of airy, x2 for diameter.
  N_airy = (D_airy / results.D_px)**2. * 0.9069 # 0.9 for hex packing.
  N_airy = np.array([1. if i < 1 else i for i in N_airy])

  results.corr_NET_array = np.sqrt(results.NET_array**2. + results.NET_bunch_all**2./results.num_bolo*(N_airy-1))
  results.corrCMB_NET_array = np.sqrt(results.NET_array**2. + results.NET_bunch_CMB**2./results.num_bolo*(N_airy-1))

  # NET_new**2 = NET_ar**2 + bose**2/N(N_airy-1)  
  yrs2sec_deg2arcmin = np.sqrt(settings.sky_area * 3600. / ( 31557600.))  # sky area is in sq. deg.
  results.corr_pol_weight = results.corr_NET_array * yrs2sec_deg2arcmin  / np.sqrt(settings.mission_length) * np.sqrt(2.)

  return results

def calc_array_noise(results):
  # calculate per band sensitivity, mapping speed, and pol weight
  results.NET_array = results.NET_total / np.sqrt(results.num_bolo)
  results.map_speed = ((results.FWHM / 60.) / results.NET_array)**2.  ## FWHM**2 / NET **2 , deg^2 / (K^2 sec)
  yrs2sec_deg2arcmin = np.sqrt(settings.sky_area * 3600. / ( 31557600.))
  results.weight = results.NET_array * yrs2sec_deg2arcmin  / np.sqrt(settings.mission_length)  ## array_NET * sqrt(sky_area / time) K * arcmin
  results.pol_weight = results.weight * np.sqrt(2.) ## extra root 2 for polarization data.

  if settings.calc_correlated_noise:
    results = add_correlated_noise(results)
  return results


def total_CMB_sense(results):
  # takes px number and NET per pixel.
  # gives final CMB sensitivity of current system.
  
  # sum 1/squares of noise, take sqrt, invert.
  total_pol = 1./np.sqrt(np.sum(np.power(results.pol_weight,-2.)))
  total_corrPol = 1./np.sqrt(np.sum(np.power(results.corr_pol_weight,-2.)))
  print 'total uK arcmin uncorrelate:  %.2f' %(total_pol*1.e6)
  print 'total uK arcmin correlated :  %.2f' %(total_corrPol*1.e6)
  return total_pol, total_corrPol

pixel_types = []
for i in bands.loc[:,('pixel')]:
  if i not in pixel_types:
    pixel_types.append(i)
pixel_types = np.array(pixel_types) ## later code works if it's an np array.

px_count = np.zeros(len(pixel_types))
px_areas = np.zeros(len(pixel_types))
for i,pixel in enumerate(pixel_types):
  px_indices = np.where(bands.pixel==pixel)[0]

  ## just use first entry for D_px, all should be equal. then get area
  px_areas[i] = np.pi*results.D_px[px_indices[0]]**2. / 4. 


if settings.calc_N_px_by_area_csv:
  print 'Estimating number of pixel which fit in focal plane.\n'

  FP_areas = pandas.read_csv(os.path.join(settings.FP_areas_path))
  area_diffs = -np.diff(FP_areas.area) * 0.9069   # area for each frequency. include hex pack factor

  ######## 3 single band pixels hardcoded!!!  ###############

  # how to deal with single band pixels?!?!??!?!?!?!?
  # not using most of the area anyway......
  # for G,H,I.  do even thirds at smallest area.
  # area_diffs = np.append(area_diffs,[FP_areas.area[FP_areas.index[-1]]*0.9069/3]*3)  # add in area for last 3 rows.
  # area_diffs = np.append(area_diffs,[FP_areas.area[FP_areas.index[-1]]*0.9069/2]*2)  # add in area for last 2 rows.
  area_diffs = np.append(area_diffs,[FP_areas.area[FP_areas.index[-1]]*0.9069/1]*1)  # add in area for last 1 rows.
  
  ##################
 
  counts = area_diffs/px_areas # pixels which fit per region if packing is perfect. 
  px_count = np.floor(counts)    # (rounded down)

  change = True

  #  some way to re-arrange pixel counts
  while change:
    used_areas = px_count*px_areas
    print_status(pixel_types, px_count, px_areas, used_areas, area_diffs)
   
    print '\nChange the px_count variable to change pixel numbers in a given band.\nThen continue (c) the program.\n'+\
          'If happy with distribution set variable change=False.\n'
  
    check_areas(px_count,px_areas,area_diffs)
 
    # write px_count to full list corresponding to all bands.
    for i,pixel in enumerate(pixel_types):
      px_index = bands.index[bands.pixel==pixel]
      results.num_bolo[px_index] = px_count[i] * 2

    results = calc_array_noise(results)

    pdb.set_trace()

 
elif settings.calc_N_px_rough: ## estimate N_px that may fit.
  print 'Estimating number of pixel which fit in focal plane.\n'
  # usable focal plane area assuming hex pack = 0.9069
  fp_area = np.pi * settings.x_fp_diameter * settings.x_fp_diameter / 4. * 0.9069

  tot_area = np.sum(px_count*px_areas)

  # scale detector count to fill focal plane.
  px_count = np.floor(fp_area/tot_area * px_count)
  # area_out = np.sum(px_count*px_areas)
  
  # write px_count to full list corresponding to all bands.
  for i,pixel in enumerate(pixel_types):
    px_index = bands.index[bands.pixel==pixel]
    results.num_bolo[px_index] = px_count[i] * 2

  results = calc_array_noise(results)

 
else: # just use N_px from bands.csv.
  print 'Using pixel count from %s.' %(settings.bands_path.split('/')[-1])

  px_count = np.zeros(len(pixel_types))
  px_areas = np.zeros(len(pixel_types))
  for i,pixel in enumerate(pixel_types):
    px_indices = np.where(bands.pixel==pixel)[0]
    ## just use first entry for D_px, all should be equal. then get area
    px_areas[i] = np.pi*results.D_px[px_indices[0]]**2. / 4. 
    px_count[i] = bands.number[px_indices[0]]

  ## get availible areas.
  FP_areas = pandas.read_csv(os.path.join(settings.FP_areas_path))
  area_diffs = -np.diff(FP_areas.area) * 0.9069   # area for each frequency. include hex pack factor
  
  ######## 3 single band pixels hardcoded!!!  ###############
  # for G,H,I.  do even thirds at smallest area.
  # area_diffs = np.append(area_diffs,[FP_areas.area[FP_areas.index[-1]]*0.9069/3]*3)  # add in area for last 3 rows.
  # area_diffs = np.append(area_diffs,[FP_areas.area[FP_areas.index[-1]]*0.9069/3]*2)  # add in area for last 2 rows.
  area_diffs = np.append(area_diffs,[FP_areas.area[FP_areas.index[-1]]*0.9069/1]*1)  # add in area for last 1 rows.
    
  used_areas = px_count*px_areas # calc areas to check fit.
  
  if any(used_areas > area_diffs):
    print '\n\t\tPossible warning!!\n Too many pixel in Band(s) %s\n' %', '.join(pixel_types[used_areas > area_diffs])
  change=True
  #else:  Just always force user input.
  #  change = False

  while change:
    used_areas = px_count*px_areas
    print_status(pixel_types, px_count, px_areas, used_areas, area_diffs)
   
    check_areas(px_count,px_areas,area_diffs)

    print '\nChange the px_count variable to change pixel numbers in a given band.\nThen continue (c) the program.\n'+\
          'If happy with distribution set variable change=False.\n'

    for i,pixel in enumerate(pixel_types):
      px_index = bands.index[bands.pixel==pixel]
      results.num_bolo[px_index] = px_count[i] * 2
    
    results = calc_array_noise(results)
    #change=False
    pdb.set_trace()    
    

## resave the thingies.
results.to_csv(os.path.join(load_path, 'All_results_out.csv'), index=True)

#pdb.set_trace()
