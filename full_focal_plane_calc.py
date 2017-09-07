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


# load bands (has px numbers)
bands = pandas.read_csv(settings.bands_path)
# path to newest data.  May need to change this.......
print 'Loading newest data for %s\n' %settings.version
load_path = os.path.join(settings.base_path,'outputs/%s/current_directory' %settings.version)
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
def check_area(fp_area,px_area):
  # both arrays,same lenght please.
  pass
  # see if anything fails
  # print failures
 
  # print extra space.
  return




if settings.calc_N_px_by_area_csv:
  print 'Estimating number of pixel which fit in focal plane.\n'

  FP_areas = pandas.read_csv(os.path.join(settings.FP_areas_path))
  area_diffs = -np.diff(FP_areas.area) * 0.9069   #extra area for each frequency. include hex pack factor
  area_diffs = np.append(area_diffs,FP_areas.area[FP_areas.index[-1]]*0.9069)  # add in last row.

  # how to deal with single band pixels?!?!??!?!?!?!?
  # not using most of the area anyway......
  # for G,H,I.  do even thirds at smallest area.

   
  pixel_types = []
  for i in bands.loc[:,('pixel')]:
    if i not in pixel_types:
      pixel_types.append(i)

  px_count = np.zeros(len(pixel_types))
  px_areas = np.zeros(len(pixel_types))
  for i,pixel in enumerate(pixel_types):
    px_indices = np.where(bands.pixel==pixel)[0]

    ## just use first entry for D_px, all should be equal. then get area
    px_areas[i] = np.pi*results.D_px[px_indices[0]]**2. / 4. 

  px_areas_2 = np.append(px_areas[0:6],np.mean(px_areas[-3:]))
  counts = area_diffs/px_areas_2
  # split last set amongst 3 bands.  
  px_count = np.floor(np.append(counts[0:-1],[counts[-1]/3.,counts[-1]/3.,counts[-1]/3.]))

  # is there some way to re-arrange pixel counts?????
  # check if area works (function)
  # print extra areas (function)
  # change values and redo (while loop with pdb inside??)
  change = True
  while change:
    print 'Current pixel counts and areas:\n'
    print tabulate([pixel_types, px_counts,px_areas,area_diffs], 
                   headers=['pixel type', 'count','area used', 'area availible']) 

    print 'Change the px_count variable to change pixel numbers in a given band.\nThen continue (c) the program.'
    pdb.set_trace()

  # write px_count to full list corresponding to all bands.
  for i,pixel in enumerate(pixel_types):
    px_index = bands.index[bands.pixel==pixel]
    results.num_bolo[px_index] = px_count[i] * 2

  pdb.set_trace()

elif settings.calc_N_px_rough: ## estimate N_px that may fit.
  print 'Estimating number of pixel which fit in focal plane.\n'
  # usable focal plane area assuming hex pack = 0.9069
  fp_area = np.pi * settings.x_fp_diameter * settings.x_fp_diameter / 4. * 0.9069

  pixel_types = []
  for i in bands.loc[:,('pixel')]:
    if i not in pixel_types:
      pixel_types.append(i)

  px_count = np.zeros(len(pixel_types))
  px_areas = np.zeros(len(pixel_types))
  for i,pixel in enumerate(pixel_types):
    px_indices = np.where(bands.pixel==pixel)[0]

    px_count[i] = np.max(bands.number[px_indices])
    ## just use first entry for D_px, all should be equal. then get area
    px_areas[i] = np.pi*results.D_px[px_indices[0]]**2. / 4. 

  tot_area = np.sum(px_count*px_areas)

  # scale detector count to fill focal plane.
  px_count = np.floor(fp_area/tot_area * px_count)
  # area_out = np.sum(px_count*px_areas)
  
  # write px_count to full list corresponding to all bands.
  for i,pixel in enumerate(pixel_types):
    px_index = bands.index[bands.pixel==pixel]
    results.num_bolo[px_index] = px_count[i] * 2

 
else: # just use N_px from bands.csv.
  print 'Using pixel count from %s.' %(settings.bands_path.split('/')[-1])
  results.num_bolo = bands.number * 2


# calculate per band sensitivity, mapping speed, and pol weight
results.NET_array = results.NET_total / np.sqrt(results.num_bolo)
results.map_speed = ((results.FWHM / 60.) / results.NET_array)**2.  ## FWHM**2 / NET **2 , deg^2 / (K^2 sec)
results.weight = results.NET_array * np.sqrt(settings.sky_area * 3600. / (settings.mission_length * 3.154e7)) ## array_NET * sqrt(sky_area / time) K * arcmin ?? may be wrong???
results.pol_weight = results.weight * np.sqrt(2.) ## ?? may be wrong??


## resave the thingies.
results.to_csv(os.path.join(load_path, 'All_results_out.csv'), index=True)

pdb.set_trace()
