# quick rough code to scrape data out of .csv files 
# and return various noise averages. to give a 1 number metric.

import os
import pandas
import numpy as np
import optparse
import pdb

from default_settings import *

# optparse it!
usage = "usage: %prog -f <folder_to_read_data_from>"
parser = optparse.OptionParser(usage)
parser.add_option('-f', dest='folder', action='store', type='str', default=None,
                  help='folder which all results_df is loaded from.  Finds All_results_out.csv in this folder.')
(option, args) = parser.parse_args()
#load_path = os.path.join(settings.base_path,'%s' %(option.folder))
load_path = option.folder

#load data
results = pandas.read_csv(os.path.join(load_path, 'All_results_out.csv'),index_col=0)

def total_CMB_sense(results):
  # takes px number and NET per pixel.
  # gives final CMB sensitivity of current system.
      
  # sum 1/squares of noise, take sqrt, invert.
  total_pol = 1./np.sqrt(np.sum(np.power(results.pol_weight,-2.)))
  total_corrPol = 1./np.sqrt(np.sum(np.power(results.corr_pol_weight,-2.)))
  print 'total uK arcmin uncorrelate:  %.2f' %(total_pol*1.e6)
  print 'total uK arcmin correlated :  %.2f' %(total_corrPol*1.e6)
  return total_pol, total_corrPol


def hifreq_CMB_sense(results):
  # takes px number and NET per pixel.
  # gives final CMB sensitivity of current system.
      
  # sum 1/squares of noise, take sqrt, invert.
  total_pol = 1./np.sqrt(np.sum(np.power(results.pol_weight[-6:],-2.)))
  total_corrPol = 1./np.sqrt(np.sum(np.power(results.corr_pol_weight[-6:],-2.)))
  print 'total uK arcmin uncorrelate:  %.2f' %(total_pol*1.e6)
  print 'total uK arcmin correlated :  %.2f' %(total_corrPol*1.e6)
  return total_pol, total_corrPol


def lofreq_CMB_sense(results):
  # takes px number and NET per pixel.
  # gives final CMB sensitivity of current system.
      
  # sum 1/squares of noise, take sqrt, invert.
  total_pol = 1./np.sqrt(np.sum(np.power(results.pol_weight[:6],-2.)))
  total_corrPol = 1./np.sqrt(np.sum(np.power(results.corr_pol_weight[:6],-2.)))
  print 'total uK arcmin uncorrelate:  %.2f' %(total_pol*1.e6)
  print 'total uK arcmin correlated :  %.2f' %(total_corrPol*1.e6)
  return total_pol, total_corrPol

print 'All bands'
total_CMB_sense(results)

print '6 highest bands'
hifreq_CMB_sense(results)

print '6 lowest bands'
lofreq_CMB_sense(results)
