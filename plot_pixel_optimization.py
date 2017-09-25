# -*- coding: utf-8 -*-
'''
Code to read in .xlsx file of data at many edge tapers and plot 
mapping speed as function of pixel size.


Karl Young, Sept. 22 2017

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
import pdb

matplotlib.rcParams.update({'font.size': 18}) # set good font sizes

# optparse it!
usage = "usage: %prog  -n name  <input_xlsx> \n plot various params of single .xlsx."
parser = optparse.OptionParser(usage)
parser.add_option('-n', '--name', dest='name', action='store', type='str', default=None)
parser.add_option('-f', '--fnum', dest='fnum', action='store', type='float', 
                  default=0.0, help='f number of current system.  Critical for plotting!!')

#parser.add_option('-i', '--import', dest='import_num', action='store', default=0,type='int')
#parser.add_option('-s', '--shift', dest='shift', action='store', 
#                  default=0,type='float',
#                  help='Adds a shift in frequency to the transmission data')
(option, args) = parser.parse_args()

dirs = args[0]

name = option.name
f_number = option.fnum
if not f_number > 0.0:
  print '\nEnter system fnumber with -f flag!  Please.\n'
  sys.exit()

# load xls.
data_dict = pandas.read_excel(dirs,sheetname=None)
keys = data_dict.keys()
keys.sort(key=float)
index = data_dict[keys[0]].index

#want at 1 frequency, dict with key=freq, list is all pixel sizes.
edge_db = {}
D_px_flambda = {}  # have to hand fnumber to this.  possibility for lots of errors... get it from settings??
D_px = {}
NET_arr = {}
corr_NET_arr = {}
corrCMB_NET_arr = {}
bands = []

for i in index:  ## step through bands.
  for key in keys:
    df = data_dict[key]
    freq = df.nu[i]
    edge_db.setdefault(freq,[]).append(df.edge_dB[i])  # odd method found on stack overflow.
    D_px.setdefault(freq,[]).append(df.D_px[i]*1e3)
    D_px_flambda.setdefault(freq,[]).append(df.D_px[i]/(f_number*0.299792458/freq))
    NET_arr.setdefault(freq,[]).append(df.NET_array[i]*1e6)
    corr_NET_arr.setdefault(freq,[]).append(df.corr_NET_array[i]*1e6)
    corrCMB_NET_arr.setdefault(freq,[]).append(df.corrCMB_NET_array[i]*1e6)
  bands.append(freq)

# plot MS vs px pitch, D
 # for base case, corr CMB, corr all.
 # probably do is some sort of loop.

def make_and_plot_MS(band_arr):
  fig,ax = plt.subplots()
  ax.set_xlabel('Pixel size, F-lambda')
  # can add edge taper on top axis??
  ax.set_ylabel(r'1/NET**2, $(uK^2/\sqrt{s})^{-1}$')
  ax.set_title(name)

  colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']*5
  for i,band in enumerate(band_arr):
    ax.plot(D_px_flambda[band],1/np.array(NET_arr[band])**2,marker='.',ls='--',alpha=.5,
            color=colors[i], label='%.1f GHz, no correlation' %band)
    ax.plot(D_px_flambda[band],1/np.array(corrCMB_NET_arr[band])**2,marker='.',ls='-',alpha=.5,
            color=colors[i], label='%.1f GHz, CMB correlated' %band)
    ax.plot(D_px_flambda[band],1/np.array(corr_NET_arr[band])**2,marker='.',ls='--',alpha=1,
            color=colors[i], label='%.1f GHz, all correlated' %band)

  ax.legend(prop={'size':8})
  fig.tight_layout()
  fig.savefig('./outputs/plots/MS_vs_px_%s.png' %'_'.join([str(i) for i in np.round(band_arr,1)]))
  plt.show()

  return

for i in range(7):
  make_and_plot_MS(bands[i*3:((i+1)*3)])

pdb.set_trace()

