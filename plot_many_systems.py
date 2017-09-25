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
import pdb

matplotlib.rcParams.update({'font.size': 18}) # set good font sizes

# optparse it!
usage = "usage: %prog  <input_csv1> <input_csv2> . . . \n plots various params of all csvs."
parser = optparse.OptionParser(usage)
#parser.add_option('-f', '--file', dest='measured_data_path', action='store', type='str')
#parser.add_option('-i', '--import', dest='import_T', action='store_true', 
#                  default=False, help='If included program looks for'
#                  'transmission data at paths, arg1, arg2 ...')
#parser.add_option('-i', '--import', dest='import_num', action='store', default=0,type='int')
#parser.add_option('-s', '--shift', dest='shift', action='store', 
#                  default=0,type='float',
#                  help='Adds a shift in frequency to the transmission data')
(option, args) = parser.parse_args()

dirs = args

#functions:
#load results df.
def load_NEPs_data(directory):
  df = pandas.read_csv(os.path.join(directory,'All_results_out.csv'))
  return df

def load_bolo_char(directory):
  df = pandas.read_csv(os.path.join(directory,'Bolo_char_out.csv'))
  return df

#add a column.  any grabbed column. choose color, style, label
def plot_col_ratio(df1,df2,col,ax,**kwargs):
  # dataframe, chosen col name, axis handle, kwargs for plot
  ax.plot(df1.nu,df1[col]/df2[col],**kwargs)
  #print col,df1[col]/df2[col]
  return ax

def plot_col(df,col,ax,scale=1e12,**kwargs):
  ax.plot(df.nu,df[col]*scale,**kwargs)
  return ax

## get names
names = ['_'.join(os.path.basename(os.path.normpath(i)).split('_')[0:2]) for i in dirs]

## load all data
data={}
for i,path in enumerate(dirs):
  data[names[i]]=load_NEPs_data(path)

names.sort(key=lambda i: int(i.split('K')[0])) # have to sort after loading data.

pdb.set_trace()

## plot power
fig1,ax1 = plt.subplots(1)
ax1.set_xlabel('Frequency, GHz')
ax1.set_ylabel('Optical power, pW')

for name in names:
  plot_col(data[name], 'total_pow',ax1,scale=1e12,marker='.',label=name)

ax1.legend(prop={'size':10})
fig1.tight_layout()
fig1.savefig('./outputs/plots/multi_system_Popt.png')

## plot NET
fig2,ax2 = plt.subplots(1)
ax2.set_xlabel('Frequency, GHz')
ax2.set_ylabel('NET, uK*rt(s)')
ax2.set_yscale('log')

for i in names:
  plot_col(data[i], 'NET_total',ax2,scale=1e6,marker='.',label=i)

ax2.legend(prop={'size':10})
fig2.tight_layout()
fig2.savefig('./outputs/plots/multi_system_NET.png')

## plot NEP
fig3,ax3 = plt.subplots(1)
ax3.set_xlabel('Frequency, GHz')
ax3.set_ylabel('total NEP, aW/rt(Hz)')

for i in names:
  plot_col(data[i], 'NEP_total',ax3,scale=1e18,marker='.',label=i)

ax3.legend(prop={'size':10})
fig3.tight_layout()
fig3.savefig('./outputs/plots/multi_system_NEP.png')


plt.show()

