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
usage = "usage: %prog  <input_csv1> <input_csv2> . . . \n Does csv1 / csv2 in whatever params you want"
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
  return ax

def plot_col(df,col,ax,scale=1e12,**kwargs):
  ax.plot(df.nu,df[col]*scale,**kwargs)
  return ax

# load data.
data1 = load_NEPs_data(args[0]) 
data2 = load_NEPs_data(args[1])

#make cool plots.
#create a plot, x-axis frequency
fig,ax = plt.subplots()
ax.set_xlabel('Frequency, GHz')
ax.set_ylim([.3,1.])

plot_col_ratio(data1,data2,'total_pow',ax, marker='o',label='optical power')
ax.set_ylabel('Open/Crossed')

fig2,ax2 = plt.subplots()
ax2.set_xlabel('Frequency, GHz')
ax2.set_ylim([.3,1.])

for col in ['NEP_poisson', 'NEP_photon', 'NEP_phonon','NEP_johnson', 'NEP_readout']:
  plot_col_ratio(data1,data2,col,ax2,marker='.',ls='-',label=col)

ax2.set_ylabel('NEP ratios, Open/Crossed') 
ax2.legend(prop={'size':10})

#fig3,ax3 = plt.subplots()
#ax3.set_xlabel('Frequency, GHz')
#ax3.set_ylim([.3,1.])

plot_col_ratio(data1,data2,'NET_total',ax,marker='o',label='NET')
#ax3.set_ylabel('NET ratio')

ax.legend(prop={'size':13})

fig.savefig('./outputs/plots/Pow_NET_ratio.png')
fig2.savefig('./outputs/plots/NEP_ratios.png')
#fig3.savefig('./outputs/plots/NET_ratios.png')


fig4,[ax4a,ax4b] = plt.subplots(2,1,sharex=True)
ax4a.set_title('Open Dragone')
ax4b.set_xlabel('Frequency, GHz')

plot_col(data1,'total_pow',ax4a, scale=1e12,marker='o',label='optical power')
ax4a.set_ylabel('P_opt, Pw')
plot_col(data1,'NET_total',ax4b, scale=1e6,marker='o',label='Total NET')
ax4b.set_yscale('log')
ax4b.set_ylabel('NET, uK rt(sec)')



fig5,ax5 = plt.subplots(1)
ax5.set_title('Open Dragone')
ax5.set_xlabel('Frequency, GHz')

for col in ['NEP_poisson', 'NEP_photon', 'NEP_phonon','NEP_johnson', 'NEP_readout']:
  plot_col(data1,col,ax5, scale=1e18,marker='o',label=col)
ax5.set_ylabel('NEPs, Aw/rt(Hz)')
ax5.legend(prop={'size':10})

fig4.savefig('./outputs/plots/Open_Popt_NET.png')
fig5.savefig('./outputs/plots/Open_NEP.png')

plt.show()







