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

matplotlib.rcParams.update({'font.size': 16}) # set good font sizes

# optparse it!
usage = "usage: %prog  <input_csv1> <input_csv2> . . . \n Does csv1 / csv2 in whatever params you want"
parser = optparse.OptionParser(usage)
parser.add_option('--t1', dest='title1', action='store', type='str', default='Open Dragone', 
                  help='Title for first set of data.')
parser.add_option('--t2', dest='title2', action='store', type='str', default='Crossed Dragone',
                  help='Title for second set of data.')
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

# load data.
data1 = load_NEPs_data(args[0]) 
data2 = load_NEPs_data(args[1])

#make cool plots.
#create a plot, x-axis frequency
fig,ax = plt.subplots()
ax.set_xlabel('Frequency, GHz')
#ax.set_ylim([.3,1.])

plot_col_ratio(data1,data2,'NET_total',ax,marker='o',label='NET')
plot_col_ratio(data1,data2,'total_pow',ax, marker='o',label='optical power')
#ax.plot(data1.nu,np.sqrt(data1['total_pow']/data2['total_pow']),marker='x',label='sqrt(optical power)')

ax.set_ylabel(option.title1+' / '+option.title2)
#ax.set_ylabel('250 mK / 100 mK')

fig2,ax2 = plt.subplots()
ax2.set_xlabel('Frequency, GHz')
#ax2.set_ylim([.3,1.])

for col in ['NEP_poisson', 'NEP_photon', 'NEP_phonon','NEP_johnson', 'NEP_readout']:
  plot_col_ratio(data1,data2,col,ax2,marker='.',ls='-',label=col)

ax2.set_ylabel('NEP ratios,\n'+ option.title1 + ' / ' + option.title2) 
#ax2.set_ylabel('NEP ratios, 250 mK / 100 mK') 
ax2.legend(prop={'size':10})

#fig3,ax3 = plt.subplots()
#ax3.set_xlabel('Frequency, GHz')
#ax3.set_ylim([.3,1.])

#ax3.set_ylabel('NET ratio')

ax.legend(prop={'size':13})

#pdb.set_trace()
fig.tight_layout()
fig2.tight_layout()
fig.savefig('./outputs/plots/Pow_NET_ratio.png')
fig2.savefig('./outputs/plots/NEP_ratios.png')
#fig3.savefig('./outputs/plots/NET_ratios.png')

## plot for each case

#fig4,[[ax4a,ax4b],[ax4c,ax4d]] = plt.subplots(2,2,sharex='col',sharey='row',figsize=(9.6,7.2))
fig4,[ax4a,ax4b] = plt.subplots(1,2,sharex='row',figsize=(9.6,4.8))
#ax4a.set_title('Open Dragone')
#ax4a.set_title('Tbath 250 mK')
#ax4b.set_title('Crossed Dragone')
#ax4b.set_title('Tbath 100 mK')

ax4a.set_xlabel('Frequency, GHz')
ax4b.set_xlabel('Frequency, GHz')

ax4a.set_ylabel('P_opt, Pw')
ax4b.set_ylabel('NET, uK rt(sec)')
ax4b.set_yscale('log')

# data 1
plot_col(data1,'total_pow',ax4a, scale=1e12,marker='.',label=option.title1)
plot_col(data1,'NET_total',ax4b, scale=1e6,marker='.',label=option.title1)

# data 2
plot_col(data2,'total_pow',ax4a, scale=1e12,marker='.',label=option.title2)
plot_col(data2,'NET_total',ax4b, scale=1e6,marker='.',label=option.title2)

ax4a.legend(prop={'size':12})
ax4b.legend(prop={'size':12})

fig5,[ax5a,ax5b] = plt.subplots(1,2,sharex=True, sharey=True, figsize=(9.6,4.8))
ax5a.set_title(option.title1)
#ax5a.set_title('T_bath 250 mK')
ax5b.set_title(option.title2)
#ax5b.set_title('T_bath 100 mK')

ax5a.set_xlabel('Frequency, GHz')
ax5b.set_xlabel('Frequency, GHz')
ax5a.set_ylabel('NEPs, aW/rt(Hz)')

for col in ['NEP_poisson', 'NEP_photon', 'NEP_phonon','NEP_johnson', 'NEP_readout']: #
  plot_col(data1,col,ax5a, scale=1e18,marker='.',label=col)
  plot_col(data2,col,ax5b, scale=1e18,marker='.',label=col)

#ax5a_sub = ax5a.axis([.1,.7,.3,.3], axisbg='y')
#ax5b_sub = ax5b.axis([.1,.7,.3,.3], axisbg='y')
#pdb.set_trace()
#ax5a_sub.set_xlim([15,60])
#ax5b_sub.set_xlim([15,60])
#for col in ['NEP_poisson', 'NEP_photon', 'NEP_phonon','NEP_johnson', 'NEP_readout']:
#  plot_col(data1,col,ax5a_sub, scale=1e18,marker='.',label=col)
#  plot_col(data2,col,ax5b_sub, scale=1e18,marker='.',label=col)


ax5a.legend(prop={'size':10})
ax5b.legend(prop={'size':10})

fig6,ax6 = plt.subplots(1,sharex=True, sharey=True, figsize=(6.4,4.8))
#ax6a.set_title('Open Dragone')
#ax6b.set_title('Crossed Dragone')

ax6.set_xlabel('Frequency, GHz')
ax6.set_ylabel('FWHM, arcmin')

plot_col(data1,'FWHM',ax6, scale=1,marker='o',label=option.title1)
plot_col(data2,'FWHM',ax6, scale=1,marker='o',label=option.title2)

ax6.legend(prop={'size':12})

fig4.tight_layout()
fig5.tight_layout()
fig6.tight_layout()
fig4.savefig('./outputs/plots/Both_systems_Popt_NET.png')
fig5.savefig('./outputs/plots/Both_systems_NEP.png')

ax5a.set_xlim([15,100])
ax5a.set_ylim([0,10])
fig5.tight_layout()
fig5.savefig('./outputs/plots/Both_systems_NEP_zoom.png')

fig6.savefig('./outputs/plots/Both_systems_FWHM.png')



try:
  fig7,[ax7a,ax7b] = plt.subplots(1,2,sharex='row',figsize=(9.6,4.8))

  ax7a.set_xlabel('Frequency, GHz')
  ax7b.set_xlabel('Frequency, GHz')

  ax7b.set_ylabel('Polarization weight, uK arcmin')
  ax7a.set_ylabel('Polarization weight ratio \n' + option.title1 + ' / ' + option.title2)
  #ax7a.set_ylabel('Polarization weight ratio \n Open/Coma corrected')
  ax7b.set_yscale('log')

  plot_col_ratio(data1,data2,'pol_weight',ax7a, marker='.',label='Pol weight')

  # data 1
  plot_col(data1,'pol_weight',ax7b, scale=1e6,marker='.',label=option.title1)

  # data 2
  #plot_col(data2,'total_pow',ax4a, scale=1e12,marker='.',label='Cross Dragone')
  plot_col(data2,'pol_weight',ax7b, scale=1e6,marker='.',label=option.title2)
  #plot_col(data2,'pol_weight',ax7b, scale=1e6,marker='.',label='Coma corrected')

  #ax7a.legend(prop={'size':12})
  ax7b.legend(prop={'size':12})

  fig7.tight_layout()
  fig7.savefig('./outputs/plots/Both_systems_pol_weight.png')

except:
  pass




pdb.set_trace()

plt.show()

SMALL_SIZE = 26
MEDIUM_SIZE = 32
BIGGER_SIZE = 44

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rc('axes',linewidth=2)

#fig4,[[ax4a,ax4b],[ax4c,ax4d]] = plt.subplots(2,2,sharex='col',sharey='row',figsize=(9.6,7.2))
fig4,ax4b = plt.subplots(figsize=(16,8))
#ax4a.set_title('Open Dragone')
#ax4a.set_title('Tbath 250 mK')
#ax4b.set_title('Crossed Dragone')
#ax4b.set_title('Tbath 100 mK')

ax4b.set_xlabel('Frequency, GHz')

ax4b.set_ylabel(r'NET, $\mu K_{CMB}$ arcmin')
ax4b.set_yscale('log')

# data 1
plot_col(data1,'corr_pol_weight',ax4b, scale=1e6,marker='o',lw=3,markersize=10,label=option.title1)

fig4.tight_layout()
fig4.savefig('./outputs/plots/pol_weight_log_big.png', dpi=600)
plt.show()





