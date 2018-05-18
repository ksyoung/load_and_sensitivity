import os, sys
import matplotlib
import pandas
#import time
import numpy as np
#import scipy.integrate as scint
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import optparse
import pdb

matplotlib.rcParams.update({'font.size': 28}) # set good font sizes

# optparse it!
usage = "usage: %prog  -b <bands.csv>"
parser = optparse.OptionParser(usage)
parser.add_option('-b', dest='bands', action='store', type='str', default=None, 
                  help='csv file which defines bands and pixel labels.') 
parser.add_option('-n', dest='height_is_num', action='store_true', default=False, 
                  help='T/F to use pixel counts to set y axis height of bars. default False') 

(option, args) = parser.parse_args()

# load bands
bands = pandas.read_csv(option.bands)

# color vector
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

# matching SW focal plane: 
colors = ['#c0c0c0', '#00668f', '#49ff37', '#67cdff', '#0003ff', '#5c57ff', '#ff47f4', '#ff47f4', '#ff47f4']

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



# initialize plot
fig, ax = plt.subplots(figsize=(16,9))


# sort out pixels
pixel_types = []
for i in bands.loc[:,('pixel')]:
  if i not in pixel_types:
    pixel_types.append(i)
pixel_types = np.array(pixel_types) ## later code works if it's an np array.

height_is_num = True

if option.height_is_num:
  for i,pixel in enumerate(pixel_types):
    if pixel in(['G','H','I']):
      ax.bar(bands.nu_low[bands.pixel==pixel], bands.number[bands.pixel==pixel]*2., bands.width[bands.pixel==pixel], align='edge', color=colors[i], alpha=.8, ec='k')
    else:
      px_band = np.array([bands.nu_low[bands.pixel==pixel],bands.nu_high[bands.pixel==pixel]]).T # a nx2 array of nu_low, nu_high for the n bands in the pixel.

      ax.bar(px_band[:,0], bands.number[bands.pixel==pixel]*2., px_band[:,1]-px_band[:,0], align='edge', fc=colors[i], alpha=.8, ec='k', fill=True)

    print 'done with: %s' %pixel
  ax.set_ylabel('Number of Detectors ')

else:
  for i,pixel in enumerate(pixel_types):
    if pixel in(['G','H','I']):
      ax.bar(bands.nu_low[bands.pixel==pixel], 1./(i%2/10.+1), bands.width[bands.pixel==pixel], align='edge', color=colors[i], alpha=.7, ec='k')
    else:
      px_band = np.array([bands.nu_low[bands.pixel==pixel],bands.nu_high[bands.pixel==pixel]]).T # a nx2 array of nu_low, nu_high for the n bands in the pixel.

      ax.bar(px_band[:,0], [1./(i%2/10.+1)]*len(px_band), px_band[:,1]-px_band[:,0], align='edge', fc=colors[i], alpha=.7, ec='k', fill=True)

    print 'done with: %s' %pixel
  ax.set_ylabel('Arbitrary')

ax.set_xscale('log')
ax.set_xticks(np.round(bands.nu[::2]))
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
plt.tick_params(axis='x', which='minor',bottom='off')
plt.tick_params(axis='y', which='minor',left='off')

# hide axis tick labels
#ax.axes.yaxis.set_ticklabels([])

#ax.set_ylim([0,1075])
ax.set_xlabel('Frequency, GHz')
fig.tight_layout()

fig.savefig('./outputs/plots/bands.png')
fig.savefig('./outputs/plots/bands.tif',format='tiff', dpi=1200)
plt.show()

# 
# write px_count to full list corresponding to all bands.
for i,pixel in enumerate(pixel_types):
  px_index = bands.index[bands.pixel==pixel]
  #results.num_bolo[px_index] = px_count[i] * 2


# other option
# just plot all is sets of 3.
#for i in range(7):
#  make_and_plot_MS_flambda(bands[i*3:((i+1)*3)])

# plot actual pixels.
pixel_bands = [bands[0:5:2],bands[1:6:2],bands[6:11:2],bands[7:12:2],bands[12:18:2],bands[13:19:2],bands[18:]]
#for i in range(7):
# need new function. make_and_plot_MS_Dpx(pixel_bands[i])

  
pdb.set_trace()
