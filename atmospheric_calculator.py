'''
Original code by Ben Westbrook, Berkeley sometime in 2013-2014.

Lines approximately 90-105 changed by Karl Young, Minnesota Sept 2014. 
Changed to using intensity output by am code (code from CFA website) instead of
going through planck temp as output.  Should work the same, but previously I was getting
number about 10x too high for 150 GHz band.

'''

#!/usr/bin/env python
import os
import pandas
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt


def load_single_file(data_path):
    return pandas.read_csv(data_path, sep=' ')

def load_transmission_data(data_frame):
    return data_frame.ix[:,0], data_frame.ix[:,3].values

def load_band_information():
    band_directory_path = '/home/westbrook/ebex/branches/leap/resources/bands/'
    bands_dict = {}
    for frequency in [150, 250, 410]:
        full_path = os.path.join(band_directory_path, str(frequency) + 'band_avespect_bin.csv')
        band_data_frame = pandas.read_csv(full_path, sep=',')
        bands_dict[str(frequency)] = [band_data_frame.ix[:,0].values, band_data_frame.ix[:,1].values]
    return bands_dict
'''
def band_vs_transmission_plotter():
    color_dict = {'150': 'r', '250': 'g', '410': 'b'}
    fig = plt.figure(figsize=(9.75, 7.5), dpi=75)
    ax = fig.add_axes([0.13, 0.23, 0.75, 0.70])
    full_path = '/home/westbrook/Python/ebex_tools/ebex6k_mapping_speed/atmospheric_transmission/simulated_atm_models/LDB_at_34p61km.out'
    data_frame = load_single_file(full_path)
    bands_dict = load_band_information()
    frequency_, transmission_ = load_transmission_data(data_frame)
    ax.plot(frequency_, transmission_, 'k', alpha=0.4)
    for frequency, data in bands_dict.iteritems():
        frequency_vector, band_pass_vector = data[0], data[1]
        label = '%s GHz band pass'  % frequency
        ax.plot(frequency_vector, band_pass_vector, color=color_dict[frequency], label=label, lw=3.5)
    ax.set_ylim([0.1, 1.01])
    ax.set_xlim([0., 500.0])
    ax.set_ylabel('Normalized Atm Tran, Band Passes', fontsize=20)
    ax.set_xlabel('Frequency (GHz)', fontsize=20)
    ax.set_title('Atmospheric Transmission and EBEX Band Passes', fontsize=17)
    #ax.legend(loc=8, borderaxespad=0., numpoints=1, fontsize=16)
    ax.legend(bbox_to_anchor=(0.500, -0.3), loc=8, borderaxespad=0., numpoints=1, fontsize=16)
    plt.savefig('/home/westbrook/Dropbox/Thesis/Figures/ebex_bands_vs_atmospheric_transmission.pdf')
    fig.show()
    import ipdb;ipdb.set_trace()

def transmission_plotter(raw_file_name_dict, frequency_low=0., frequency_high=500.):
    location_label_dict = {'Chajnantor': 'Chajnantor, 60 deg, 1.0mm (pwv)',
                          'LDB': 'LDB, 60 deg, 34.61 km (altitude)',
                          'SouthPole': 'South Pole, 60 deg, 1.0mm (pwv)'}

    fig = plt.figure(figsize=(9.75, 7.5), dpi=75)
    ax = fig.add_axes([0.13, 0.23, 0.75, 0.70])
    for raw_file_name, plot_color in raw_file_name_dict.iteritems():
        full_path = os.path.join('/home/westbrook/Python/ebex_tools/ebex6k_mapping_speed/atmospheric_transmission/simulated_atm_models',
                                 raw_file_name)
        data_frame = load_single_file(full_path)
        frequency_vector, transmission_vector = load_transmission_data(data_frame)
        #normalized_transmission_vector = transmission_vector / max(transmission_vector)
        location = raw_file_name.split('_')[0]
        label = location_label_dict[location]
        ax.plot(frequency_vector, transmission_vector, color=plot_color, alpha=0.6, label=label)
    ax.set_ylim([0.1, 1.01])
    ax.set_xlim([0., 500.0])
    ax.set_ylabel('Atmospheric Transmission', fontsize=20)
    ax.set_xlabel('Frequency (GHz)', fontsize=20)
    ax.set_title('Atmospheric Transmission for the South Pole, McMurdo LDB, and Chajnantor', fontsize=17)
    ax.legend(bbox_to_anchor=(0.81, -0.20), loc=5, borderaxespad=0., numpoints=1, fontsize=16)
    plt.savefig('/home/westbrook/Dropbox/Thesis/Figures/Atmospheric_Comparison.pdf')
    fig.show()
    import ipdb;ipdb.set_trace()
'''
def atmospheric_calculator(atmospheric_window_data_path, frequency_low, frequency_high):
    '''
    This function returns the effetive emissivity and transmition of the atmpohsere
    given an atmpospheric transmission data file and a frequency range
    '''

    #---CONSTANT---------------------
    h_ = 6.626e-34 # %planck's constant
    k_b = 1.38e-23 # %boltzmann constant J/ K
    c_ = 3e8 # speed of light meters / s

    # Load in the atmospheric transimissino data file
    atmosphere_data_frame = load_single_file(atmospheric_window_data_path)
    frequency = atmosphere_data_frame.iloc[:,0] * 1e9  # in Hz
    optical_depth = atmosphere_data_frame.iloc[:,1]
    planck_temp = atmosphere_data_frame.iloc[:,2]
    transmission = atmosphere_data_frame.iloc[:,3]
    intensity = atmosphere_data_frame.iloc[:,4]  #adjusting units? testing...  #check units!  output of am code should be, watt*m-2*Hz-1*sr-1
    '''
    Ben's old code.  File seems to be wrong format for this to work.
    frequency = atmosphere_data_frame['Frequency'].values * 1e9  # in Hz
    optical_depth = atmosphere_data_frame['Optical_Depth'].values
    planck_temp = atmosphere_data_frame['Planck_Temp'].values
    transmission = atmosphere_data_frame['Transmission'].values
    '''
    #black_body_ATM = 2*(h_ * frequency ** 3 / (c_ ** 2)) * (1 / (np.e ** ((h_ * frequency) / (planck_temp * k_b)) - 1))
    black_body_ATM2 = intensity  #just use intensity from am code, compare to B_nu(277), get effective emiss.
    black_body_277 = 2*(h_ * frequency ** 3 / (c_ ** 2)) * (1/ (np.e ** ((h_ * frequency) / (277 * k_b)) - 1))

    #print 'max of intensity', max(intensity), np.mean(intensity)
    #print 'max of 277atm', max(black_body_277), np.mean(black_body_277)
    #Calculate PWV equivalent of emissivity
    in_band_flag = np.logical_and(frequency >= frequency_low, frequency <= frequency_high)
    frequency_band = np.extract(in_band_flag, frequency)
    #blackbody_band_atm = np.extract(in_band_flag, black_body_ATM)
    blackbody_band_atm2 = np.extract(in_band_flag, black_body_ATM2)
    blackbody_band_277 = np.extract(in_band_flag, black_body_277)
    transmission_band = np.extract(in_band_flag, transmission)

    #integrated_blackbody_band_atm = np.trapz(blackbody_band_atm, frequency_band)
    integrated_blackbody_band_277 = np.trapz(blackbody_band_277, frequency_band)
    integrated_blackbody_band_atm2 = np.trapz(blackbody_band_atm2, frequency_band)  ###
    #print 'atm blackbody, integrated', integrated_blackbody_band_atm
    #print 'atm2 blackbody, integrated', integrated_blackbody_band_atm2
    #print '277 blackbody, integrated', integrated_blackbody_band_277

    effective_emissivity = integrated_blackbody_band_atm2 / integrated_blackbody_band_277
    mean_transmission = np.mean(transmission_band)
    return (effective_emissivity, mean_transmission)




if __name__ == '__main__':
    raw_file_name_dict = {'Chajnantor_at_60deg_1p0mmH20.out': 'r',
                          'SouthPole_60Deg1p0mm.out': 'b',
                          'LDB_at_34p61km.out': 'k'}
    band_vs_transmission_plotter()
    #transmission_plotter(raw_file_name_dict)

