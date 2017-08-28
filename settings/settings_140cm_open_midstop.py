# -*- coding: utf-8 -*-
import os
from pylab import sqrt
from pywtl.core.wtl_ConvertUtils import convert_squid
import pywtl.common.analysis.noise.analysis.NoisePred as NP
import pywtl.common.analysis.noise.analysis.ParameterLib as PL
import numpy as np

class Class(object):
    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__, str(self.__dict__))

settings = Class()
# Run time paramters
settings.freq = 'All_GHz'
settings.version = 'open_midstop'
settings.name = '1.4m_open_midstop'
settings.verbose = True

# Telescope/ Receiver Optical Parameters
settings.mult_bands = True
#settings.band = [133.,167.]  # GHz, lower, upper.  150, width = 34.
#settings.band = np.array(settings.band)*1e9  # Hz

settings.aperture_radius = 0.7  # aperture radius in meters (2.5 meter primary = 1.25m radius)
settings.f_number = 1.5  # 

settings.edge_db = 10   # edge taper on primary mirror in dB.  May later be calculated from pixel sizes.

# Bolo parameters
settings.t_bath = 0.25  # Kelvin
settings.safety_factor = 2.5  # Unitless, ratio of p_sat to p_opt
settings.n = 3.0            # thermal power law exponent (EBEX was ~2)
settings.bolo_Rn = 1.33  # Ohms.  TES resistance warm.
settings.bias_point = 0.75  # depth in transition assumed
settings.bolo_resistance = settings.bolo_Rn*settings.bias_point  # Ohms
## old readout noise method.
settings.readout_noise_amps = 7e-12 # Amps*rt(sec), number from Franky for scaling readout noise.

# More bolo parameters for franky's noise code.
settings.conv = convert_squid('Normalized_16ch')

settings.prediction_type = 'theoretical'

dfmux_settings = {}
squid_settings = {}
bolo_char = {}

settings.noise_type = "transition"

# boost factors for noise.
settings.johnson_and_readout_factor = None  # use un-modified noise theory.

# DfMUX general setup
dfmux_settings['DAN_firmware'] = True
dfmux_settings['DAN_parser'] = False
dfmux_settings['bitshift'] = 5  # can be 8 for 24 bit.
dfmux_settings['fir'] = 6
dfmux_settings['fsamp'] = 25e6/2**(11+dfmux_settings['fir'])

# DAC/ADC settings
dfmux_settings['fc'] = 500000.
dfmux_settings['Gc'] = 2
dfmux_settings['fn'] = dfmux_settings['fc']
dfmux_settings['Gn'] = 1
dfmux_settings['fd'] = dfmux_settings['fc']
dfmux_settings['Gd'] = 0

# SQUID/cryostat settings
squid_settings['R_FB'] = 5000.
##calced in code   # bolo_char['nu'] = 150e9
##calced in code   # bolo_char['dnu'] = 34e9
bolo_char['Zt'] = 320.26
bolo_char['L_fll'] = PL.LoopGain(bolo_char['Zt'])
settings.R_wire = 10.  # is warm wire, squid board to squid controller.
bolo_char['Tbath'] = settings.t_bath

# bolometer characteristics
##calced in code   # bolo_char['Tc'] = 0.42625 
##calced in code   # bolo_char['Tbolo'] = bolo_char['Tc']
bolo_char['Rn'] = settings.bolo_Rn
bolo_char['R'] = settings.bolo_resistance
bolo_char['tau_etf'] = 0.010
##calced in code   # bolo_char['Popt'] = .124127e-12
bolo_char['L'] = 25. ##### something reasonable, but could be any number.
bolo_char['xi'] = 1  ## assume all correlation noise.

# other bolometer characteristics
##calced in code   # bolo_char['Psat'] = 2.5 * bolo_char['Popt']
bolo_char['n'] = settings.n
settings.cryo = 'EBEX'

# derived values from Psat, Tc, Tbath and n
##calced in code  bolo_char['Gbar'] = bolo_char['Psat']/(bolo_char['Tc']-bolo_char['Tbath'])
##calced in code  bolo_char['G'] = PL.G_dyn(bolo_char['Gbar'], bolo_char['Tbath'], bolo_char['Tc'], bolo_char['n'])
##calced in code  bolo_char['gamma'] = round(PL.CalcGamma(bolo_char['Tc'], bolo_char['Tbath'], bolo_char['n']), 3)

# derived DfMUX settings
##calced in code  dfmux_settings['Vb'] = bolo_char['R'] * (bolo_char['Gbar']*(bolo_char['Tbolo']-bolo_char['Tbath']) - bolo_char['Popt'])
##calced in code  dfmux_settings['Vb'] = sqrt(dfmux_settings['Vb'])
##calced in code  dfmux_settings['Ac'] = dfmux_settings['Vb'] / conv.DDStoVbias(Carrier_amplitude=1,
#                                                              Carrier_gain=dfmux_settings['Gc'],
#                                                              firmware_version='16ch')

##calced in code  R_gain = [2000., 820., 200., 0.]
##calced in code  dfmux_settings['An'] = dfmux_settings['Ac']/3. * \
#    (R_gain[dfmux_settings['Gn']] + 100.) / (R_gain[dfmux_settings['Gc']] + 100.)

##calced in code  bolo_char['Si'] = NP.Si(dfmux_settings['Vb'])  # assumes deep in transition

# set a value to have noise in counts or Kcmb (from franky's code)
settings.A_per_count = None
settings.Kcmb_per_cnt = None


# copying those into settings
settings.dfmux_settings = dfmux_settings
settings.squid_settings = squid_settings
settings.bolo_char = bolo_char

# Paths
settings.base_path = '/home/astro/kyoung/Documents/load_and_sensitivity/'
settings.elements_path = os.path.join(settings.base_path,
                 'inputs/140cm_open_dragone_midstop.csv')  # all telescope surfaces, lenses, etc.

# now being defined in code.
#settings.elements_out_path = os.path.join(settings.base_path,
#                 'outputs/%s_%s_elements_out.csv ' %(settings.freq, settings.version))  # data that gets saved.

settings.bands_path = os.path.join(settings.base_path,
                 'inputs/CMBP_bands.csv')  # csv of bands.


# unneeded stuff below this line.
'''
# Run time paramters
settings.name = '150'
settings.version = 'greg_f3'
settings.design = '1m_EBEX10K_f3_2016'
settings.do_point_source_analysis = False
settings.footer = True
# Telescope/ Receiver Optical Parameters
settings.frequency = 150.0  # Ghz
settings.frequency_range = (133, 167.0, 150.)
settings.bandwidth = settings.frequency_range[1] - settings.frequency_range[0]  # Ghz
# Optics Parameters
settings.diameter_to_waist_ratio = 2.95  #
# Pixel Parameters
settings.pixel_diameter_step = 0.1
settings.pixel_diameter_range = [settings.pixel_diameter_step, 10 + settings.pixel_diameter_step]  # in mm 0 to 20 mm
settings.default_pixel_diameter_m = 0.0042 # In meters 
settings.default_pixel_diameter_mm = settings.default_pixel_diameter_m * 1000 # In millimeters 
settings.lens_loss_tangent = 9e-5
settings.lens_index = 3.2
settings.lens_thickness = 0.050  # In meters
#Fmux parameters, Readout Noise Contribution
settings.readout_contribution = 0.10  # readout noise should increasee noise by 10%
settings.tes_accuracy = 1.15  # TES resistance value across wafer +/- 15%
settings.do_fmux = True  # 1 to calculate fMUx params such as L value, cross talk etc
settings.fmux_max_freq = 1.0e6  # maximum frequency for fMux readout
settings.fmux_max_mux = 36  # mux factor
settings.capacitor_tan_d = 0.0002  # expected tand for interdigitated capacitor
settings.cross_talk_level = 0.01  # allowed cross-talk level
settings.capacitor_accuracy = 0.005  # fractional accuracy of capacitor
settings.f_delta_factor = 1.25  # factor to increase frequency spacing look at crosstalk_off to decide
settings.esr_contribution = 0.10  # allowed fractional ESR contribution to total R (10% = 0.01)
settings.readout_noise = 7e-12  # readout noise in A*sqrt(s)
settings.v_bias  = 1e-6  #4 microVolt voltage bias
settings.nep_readout_fixed = 9.3e-18  # in Watts/root(Hz) this number from Franky July 7th, noise email to Shaul. 
# Computatonal Parameters
settings.spatial_integration_accuracy = 10000  # numerical integration accuracy. ex: 100 splits integration space in 100 rectangles.
settings.frequency_bin_size = 0.384  # GHz
settings.integration_time = 0.5  # seconds
#Bolometer Parameters
settings.num_bolo_legs = 4
settings.bath_temp = 0.275  # Kelvin
settings.bolo_R_normal = 1.333333333  # Ohms.  TES resistance warm.
settings.bias_fraction_rnormal = 0.75
settings.bolo_resistance = settings.bolo_R_normal*settings.bias_fraction_rnormal  # Ohms
settings.thermal_carrier_exponent = 3.0
settings.a_over_l_bolo_leg = 149.2*1e-3*1e-12 # From Toki's script
settings.alpha = 250  # Measure for AlTi R v T curve d(logR)/dT
settings.tau_safety_factor = 5.8  # How much slower than readout bolo should be 
#Design parameters/Goals
settings.psat_to_optical_power = 2.5  # Unitless
settings.bias_fraction_rnormal = 0.6
settings.target_num_pixels =  924
settings.num_sub_arrays = 6
settings.max_focal_plane_diameter_x = 0.25  # in m
settings.max_focal_plane_diameter_y = 0.25  # in m
settings.max_num_pixels = 10000
#settings.max_focal_plane_diameter = FOV_deg*(np.pi/180)*2.0*settings.aperture_radius*settings.f_number
settings.hex_outer_radius = 0.096  # in m
#Observing Parameters
settings.fract_sky = 0.5
settings.obs_time = 10  #flight time in days.
settings.obs_efficiency =.8  #percent of flight time when the telescope is acutally observing.
# Input Paths
settings.base_path = '/home/astro/kyoung/Documents/35cm-xdragone/Mapping_speed_code_py'
settings.atmospheric_windows_data_file = os.path.join(settings.base_path, 'input_data_products', 'LDB_34km_30deg_el.out')  # Atmospheric Windows 
settings.receiver_throughput_data_path = os.path.join(settings.base_path, 'input_data_products',
                                                      '1m_EBEX10K_2016_Throughput_%sGHz.csv' % str(int(settings.frequency)))  # Aperture Illumination

#Aperture Output paths
settings.effective_aperture_data_path = os.path.join(settings.base_path,
                                        'output/effective_aperture',
                                        '1m_EBEX10K_f3_2016_Effective_Aperture_Output_%s_GHz_%s.dat' %(str(int(settings.frequency)),
                                        settings.version))  # Aperture Illumination

settings.effective_aperture_png_path = os.path.join(settings.base_path,
                                        'output/effective_aperture',
                                        '1m_EBEX10K_f3_2016_Effective_Aperture_Output_%s_GHz_%s.png' %(str(int(settings.frequency)),
                                        settings.version))  # Aperture Illumination

#Output Paths
settings.base_path = '/home/astro/kyoung/Documents/35cm-xdragone/Mapping_speed_code_py/output/mapping_speed/EBEX10K_f3_2016_gregorian'
settings.mapping_speed_output_png_path = os.path.join(settings.base_path, '1m_EBEX10K_f3_2016_mapping_speed_%sGHz_%s' % (str(int(settings.frequency)), settings.version))
settings.mapping_speed_output_data_path = os.path.join(settings.base_path, 'mapping_speed_%sGHz_%s.csv' % (str(int(settings.frequency)), settings.version))
settings.num_pixels_output_png_path = os.path.join(settings.base_path, 'num_pixels_%sGHz_%s.png' % (str(int(settings.frequency)), settings.version))
settings.throughput_path_out = os.path.join(settings.base_path, 'throughput_%sGHz_%s.csv' % (str(int(settings.frequency)), settings.version))

'''
