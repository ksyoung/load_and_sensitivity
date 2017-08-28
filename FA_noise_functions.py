#! /usr/bin/env python
'''
noise functions from Francios's code in leap, 
     leap/apps/bolometer_noise/noise_spectrum_analysis/noise_prediction/noise_prediction.py

These are just straight copies as of August 25 2017 so I can call parts of his
code in my noise calculations.

Karl Young, UMN, 8-25-2017

'''

import pdb
import copy
import numpy as np
import pylab as pl
from operator import itemgetter
from leap.lib.leap_app import leap_app
import pywtl.core.wtl_ConvertUtils as conv
import pywtl.common.analysis.noise.analysis.NoisePred as NP
import pywtl.common.analysis.noise.analysis.DataHandlingLib as DataHand
from pywtl.common.PyPolarbear.Mapping import DeviceManager
import pywtl.common.analysis.noise.analysis.ParameterLib as par
#from leap.lib.io_management import dirfile_loading

def extract_noise_from_text(settings, text):
        """
        Given the noise text detailed result, it is split up in some types of noise
        """
        text = text.split('\n')
        text.pop()

        noise = []
        for i in range(len(text)):
            noise.append(text[i].split(': ')[1])
            noise[-1] = float(noise[-1].split(' ')[-2])

        noise_values = {}
        # the order of noise sources is dfferent for SQUID noise since there are less sources
        if settings.noise_type == "darksquid":
            noise_values['cold'] = pl.sqrt(sum(pl.array([noise[0]])**2))
            noise_values['warm'] = pl.sqrt(sum(pl.array(itemgetter(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)(noise))**2))
        else:
            noise_values['cold'] = pl.sqrt(sum(pl.array(itemgetter(0, 10, 11)(noise))**2))
            noise_values['warm'] = pl.sqrt(sum(pl.array(itemgetter(1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17,
                                                                   18)(noise))**2))
        if len(noise) >= 20:
            noise_values['johnson'] = noise[19]
        else:
            noise_values['johnson'] = pl.nan
        if len(noise) >= 21:
            noise_values['phonon'] = noise[20]
        else:
            noise_values['phonon'] = pl.nan
        if len(noise) >= 22:
            noise_values['photon'] = noise[21]
        else:
            noise_values['photon'] = pl.nan
        if len(noise) >= 23:
            noise_values['bunch'] = noise[22]
        else:
            noise_values['bunch'] = pl.nan
        noise_values['total'] = pl.sqrt(sum(pl.array(noise)**2))

        return noise_values

def summarize_noise(settings,noise):
    """
    Given a dictionary with the proper noise types, a summary txt is written
    """
    summary = [pl.nan, pl.nan, pl.nan, pl.nan]
    if settings.johnson_and_readout_factor != None:
      text = "{0:20s} : {1:.2f}\n".format('Johnson and ro boost', settings.johnson_and_readout_factor)
      text += "\n"
    else:
      text = ""

    # current at SQUID
    value = pl.sqrt(noise['warm']**2 + noise['cold']**2)
    text += "{0:20s} : {1:7.3f} pA/sqrt(Hz)\n".format('readout', value)
    text += "{0:20s} : {1:7.3f} pA/sqrt(Hz)\n".format('Johnson', noise['johnson'])
    text += "{0:20s} : {1:7.3f} pA/sqrt(Hz)\n".format('phonon', noise['phonon'])
    value = pl.sqrt(noise['photon']**2 + noise['bunch']**2)
    text += "{0:20s} : {1:7.3f} pA/sqrt(Hz)\n".format('photon', value)
    text += "{0:20s} : {1:7.3f} pA/sqrt(Hz)\n".format('total', noise['total'])
    text += '\n'
    summary[0] = copy.deepcopy(noise['total']*1e-12)

    # readout power at bolo
    value = pl.sqrt(noise['warm']**2 + noise['cold']**2) * 1e6 / NP.Si(settings.dfmux_settings['Vb'],
                                                                       False, settings.bolo_char['L'])
    text += "{0:20s} : {1:7.3f} aW/sqrt(Hz)\n".format('readout', value)
    value = noise['johnson'] * 1e6 / NP.Si(settings.dfmux_settings['Vb'], False, settings.bolo_char['L'])
    text += "{0:20s} : {1:7.3f} aW/sqrt(Hz)\n".format('Johnson', value)
    value = noise['phonon'] * 1e6 / NP.Si(settings.dfmux_settings['Vb'], False, settings.bolo_char['L'])
    text += "{0:20s} : {1:7.3f} aW/sqrt(Hz)\n".format('phonon', value)
    value = pl.sqrt(noise['photon']**2 + noise['bunch']**2) * 1e6 / NP.Si(settings.dfmux_settings['Vb'],
                                                                          False, settings.bolo_char['L'])
    text += "{0:20s} : {1:7.3f} aW/sqrt(Hz)\n".format('photon', value)
    value = noise['total'] * 1e6 / NP.Si(settings.dfmux_settings['Vb'], False, settings.bolo_char['L'])
    text += "{0:20s} : {1:7.3f} aW/sqrt(Hz)\n".format('total', value)
    text += '\n'
    summary[1] = copy.deepcopy(value*1e-18)

    # counts
    if settings.A_per_count is not None:
        value = pl.sqrt(noise['warm']**2 + noise['cold']**2) * 1e-12 / settings.A_per_count
        text += "{0:20s} : {1:7.3f} counts/sqrt(Hz)\n".format('readout', value)
        value = noise['johnson'] * 1e-12 / settings.A_per_count
        text += "{0:20s} : {1:7.3f} counts/sqrt(Hz)\n".format('Johnson', value)
        value = noise['phonon'] * 1e-12 / settings.A_per_count
        text += "{0:20s} : {1:7.3f} counts/sqrt(Hz)\n".format('phonon', value)
        value = pl.sqrt(noise['photon']**2 + noise['bunch']**2) * 1e-12 / settings.A_per_count
        text += "{0:20s} : {1:7.3f} counts/sqrt(Hz)\n".format('photon', value)
        value = noise['total'] * 1e-12 / settings.A_per_count
        text += "{0:20s} : {1:7.3f} counts/sqrt(Hz)\n".format('total', value)
        text += '\n'
        summary[2] = copy.deepcopy(value)

    # Kcmb
    if settings.A_per_count is not None and settings.Kcmb_per_cnt is not None:
        fact = abs(settings.Kcmb_per_cnt) * 1e6 / settings.A_per_count
        value = pl.sqrt(noise['warm']**2 + noise['cold']**2) * 1e-12 * fact
        text += "{0:20s} : {1:7.3f} uKcmb/sqrt(Hz)\n".format('readout', value)
        value = noise['johnson'] * 1e-12 * fact
        text += "{0:20s} : {1:7.3f} uKcmb/sqrt(Hz)\n".format('Johnson', value)
        value = noise['phonon'] * 1e-12 * fact
        text += "{0:20s} : {1:7.3f} uKcmb/sqrt(Hz)\n".format('phonon', value)
        value = pl.sqrt(noise['photon']**2 + noise['bunch']**2) * 1e-12 * fact
        text += "{0:20s} : {1:7.3f} uKcmb/sqrt(Hz)\n".format('photon', value)
        value = noise['total'] * 1e-12 * fact
        text += "{0:20s} : {1:7.3f} uKcmb/sqrt(Hz)\n".format('total', value)
        text += '\n'
        summary[3] = copy.deepcopy(value*1e-6)

    return text, summary
