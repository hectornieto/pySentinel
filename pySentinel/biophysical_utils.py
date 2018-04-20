# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 08:56:35 2018

@author: hnieto
"""

import numpy as np
from pySentinel import gdal_utils as gu
from pyTSEB import TSEB as tseb

def fg_from_fapar(fapar, LAI, sza, f_c = 1, x_LAD = 1):
    
    # Sentinel 2 BioPhysical product estimates FPAR assuing an homogeneous canopy
    # Omega0 = 1 for coherence between fapar and fipar
    fipar = tseb.calc_F_theta_campbell(sza, 
                                         LAI, 
                                         w_C  = 1, 
                                         Omega0 = 1, 
                                         x_LAD = x_LAD)
    
    f_g = fapar / fipar
    f_g = np.clip(f_g, 0., 1.)
    return f_g

def leaf_spectra_from_biophysical(cab, cw):
    
    cab = np.clip(cab, 0., 140.)
    cw = np.clip(cw, 0., 0.1)
    
    rho_leaf_vis, tau_leaf_vis = cab_to_vis_spectrum(cab)
    rho_leaf_nir, tau_leaf_nir = cw_to_nir_spectrum(cw)
    
    return rho_leaf_vis, tau_leaf_vis, rho_leaf_nir, tau_leaf_nir
    
def cab_to_vis_spectrum(cab, 
                        coeffs_wc_rho_vis = [0.14096573, -0.09648072, -0.06328343],
                        coeffs_wc_tau_vis = [0.08543707, -0.08072709, -0.06562554]):
    
    
    
    rho_leaf_vis = watercloud_model(cab, *coeffs_wc_rho_vis)
    tau_leaf_vis = watercloud_model(cab, *coeffs_wc_tau_vis)

    rho_leaf_vis = np.clip(rho_leaf_vis, 0, 1)
    tau_leaf_vis = np.clip(tau_leaf_vis, 0, 1)
    
    return rho_leaf_vis, tau_leaf_vis

def cw_to_nir_spectrum(cw,                        
                       coeffs_wc_rho_nir = [0.38976106, -0.17260689, -65.7445699],
                       coeffs_wc_tau_nir = [0.36187620, -0.18374560, -65.3125878]):

    
    rho_leaf_nir = watercloud_model(cw, *coeffs_wc_rho_nir)
    tau_leaf_nir = watercloud_model(cw, *coeffs_wc_tau_nir)
    
    rho_leaf_nir = np.clip(rho_leaf_nir, 0, 1)
    tau_leaf_nir = np.clip(rho_leaf_nir, 0, 1)
    
    return rho_leaf_nir, tau_leaf_nir
    
def watercloud_model(param, a, b, c):
    
    result = a + b * (1.0 - np.exp(c * param))
    
    return result    