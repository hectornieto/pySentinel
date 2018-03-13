# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 13:32:01 2017

@author: radoslaw guzinski
"""

import os
import numpy as np
from osgeo import gdal, gdalconst
from pyDMS import pyDMSUtils as u

def get_landcover_data(landcoverSource, LAI, NDVI, EVI, dataFile, landcoverResamplingMode = gdalconst.GRA_NearestNeighbour, cropSensence = False):
    
    # land cover specific data
    # Crop land growth limits are fitted to barley based on
    # https://cereals.ahdb.org.uk/media/186381/g67-barley-growth-guide.pdf
    
    # IGBP classification
    #                       Everg   Everg.  Decid.  Decid.                                                                          Crop
    #                       Needle  Broad   Needle  Broad   Mixed   Closed  Open    Woody           Grass   Wet     Crop            Veg     Snow
    #                 Water Forest  Forest  Forest  Forest  Forest  shrubs  shrubs  Savanna Savanna land    lands   lands  Urban    Mosaic  Ice     Baren   
    landCoverClasses= [0,   1,      2,      3,      4,      5,      6,      7,      8,      9,      10,     11,     12,    13,      14,     15,     16 ];  # land cover class, from MOD12Q1 
    C_h_C =           [0.0, 5.0,    5.0,    5.0,    5.0,    5.0,    1.5,    1.5,    10.0,   8.0,    0.5,    1.0,    1.2,   5.0,     1.2,   1.0,    1.0];  # max height
    C_defaultLAI =    [0.01,0.01,   0.01,   0.01,   0.01,   0.01,   0.01,   4.0,   0.01,   0.01,    4.0,    0.01,   5.0,   0.01,    5.0,   0.01,   1.0];  # fully developed LAI
    #C_f_C =           [1.0, 0.8,    0.9,    0.8,    0.9,    0.9,    0.9,    0.9,    0.9,    1.0,    1.0,    1.0,    0.9,   1.0,     0.9,    1.0,    1.0];  # Aparent (large scale) fractional vegetation cover
    C_w_C =           [0.0, 2.0,    1.0,    2.0,    1.0,    1.5,    1.5,    1.2,    1.5,    1.0,    1.0,    1.0,    1.0,   0.0,     1.0,    0.0,    1.0];  # ratio of canopy height to width, forest values taken from K.J. Schaudt, R.E. Dickinson 2000 (table 1)
    C_s =             [0.0, 0.05,   0.10,   0.05,   0.10,   0.07,   0.10,   0.05,    0.02,   0.1,   0.02,   0.02,   0.2,   0.0,     0.2,    0.0,    0.02]; # leaf size, based on Houbourg 2009 paper
    #C_z_T =           [2.0, 24.0,   19.0,   24.0,   19.0,   24.0,   4.0,    4.0,    15.0,   50.0,   4.0,   50.0,   4.0,  50.0,    4.0,   50.0,   50.0]; # modelled temp 50m above ground/canopy
    #C_z_U =           [10.0,26.0,   21.0,   26.0,   21.0,   26.0,   6.0,   6.0,   20.0,   13.0,    6.0,   10.0,   6.0,  15.0,    6.0,   10.0,   10.0]; # modelled wind 10m above ground/canopy
    C_x_LAD =         [0.0, 1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    0.2,    0.5,    0.8,    0.5,    0.5,    0.5,   0.0,     0.5,    0.0,    0.0];  # Cambpbell 1990 leaf inclination distribution parameter:[x_LAD=1 for spherical LIDF, x_LAD=0 for vertical LIDF, x_LAD=float(inf) for horzontal LIDF]
    C_tau_NIR =       [0.0, 0.3,    0.4,    0.3,    0.4,    0.35,   0.4,    0.4,    0.4,    0.4,    0.4,    0.4,    0.4,   0.0,     0.4,    0.0,    0.0] # leaf transmittance in NIR. Confier taken from Mesarch et al. 1999, other from standard leaf spectra.
    C_tau_VIS =       [0.0, 0.05,   0.15,   0.05,   0.15,  0.10,   0.15,   0.15,   0.15,   0.15,   0.15,   0.15,   0.15,  0.0,     0.15,   0.0,    0.0] # leaf transmittance in VIS. Confier taken from Mesarch et al. 1999, other from standard leaf spectra.

    data = {}

    
    # Read the landcover data from file if it looks like landcoverSource
    # is a raster 
    if os.path.exists(landcoverSource) or landcoverSource.startswith("HDF5:") \
      or landcoverSource.startswith("HDF4_EOS:") \
      or landcoverSource.startswith("NETCDF:"):
        im = gdal.Open(landcoverSource, gdal.GA_ReadOnly)
        #landcoverData = im.GetRasterBand(1).ReadAsArray()
    else:
        raise AttributeError
        
    # Inputs created by SNAP do not always have the same extent. Therfore, 
    # resample and subset the land cover raster to the correct extent.
    #resampled = u.resampleWithGdal(landcoverData, im.GetGeoTransform(), 
    #                               dataFile.GetGeoTransform(), im.GetProjection(), 
    #                               dataFile.GetProjection(), LAI.shape, 
    #                               resampling = landcoverResamplingMode)
    resampled = u.resampleWithGdalWarp(landcoverSource, dataFile, resampleAlg = landcoverResamplingMode)
    landcoverData = resampled.GetRasterBand(1).ReadAsArray()
    im = None
    resampled = None
    
    # Create arrays for the landcover dependent parameters
    data['F_G'] = np.zeros(landcoverData.shape)
    data['H_C'] = np.zeros(landcoverData.shape)
    data['W_C'] = np.zeros(landcoverData.shape)
    data['S'] = np.zeros(landcoverData.shape)
    data['Z_T'] = np.zeros(landcoverData.shape)
    data['Z_U'] = np.zeros(landcoverData.shape)
    data['X_LAD'] = np.zeros(landcoverData.shape)
    data['TAU_NIR'] = np.zeros(landcoverData.shape)
    data['TAU_VIS'] = np.zeros(landcoverData.shape)
        
    # Set the parameters for each present land cover class 
    for lcClass in np.unique(landcoverData[~np.isnan(landcoverData)]):
        lcPixels = np.where(landcoverData == lcClass)
        lcIndex = landCoverClasses.index(lcClass)
        
        data['W_C'][lcPixels] = C_w_C[lcIndex]
        data['S'][lcPixels] = C_s[lcIndex]
        data['Z_T'][lcPixels] = 100.0 #30.0 # C_z_T[lcIndex] #  
        data['Z_U'][lcPixels] = 10.0 + C_h_C[lcIndex] #38.0 #C_z_U[lcIndex] #  
        # In croplands use EVI as proxy of green vegetation fraction only during
        # senescence. In other land covers use the equation from Fisher et al. 2008
        if lcClass == 12 or lcClass == 14:
            if cropSensence:
                data['F_G'][lcPixels] = EVI[lcPixels]
            else:
                data['F_G'][lcPixels] = 1.0
        else:
            data['F_G'][lcPixels] = np.minimum(1.2 * EVI[lcPixels] / NDVI[lcPixels], 1.0)
        # Scale height with PAI for grasslands, croplands, etc.
        if lcClass == 10 or lcClass == 12 or lcClass == 14:
            PAI = LAI[lcPixels] / data['F_G'][lcPixels] 
            data['H_C'][lcPixels] = 0.1*C_h_C[lcIndex] + 0.9*C_h_C[lcIndex]*np.minimum((PAI/C_defaultLAI[lcIndex])**3.0, 1.0)
        else:
            data['H_C'][lcPixels] = C_h_C[lcIndex]
        
        data['X_LAD'][lcPixels] = C_x_LAD[lcIndex]
        data['TAU_NIR'][lcPixels] = C_tau_NIR[lcIndex]
        data['TAU_VIS'][lcPixels] = C_tau_VIS[lcIndex]
        
    data['lc_IGBP'] = landcoverData
        
    return data
  