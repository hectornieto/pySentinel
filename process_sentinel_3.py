# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 18:03:09 2017

@author: hnieto
"""
from __future__ import print_function
import os.path as pth
import os
import shutil
import glob
import gdal
import numpy as np
import pySentinel.sentinel_process as sen
import sys

VALID_PIXEL=1
emis_bands=([0.98,0.95],[0.97,0.94])
workdir=pth.join(os.getcwd(),'Sentinel-3','SLSTR')
lleida_latlon=(0.32,42.40, 1.70, 41.21)
resolution=[1000,1000]
inputEPSG=4326
outputEPSG=32631
band_list=('S8_BT_in','S9_BT_in','S8_BT_io','S9_BT_io','cloud_in','cloud_io','sat_zenith_tn','sat_zenith_to','total_column_water_vapour_tx','temperature_tx','surface_pressure_tx')

ul_X,ul_Y,Z_out=sen.coordinate_convert(lleida_latlon[0],lleida_latlon[1],inputEPSG,outputEPSG,Z_in=0)
lr_X,lr_Y,Z_out=sen.coordinate_convert(lleida_latlon[2],lleida_latlon[3],inputEPSG,outputEPSG,Z_in=0)
extent=ul_X,ul_Y,lr_X,lr_Y


commands=sys.argv

if len(commands)<2:
    dates=[input('Type the date you want to process\nYYYYMMDD:')]
else:
    dates=commands[1:]
for date in dates:
    file_list=glob.glob(pth.join(workdir,'S3A_SL_1_RBT____%s*.zip'%(date)))
    for input_file in file_list:
        input_file=input_file[:-4]
        reproject_file=input_file+'_reproject'
        if not pth.isdir(reproject_file+'.data'):
            reproject_file=sen.reproject_S3_SNAP(input_file,outputEPSG,resolution,extent=extent,output_file=None)
            if reproject_file:
                sen.delete_DIMAP_bands(reproject_file,band_list)
                shutil.rmtree(input_file)
        
        # Read total column water vapour
        reproject_file+='.data'
        cloudFile=pth.join(reproject_file,'cloud_in.img')
        fid=gdal.Open(cloudFile,gdal.GA_ReadOnly)
        cloud = fid.GetRasterBand(1).ReadAsArray()
        fid = None
        valid=cloud<=VALID_PIXEL
        if not valid.any():
            print('Not valid data found')
            continue
            
        # Read total column water vapour
        waterVapourFile=pth.join(reproject_file,'total_column_water_vapour_tx.img')
        fid=gdal.Open(waterVapourFile,gdal.GA_ReadOnly)
        totalColumnWaterVapour = fid.GetRasterBand(1).ReadAsArray()
        # Convert from kg/m**2 to g/cm**2
        totalColumnWaterVapour[valid] = totalColumnWaterVapour[valid] * 0.1
        totalColumnWaterVapour[~valid]=np.nan
        fid = None
        
        # Read view zenith angles
        viewZenithAngleFile=pth.join(reproject_file,'sat_zenith_tn.img')
        fid=gdal.Open(viewZenithAngleFile,gdal.GA_ReadOnly)
        viewZenithAngle = fid.GetRasterBand(1).ReadAsArray()
        viewZenithAngle[~valid]=np.nan
        fid = None    
        
        # Read brightness temperatures
        BT11_File=pth.join(reproject_file,'S8_BT_in.img')
        fid=gdal.Open(BT11_File,gdal.GA_ReadOnly)
        bt11 = fid.GetRasterBand(1).ReadAsArray()
        bt11[~valid]=np.nan
        fid = None    
       
        # Read brightness temperatures
        BT12_File=pth.join(reproject_file,'S9_BT_in.img')
        fid=gdal.Open(BT12_File,gdal.GA_ReadOnly)
        bt12 = fid.GetRasterBand(1).ReadAsArray()
        bt12[~valid]=np.nan
        
        # Resample Sentinel2 LAI file
        inFile=pth.join(os.getcwd(),'Sentinel-2','L2','S2A_MSIL2A_%s_biophysical_30m.bsq'%(date))
        outFile = pth.join(workdir,'L2','S2A_MSIL2A_%s_LAI_1000m.tif'%(date))
        file_LAI=sen.resample_GDAL(inFile, fid.GetGeoTransform() , fid.GetProjection(),outFile =outFile, band_id=1, resampling = gdal.gdalconst.GRA_Average)
        LAI = file_LAI.GetRasterBand(1).ReadAsArray()
        LAI[~valid]=np.nan
        file_LAI=None
    
        alpha=57.3
        emis_veg=[0.98,0.975]
        emis_soil=[0.95,0.945]
        emissivity=np.zeros((LAI.shape[0],LAI.shape[1],len(emis_bands)))
        for i,emis_band in enumerate(emis_bands):
            emissivity[:,:,i]=sen.calc_emissivity_4SAIL(LAI,viewZenithAngle,alpha,emis_veg=emis_band[0],emis_soil=emis_band[1])
        brightnessTemperature=np.dstack((bt11,bt12))
        LST=sen.calc_LST_Sobrino(brightnessTemperature, emissivity, totalColumnWaterVapour, viewZenithAngle)
        LST[~valid]=np.nan
        output_filename=pth.basename(input_file)
        outPath=pth.join(workdir,'L2','%s_LST_n.tif'%(output_filename))
        sen.saveImg (np.dstack((LST,emissivity[:,:,0],emissivity[:,:,1])), fid.GetGeoTransform(), fid.GetProjection(), outPath, noDataValue =np.nan)
        fid = None
        #shutil.rmtree(reproject_file)
        os.remove(reproject_file.replace('.data','.dim'))