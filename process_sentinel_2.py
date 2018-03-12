# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:27:01 2017

@author: hnieto
"""
import os
import os.path as pth
import pySentinel.sentinel_process as sen
#import sentinel_process as sen
import glob
import sys
import shutil
import gdal

biophysical_variables=('lai','fapar')
# We remove aerosol, water vapour and cirrus bands
reflectance_variables=('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
                       'view_zenith_mean','view_azimuth_mean','sun_zenith',
                       'sun_azimuth','quality_snow_confidence',
                       'quality_cloud_confidence','quality_scene_classification')
workdir=os.getcwd()
# Process Lleida area   
resolution=30
lleida_latlon=(0.32,42.40,1.70,41.21)
inputEPSG=4326
outputEPSG=32631
QC_VALID_VALUES=[4,5]

ul_X,ul_Y,Z_out=sen.coordinate_convert(lleida_latlon[0],lleida_latlon[1],inputEPSG,outputEPSG,Z_in=0)
lr_X,lr_Y,Z_out=sen.coordinate_convert(lleida_latlon[2],lleida_latlon[3],inputEPSG,outputEPSG,Z_in=0)
para_dict={'crs':'EPSG:%s'%(outputEPSG),'northBound':float(lleida_latlon[3]),'southBound':float(lleida_latlon[1]),
'eastBound':float(lleida_latlon[2]),'westBound':float(lleida_latlon[0]),'pixelSizeX':float(resolution),'pixelSizeY':float(resolution)}


biophysical_variables_dict={}
for var in biophysical_variables:
    biophysical_variables_dict[var]=var
reflectance_variables_dict={}
for var in reflectance_variables:
    reflectance_variables_dict[var]=var
                              
commands=sys.argv
if len(commands)<2:
    dates=[input('Type the date you want to process\nYYYYMMDD:')]
else:
    dates=commands[1:]
for date in dates:
    l2b_lleida_biophysical=[]
    l2b_lleida_resample=[]
    l1c_files=glob.glob(pth.join(workdir,'Sentinel-2','S2?_MSIL1C_%sT*.SAFE.zip'%(date)))
    if len(l1c_files)>0:   
        for l1c_file in l1c_files:
            l1c_file=l1c_file.replace('.zip','')
            l2b_file=l1c_file.replace('.SAFE','_biophysical_%sm'%resolution)
            l2b_file=l2b_file.replace('_MSIL1C_','_MSIL2A_')
            if not pth.isdir(l2b_file+'.data'):
                l2b_file=sen.sen2lai(l1c_file,resolution=resolution)
        l2b_lleida_biophysical=glob.glob(pth.join(workdir,'Sentinel-2','S2?_MSIL2A_%sT*_biophysical_%sm.data'%(date,resolution)))
        
        for i,l2b_file in enumerate(l2b_lleida_biophysical):
            sen.delete_DIMAP_bands(l2b_file,biophysical_variables)
            l2b_lleida_biophysical[i]=l2b_file.replace('.data','')
        
        l2b_lleida_resample=glob.glob(pth.join(workdir,'Sentinel-2','S2?_MSIL2A_%sT*_resample_%sm.data'%(date,resolution)))
        for i,l2b_file in enumerate(l2b_lleida_resample):
            sen.delete_DIMAP_bands(l2b_file,reflectance_variables)
            l2b_lleida_resample[i]=l2b_file.replace('.data','')
            
        # Mosaic tiles
        mosaicfile=pth.join(workdir,'Sentinel-2','L2','S2_MSIL2A_%s_resample_%sm.bsq'%(date,resolution))  
        if pth.exists(mosaicfile):
            os.remove(mosaicfile)
        sen.mosaic_sentinel_SAGA(l2b_lleida_resample,mosaicfile, reflectance_variables, resolution, extent=(ul_X,ul_Y,lr_X,lr_Y))
        for file in l2b_lleida_resample:
            try:
                os.remove(file+'.dim')
                shutil.rmtree(file+'.data')
            except:
                print('%s not found, delete skipped'%file)
            
        # Create binary mask
        # Open quality_scene_classification band
        band=reflectance_variables.index('quality_scene_classification')
        fid=gdal.Open(mosaicfile,gdal.GA_ReadOnly)
        qc_array=fid.GetRasterBand(band+1).ReadAsArray()
        mask=sen.create_binary_mask(qc_array, valid_values=QC_VALID_VALUES)
        maskfile=pth.join(workdir,'Sentinel-2','L2','S2_MSIL2A_%s_mask_%sm.bsq'%(date,resolution))
        sen.saveImg (mask, fid.GetGeoTransform(), fid.GetProjection(), maskfile, dtype=gdal.GDT_Byte)

        # Mosaic tiles
        mosaicfile=pth.join(workdir,'Sentinel-2','L2','S2_MSIL2A_%s_biophysical_%sm.bsq'%(date,resolution))  
        if pth.exists(mosaicfile):
            os.remove(mosaicfile)
        sen.mosaic_sentinel_SAGA(l2b_lleida_biophysical,mosaicfile, biophysical_variables, resolution, extent=(ul_X,ul_Y,lr_X,lr_Y))
        for file in l2b_lleida_biophysical:
            try:
                os.remove(file+'.dim')
                shutil.rmtree(file+'.data')
            except:
                print('%s not found, delete skipped'%file)

            
        shutil.rmtree(pth.join(workdir,'Sentinel-2','L2','tmp'))
            
    else:
        print('No files found for date  %s'%(date))
