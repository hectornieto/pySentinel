# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:27:01 2017

@author: hnieto
"""
import os
import os.path as pth
import pySentinel.sentinel_process as sen
import pySentinel.gdal_utils as gu
import glob
import sys
import shutil
import gdal
import multiprocessing
import re

nproc=8
biophysical_variables=('lai','fapar','Cab','Cw')
# We remove aerosol, water vapour and cirrus bands
reflectance_variables=('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
                       'view_zenith_mean','view_azimuth_mean','sun_zenith',
                       'sun_azimuth','quality_snow_confidence',
                       'quality_cloud_confidence','quality_scene_classification')
workdir='/home/eouser/Data/'
cloudbasedir='/eodata/Sentinel-2/MSI/'
# Process Lleida area   
resolution=20
latlon=(0.32,42.40,1.70,41.21)
inputEPSG=4326
outputEPSG=32631
QC_VALID_VALUES=[4,5]
valid_tiles=['T31TBG','T31TBF','T31TCG','T31TCF','T31TCH']
ul_X,ul_Y,Z_out=gu.coordinate_convert((latlon[0],latlon[1]),
                                       inputEPSG,
                                       outputEPSG,
                                       Z_in=0)
                                       
lr_X,lr_Y,Z_out=gu.coordinate_convert((latlon[2],latlon[3]),
                                        inputEPSG
                                        ,outputEPSG,
                                        Z_in=0)
para_dict={'crs':'EPSG:%s'%(outputEPSG),'northBound':float(latlon[3]),'southBound':float(latlon[1]),
'eastBound':float(latlon[2]),'westBound':float(latlon[0]),'pixelSizeX':float(resolution),'pixelSizeY':float(resolution)}


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

dates=list(dates)

print("S2 atmospheric correction...")
for date in dates:
    l2b_files=[]
    year=date[0:4]
    month=date[4:6]
    day=date[6:8]
    clouddir=pth.join(cloudbasedir,'L1C',year,month,day)
    print("S2 atmospheric correction...")
    pool = multiprocessing.Pool(nproc)
    jobArgs=[]        
    sen2CorGippFile = pth.join(workdir,'Sen2Cor_L2A_GIPP.xml')
    sen.prepareSen2CorGippFile(sen2CorGippFile, workdir)
    l2a_files=[]
    for tile in valid_tiles:
       	l1c_file=glob.glob(pth.join(clouddir,'S2?_MSIL1C_%sT??????_N????_R???_%s_*.SAFE'%(date,tile)))
        print(l1c_file)
        if len(l1c_file)==1:   
            # For each tile call Sen2Cor in parallel processes
            # for all the S2 tiles in that date
            jobArgs.append((l1c_file[0],sen2CorGippFile,resolution))
            l2a_file=pth.join(workdir,re.sub('MSIL1C','MSIL2A',pth.basename(l1c_file[0])))
            l2a_files.append(l2a_file)
    
    print(jobArgs)       
    target = sen.sen2cor
    pool.map(target, jobArgs)
    l2a_files=glob.glob(pth.join(workdir,'S2?_MSIL2A_%sT??????_*.SAFE'%date))
    print("S2 SNAP processing...")
    l2b_files=[]
    if len(l2a_files)>0:
        for l2a_file in l2a_files:
            l2b_file=sen.sen2lai(l2a_file, resolution=resolution, calcLAI=True, calcFAPAR=True, calcCab=True, CalcCw=True, CalcFVC=False, outdir=workdir,paralellism=120)
            if l2b_file:
                l2b_files.append(l2b_file)
        
        for i,l2b_file in enumerate(l2b_files):
            sen.delete_DIMAP_bands(l2b_file,biophysical_variables)
            sen.delete_DIMAP_bands(l2b_file,reflectance_variables)
            l2b_files[i]=l2b_file.replace('.data','')
            
        # Mosaic tiles
        mosaicfile=pth.join(workdir,'S2_MSIL2A_%s_reflectance_%sm.bsq'%(date,resolution))  
        if pth.exists(mosaicfile):
            os.remove(mosaicfile)
        sen.mosaic_sentinel_SAGA(l2b_files,mosaicfile, reflectance_variables, resolution, extent=(ul_X,ul_Y,lr_X,lr_Y))
            
        # Create binary mask
        # Open quality_scene_classification band
        band=reflectance_variables.index('quality_scene_classification')
        fid=gdal.Open(mosaicfile,gdal.GA_ReadOnly)
        qc_array=fid.GetRasterBand(band+1).ReadAsArray()
        mask=sen.create_binary_mask(qc_array, valid_values=QC_VALID_VALUES)
        maskfile=pth.join(workdir,'Sentinel-2','L2','S2_MSIL2A_%s_mask_%sm.bsq'%(date,resolution))
        gu.save_img (mask, 
                     fid.GetGeoTransform(), 
                     fid.GetProjection(), 
                     maskfile, 
                     dtype=gdal.GDT_Byte)
    
        # Mosaic tiles
        mosaicfile=pth.join(workdir,'Sentinel-2','L2','S2_MSIL2A_%s_biophysical_%sm.bsq'%(date,resolution))  
        if pth.exists(mosaicfile):
            os.remove(mosaicfile)
        sen.mosaic_sentinel_SAGA(l2b_files,mosaicfile, biophysical_variables, resolution, extent=(ul_X,ul_Y,lr_X,lr_Y))
        for file in l2b_files:
            try:
                os.remove(file+'.dim')
                shutil.rmtree(file+'.data')
            except:
                print('%s not found, delete skipped'%file)
            
        shutil.rmtree(pth.join(workdir,'Sentinel-2','L2','tmp'))          
