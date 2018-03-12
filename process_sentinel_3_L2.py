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

workdir=pth.join(os.getcwd(),'Sentinel-3','SLSTR')
lleida_latlon=(0.32,42.40, 1.70, 41.21)
resolution=[1000,1000]
inputEPSG=4326
outputEPSG=32631
band_list=('LST','TCWV','NDVI','biome','fraction')

ul_X,ul_Y,Z_out=sen.coordinate_convert(lleida_latlon[0],lleida_latlon[1],inputEPSG,outputEPSG,Z_in=0)
lr_X,lr_Y,Z_out=sen.coordinate_convert(lleida_latlon[2],lleida_latlon[3],inputEPSG,outputEPSG,Z_in=0)
extent=ul_X,ul_Y,lr_X,lr_Y


commands=sys.argv

if len(commands)<2:
    dates=[input('Type the date you want to process\nYYYYMMDD:')]
else:
    dates=commands[1:]
for date in dates:
    file_list=glob.glob(pth.join(workdir,'S3A_SL_2_LST____%s*.zip'%(date)))
    for input_file in file_list:
        input_file=input_file[:-4]
        reproject_file=input_file+'_reproject'
        if not pth.isdir(reproject_file+'.data'):
            reproject_file=sen.reproject_S3_SNAP(input_file,outputEPSG,resolution,extent=extent,output_file=None)
            if reproject_file:
                sen.delete_DIMAP_bands(reproject_file,band_list)
                shutil.rmtree(input_file)
        
