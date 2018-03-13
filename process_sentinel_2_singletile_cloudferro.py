# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:27:01 2017

@author: hnieto
"""
import os
import os.path as pth
import pySentinel.sentinel_process as sen
import glob
import sys
import shutil
import gdal
import multiprocessing
import re
import numpy as np


nproc=7 # Set the number of CPUs to allow multiprocessing
biophysical_variables=('lai','fapar','lai_cab','lai_cw','fcover') # Set the biophysical variables to keep
# We remove aerosol, water vapour and cirrus bands
reflectance_variables=('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
                       'view_zenith_mean','view_azimuth_mean','sun_zenith',
                       'sun_azimuth','quality_snow_confidence',
                       'quality_cloud_confidence',
			'quality_scene_classification',
			'quality_aot','quality_wvp')
# List of SEN-ET sites
sites={'Voulund':   (56.0376444,9.1593833),
        'Gludsted':  (56.0736389,9.3339972),
        'Hyltemossa':(56.1,     13.4166667),
        'Dahra':    (15.40278,  -15.43222000),
        'Borden':   (38.289355, -121.11779),
        'Choptank': (38.8661,   -75.8903),
        'Kairouan': (35.7055556,9.6958333),
        'Aurade':   (43.559,    1.06434000),
        'Majadas':  (39.9403315,-5.7746471),
        'Harvard':  (42.5378,   -72.1715),
        'Burdette': (35.8089,   -90.0327),
        'Walnut_Gulch':(31.7365,-109.9419),
        'Stordalen':(68.35,     19.05),
        'Degero':   (64.1833333,19.55)}




workdir=os.getcwd()
site='Majadas'# set which site to process
image_list_file=pth.join(workdir,'S2','product_list_%s'%site) # point to the list of image names to process

cloudbasedir='/eodata/Sentinel-2/MSI/' # Base folder where the sentinel-2 L1c images are located

resolution=20 # Set the output desired resoltuion
QC_VALID_VALUES=[4,5]

# Create the output dir, which will be a subfolder with the name of the site
outdir=pth.join(workdir,'S2',site)
if not pth.isdir(outdir):
    os.mkdir(outdir)

# Read the file list and extrac a list with all the files to process
file_list=np.genfromtxt(image_list_file,dtype=None, names=True)
file_list=file_list['filename']

# Create the processing dictionary options
biophysical_variables_dict={}
for var in biophysical_variables:
    biophysical_variables_dict[var]=var
reflectance_variables_dict={}
for var in reflectance_variables:
    reflectance_variables_dict[var]=var

# Create the list of variables to keep
keep_variables=list(biophysical_variables)
for var in reflectance_variables:
    keep_variables.append(var) 

print("S2 atmospheric correction...")

# Get the list of already processed files
processed_files=[]
if not pth.isfile(pth.join(outdir,'processed_%s.txt'%site)):
    logfid=open(pth.join(outdir,'processed_%s.txt'%site),'w')
    logfid.close()
else:
    logfid=open(pth.join(outdir,'processed_%s.txt'%site),'r')
    for line in logfid.readlines():
        processed_files.append(line.rstrip('\n').rstrip('\r'))
    logfid.close()
    
    
jobArgs=[]
print("S2 atmospheric correction...")
for filename in file_list:
    filename=str(filename.decode())
    test_filename=filename.rstrip('.SAFE')[11:]
    
    if test_filename not in processed_files:
            
        # Extract day, year month from the filename
        year=filename[11:15]
        month=filename[15:17]
        day=filename[17:19]
        # Cloudferro Cloud is structured in dates
        clouddir=pth.join(cloudbasedir,'L1C',year,month,day)
        
        # Setting-up multiprocessing of several Sen2Cor instances
        pool = multiprocessing.Pool(nproc)
        # Prepare sen2cor configuration, L2A will be produced at the output dir
        sen2CorGippFile = re.sub("_template", "", sen.sen2CorTemplateFile)
        sen.prepareSen2CorGippFile(sen.sen2CorTemplateFile, sen2CorGippFile, outdir)
        l1c_file=pth.join(clouddir,filename) # Full path of the input L1c file
        print(test_filename)
        if pth.exists(l1c_file) :  
            # For each tile call Sen2Cor in parallel processes
            # for all the S2 tiles in that date
            jobArgs.append((l1c_file,sen2CorGippFile,resolution))


# Atomspheric correction of all files in the list via multiprocessing    
target = sen.sen2cor
pool.map(target, jobArgs)

#Once done, find all L2A files produced
l2a_files=glob.glob(pth.join(outdir,'S2?_MSIL2A_*.SAFE'))
print("S2 SNAP processing...")

for l2a_file in l2a_files:
    filename=pth.basename(l2a_file).replace('.SAFE','_L2B_%sm'%resolution)
    test_filename=pth.basename(l2a_file).rstrip('.SAFE')
    test_filename=test_filename[11:]
    print(test_filename)
    if test_filename not in processed_files:
        l2b_file=sen.sen2lai(l2a_file, resolution=resolution, calcLAI=True, calcFAPAR=True, calcCab=True, CalcCw=True, CalcFVC=True, outdir=outdir,paralellism=120)
        if l2b_file:
            sen.delete_DIMAP_bands(l2b_file,keep_variables)
        	
        QC_image=pth.join(l2b_file+'.data','quality_scene_classification.img')
        fid=gdal.Open(QC_image,gdal.GA_ReadOnly)
        qc_array=fid.GetRasterBand(1).ReadAsArray()
        mask=sen.create_binary_mask(qc_array, valid_values=QC_VALID_VALUES)
        maskfile=pth.join(l2b_file+'.data','mask.tif')
        sen.saveImg (mask, fid.GetGeoTransform(), fid.GetProjection(), maskfile, dtype=gdal.GDT_Byte)
        
        logfid=open(pth.join(outdir,'processed_%s.txt'%site),'a')
        logfid.write('\n'+test_filename)
        logfid.flush()
        logfid.close()
        
    else:
        print('Image %s already processed'%filename)
