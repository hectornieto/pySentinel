# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:30:19 2018

@author: hnieto
"""

import subprocess as sp
import os.path as pth
import glob

def build_gdal_vrt(input_files, output_file):
    
    command='gdalbuildvrt -separate "%s.vrt" '%output_file
    
    for file in input_files:
        command+='"%s" '%file
    
    proc=sp.Popen(command,shell=True,stdout=sp.PIPE,stdin=sp.PIPE,
                  stderr=sp.STDOUT,universal_newlines=True)
    
    proc.wait()

if __name__=='__main__':

    valid_bands=('B2', 'B3','B4','B5','B6','B7','B8A','B11','B12',
                 'view_azimuth_mean','view_zenith_mean')
     
    input_folder=pth.join('G:/','Sentinel','Sentinel-2','Borden')
    
    images=glob.glob(pth.join(input_folder,'S2?_MSIL2A_*.data'))
    for image in images:
        filename=pth.basename(image).strip('.data')
        output_file=pth.join(input_folder,filename)
        
        input_files=[]
        for band in valid_bands:
            input_files.append(pth.join(image,'%s.img'%band))
        
        build_gdal_vrt(input_files, output_file)
            
    
    


        
        