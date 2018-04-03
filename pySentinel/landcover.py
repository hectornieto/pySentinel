# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 10:29:31 2018

@author: hector
"""

import numpy as np
from pySentinel import gdal_utils as gu
import os.path as pth
import gdal


LC_FILENAME = 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif'

ESACCI_LC = (0, 10, 11, 12, 20, 30, 40, 50, 60, 61, 62, 70, 71,72, 80, 81, 82, 
             90, 100, 110, 120, 121, 122, 130, 140, 150, 151, 152, 153, 160, 
             170, 180, 190, 200, 201, 202, 210, 220)
                           
ESACCI_LUT = {'IGBP': (255, 12, 12, 12, 12, 14, 14, 2, 4, 4, 4, 1, 1, 1, 3, 3, 
                       3, 5, 8, 9, 6, 6, 6, 10, 16, 9, 9, 7, 10, 11, 11, 11, 
                       13, 16, 16, 16, 0, 15),
                       
              'H_C_MAX': (0, 1.2, 1, 2, 1.2, 1.2, 2, 10, 10, 10, 10, 20, 20, 
                          20, 20, 20, 20, 15, 8, 8, 1.5, 1.5, 1.5, 0.5, 0.05, 
                          2, 10, 1.5, 0.5, 10, 10, 1, 20, 0, 0, 0, 0, 0),
                          
              'LAI_MAX': (0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
                          5, 5, 5, 4, 4, 4, 4, 1, 2, 5, 4, 4, 5, 5, 5, 0, 0, 
                          0, 0, 0, 0),
                          
              'F_C': (0, 1, 1, 0.5, 1, 0.5, 0.5, 1, 1, 1, 0.4, 1, 1, 0.4, 1, 
                      1, 0.4, 1, 0.75, 0.25, 1, 1, 1, 1, 1, 0.15, 0.15, 0.15, 
                      0.15, 1, 1, 1, 0, 0, 0, 0, 0, 0),
                      
              'W_C': (0, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1.5, 
                      1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 
                      0, 0),
                      
              'LEAF_WIDTH': (0, 0.02, 0.02, 0.1, 0.02, 0.05, 0.1, 0.15, 0.15, 
                             0.15, 0.15, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
                             0.1, 0.15, 0.02, 0.05, 0.05, 0.05, 0.02, 0.001, 
                             0.05, 0.1, 0.05, 0.02, 0.1, 0.1, 0.02, 0, 0, 0, 
                             0, 0, 0),
                             
              'X_LAD': (0, 0.5, 0.5, 1, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                        1, 1, 1, 0.8, 0.5, 1, 1, 1, 0.5, 1, 1, 1, 1, 0.5, 1, 
                        1, 0.5, 0, 0, 0, 0, 0, 0)
                }

                  
def create_static_maps(land_cover_file, output_basename):
    # Open landcover image
    infid =  gdal.Open(land_cover_file, gdal.GA_ReadOnly)
    prj = infid.GetProjection()
    geo =  infid. GetGeoTransform()
    lc_data = infid.GetRasterBand(1).ReadAsArray().astype(np.uint8)
    del infid
    dims = lc_data.shape
    # Retrieve each variable static value based on landcover LUT
    for variable, values in ESACCI_LUT.items():
        print('Generating static Global map for %s'%variable)
        outfile = '%s_%s.tif'%(output_basename,variable)
        image = np.zeros(dims)
        # Loop all land cover types and assing the LUT value
        for i, lc_class in enumerate(ESACCI_LC):
            image[lc_data == lc_class] = values[i]
        dtype = gdal.GDT_Float32
        if variable == 'IGBP':
            dtype = gdal.GDT_Byte
        gu.save_img(image, 
                    geo, 
                    prj, 
                    outfile, 
                    noDataValue = -99999, 
                    dtype = dtype)       
        
        del image
 
def process_landcover_for_s2_tile(s2_file_path, 
                                   land_cover_file, 
                                   output_dir = None):
   
    if not output_dir:
        output_dir = pth.dirname(land_cover_file)
    
    tile = pth.basename(s2_file_path)[38:44]    
    outfile_basename = pth.join(output_dir,'%s_ESACCI_HR.tif'%tile)
    if pth.isfile(outfile_basename):
        print('%s already processed'%outfile_basename)
        return 
        
    # Open base file
    fid = gdal.Open(s2_file_path, gdal.GA_ReadOnly)    
    gt_out = fid.GetGeoTransform()
    prj_out = fid.GetProjection()
    shape_out = fid.RasterYSize, fid.RasterXSize
    del fid
    
    # Reproject the land cover
    gu.reproject_file(land_cover_file,
                       gt_out = gt_out,
                       prj_out = prj_out,
                       shape_out = shape_out,
                       outname = outfile_basename,
                       resampling = gdal.gdalconst.GRA_NearestNeighbour)    
    
    # Create statick LUT values
    create_static_maps(outfile_basename, pth.splitext(outfile_basename)[0])
    
    # Resample to 1km
    gtNew = [gt_out[0], 1000.0, 0, gt_out[3], 0, -1000.0]
    outFile_LR = pth.join(output_dir, 
                          outfile_basename.replace('HR.tif', 'LR.tif'))
                          
    gu.resample_file(outfile_basename, 
                      gtNew , 
                      prj_out, 
                      noDataValue=-99999,
                      outFile = outFile_LR, 
                      band_id=1, 
                      resampling = gdal.gdalconst.GRA_Mode)
              
    for variable in ESACCI_LUT.keys():
        land_cover_file = outfile_basename.replace('.tif','_%s.tif'%variable)
                              
        outFile_LR = land_cover_file.replace('_HR_', '_LR_')
        
        resampling = gdal.gdalconst.GRA_Average
        if variable == 'IGBP':
            resampling = gdal.gdalconst.GRA_Mode
                     
        gu.resample_file(land_cover_file, 
                          gtNew , 
                          prj_out, 
                          noDataValue=-99999,
                          outFile = outFile_LR, 
                          band_id=1, 
                          resampling = resampling) 