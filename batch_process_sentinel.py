# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 14:31:25 2018

@author: hector
"""

from pySentinel import landcover, gdal_utils, dem, ecmwf, solar_irradiance
import os.path as pth
import gdal
import datetime as dt
import glob

site = 'Borden'    

workdir = '/home/hector/Data/Projects/SEN4ET/WP3'    

S2_dir = pth.join(workdir, 'S2', site)    
S3_dir = pth.join(workdir, 'S3', site , 'L2')    
lc_dir = pth.join(workdir, 'CCI_LC')
dem_dir = pth.join(workdir, 'DEM')
Synergy_dir = pth.join(workdir, 'S3S2', site)    
ecmwf_dir = pth.join(workdir, 'ECMWF', site)   
land_cover_file = pth.join(lc_dir, landcover.LC_FILENAME)
ecmwf_variables = [ecmwf.TA_CODE, ecmwf.TDEW_CODE, ecmwf.UWIND_CODE, ecmwf.VWIND_CODE]

s2_file_list = glob.glob(pth.join(S2_dir, '*.vrt'))

for s2_image in s2_file_list: 
    s2_name = pth.basename(s2_image)
    date_str = s2_name[11:19] 
    s3_file_list = glob.glob(pth.join(S3_dir, 'S3?_SL_2_RBT____%sT*_LST_n.tif'%date_str))
    
    for s3_lst_image in s3_file_list:
        S3_name = pth.basename(s3_lst_image)
        
        datetime_str = S3_name[16:31]
        acquisition_date = dt.datetime.strptime(datetime_str, '%Y%m%dT%H%M%S')
        
        # Process topography data
        dem.process_topography_for_s2_tile(s2_image, dem_folder = dem_dir)
        # Process land cover data
        landcover.process_landcover_for_s2_tile(s2_image, 
                                               land_cover_file, 
                                               output_dir = lc_dir)
        
        
        
        # Download ECMWF data
        fid = gdal.Open(s3_lst_image, gdal.GA_ReadOnly)
        prj_s3 = fid.GetProjection()
        gt_s3 = fid.GetGeoTransform()
        shape_s3 = fid.RasterYSize, fid.RasterXSize
        ecmwf_file_basename = pth.join(ecmwf_dir,S3_name.replace('_LST_n.tif', ''))
        
        ecmwf.get_ecmwf_resample_and_interpolate(acquisition_date, 
                                           gt_s3, 
                                           prj_s3, 
                                           shape_s3,
                                           ecmwf_file_basename, 
                                           variables = ecmwf_variables)
        
        
        
        # Get irradiance and incidence angle
        Sdn_out_file = s3_lst_image.replace('LST_n.tif', 'Sdn.tif')
        incidence_out_file= pth.join(Synergy_dir,
                                     pth.basename(s3_lst_image).replace('reproject_LST_n.tif', 
                                                                     'SEN2_incidence.tif'))    
        
        solar_irradiance.get_solar_irradiance_and_incidence(s3_lst_image, 
                                                       s2_image.replace('.vrt', 
                                                                        '.data'), 
                                                       dem_dir,
                                                       Sdn_out_file,
                                                       incidence_out_file)