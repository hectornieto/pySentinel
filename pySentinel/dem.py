# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 16:48:56 2018

@author: hector
"""

import numpy as np
import os
import os.path as pth
import subprocess as sp
import zipfile
import pySentinel.sentinel_process as sen
import gdal
import osr
import glob

LONS = np.arange(-180,180, 5)
LON_ID = range(len(LONS))
LATS = np.arange(90, -90, -5)
LATS_ID = range(len(LATS))

SRTM_TILE_GEO = [-180., 5., 0., 60., 0., -5.]

GDAL_BIN = ''

URL = 'http://data.cgiar-csi.org/srtm/tiles/GeoTIFF/'
USER = 'data_public'
PSSWD = 'GDdci'
MAX_ATTEMPTS = 2

def get_map_coordinates(row,col,geoTransform):
    X=geoTransform[0]+geoTransform[1]*col+geoTransform[2]*row
    Y=geoTransform[3]+geoTransform[4]*col+geoTransform[5]*row
    return X,Y
    
def get_pixel_coordinates(X, Y, geoTransform):
    row = (Y - geoTransform[3]) / geoTransform[5]
    col =( X - geoTransform[0]) / geoTransform[1]
    return int(row), int(col)

def convert_cooordinates(input_coordinate, input_EPSG, output_EPSG=4326): 

    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(input_EPSG)

     
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(output_EPSG)
    
    ct = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    
    # transform point
    X_0, Y_0, _ = ct.TransformPoint(input_coordinate[0], input_coordinate[1], 0)

    return X_0, Y_0
    
def get_srtm(ul, lr, output_path):
    row_ul, col_ul = get_pixel_coordinates(ul[0], ul[1], SRTM_TILE_GEO)
    row_lr, col_lr = get_pixel_coordinates(lr[0], lr[1], SRTM_TILE_GEO)
    
    srtm_tiles = []    
    
    for row in range(row_ul, row_lr+1):
        row = str(row+1).zfill(2)
        for col in range(col_ul, col_lr+1):
            col = str(col+1).zfill(2)
            file_id = 'srtm_%s_%s.zip'%(col,row)
            outfile = pth.join(output_path,file_id)
            success = wget_download(pth.join(URL,file_id), 
                                    pth.join(output_path,file_id))
            
            if success:
                srtm_tiles.append(outfile)
    
    return srtm_tiles   
    
            
def safe_unzip(zip_file, extractpath='.'):
    with zipfile.ZipFile(zip_file, 'r') as zf:
        for member in zf.infolist():
            abspath = os.path.abspath(os.path.join(extractpath, member.filename))
            if abspath.startswith(os.path.abspath(extractpath)):
                zf.extract(member, extractpath)            

        
def wget_download(url,output):
    wget_command=r'wget --user=%s --password=%s'%(USER, PSSWD)
    #wget_command=r'wget '
    # Try to download the product to a maximum number of attemps
    attempt=1
    if pth.isfile(output):
        return True
        
    while attempt <= MAX_ATTEMPTS:
        wget_file=wget_command+' --output-document=%s "%s"'%(output,url)
        print('Downloading file %s'%(output))
        print(wget_file)
        proc=sp.Popen(wget_file,shell=True,stdout=sp.PIPE,
                               stdin=sp.PIPE,stderr=sp.STDOUT,universal_newlines=True)
        
        for line in iter(proc.stdout.readline, ''):
            print(line.rstrip('\r\n'))
        proc.stdout.close()
        proc.wait()
        attempt += 1

    if pth.exists(output):
        return True
    else:
        return False 
        
def reproject_dem(input_dem,
                   gt_out = None,
                   prj_out = None,
                   shape_out = None,
                   outname = 'MEM',
                   resampling = gdal.gdalconst.GRA_NearestNeighbour):


    infid =  gdal.Open(input_dem, gdal.GA_ReadOnly)
    prj_in = infid.GetProjection()
    gt_in =  infid. GetGeoTransform()
    data = infid.GetRasterBand(1).ReadAsArray().astype(np.uint8)

    if not gt_out:
        gt_out = gt_in
        
    if not prj_out:
        prj_out = prj_in
    
    if not shape_out:
        shape_out = data.shape
    
    
    outfid = sen.saveImg(np.empty(shape_out)*np.nan, 
                         gt_out, 
                         prj_out, 
                         outname,
                         dtype = gdal.GDT_UInt16)
                         
    gdal.ReprojectImage(infid, outfid, prj_in, prj_out, resampling)
    del infid
    del outfid

def process_dem_for_s2_tile(s2_file_path, dem_folder = None):
    
    input_base_file = pth.basename(s2_file_path)
    
    tile = input_base_file[38:44]    

    if not dem_folder:
        # Use the same folder of input s2 image
        dem_folder = pth.dirname(s2_file_path)

    outname = pth.join(dem_folder, '%s_DEM_HR.tif'%tile)
    
    if pth.isfile(outname):
        print('%s already processed'%outname)
        return outname
        
    # Open base file
    fid = gdal.Open(s2_file_path, gdal.GA_ReadOnly)    
    gt_out = fid.GetGeoTransform()
    prj_out = fid.GetProjection()
    input_src = osr.SpatialReference()
    input_src.ImportFromWkt(prj_out)
    input_epsg = int(input_src.GetAttrValue('AUTHORITY',1))
    shape_out = fid.RasterYSize, fid.RasterXSize
    del fid
    
    # Convert the image extent into geographic coordinates
    ul_map = gt_out[0], gt_out[3]
    lr_map = get_map_coordinates(shape_out[0], shape_out[1], gt_out)
    ul = convert_cooordinates(ul_map, input_epsg, output_EPSG=4326)
    lr = convert_cooordinates(lr_map, input_epsg, output_EPSG=4326)
    
    srtm_files = get_srtm(ul, lr, dem_folder)
    mosaic_list = []
    for zip_file in srtm_files:
        safe_unzip(zip_file, extractpath = dem_folder)
        mosaic_list.append(zip_file.replace('.zip','.tif'))
    
    if len(mosaic_list) > 1:
        raw_dem = pth.join(dem_folder, 'temp.tif')
        sen.gdal_merge(raw_dem, 
                   mosaic_list,
                   separate = False, 
                   nodata = np.nan, 
                   a_nodata = -99999)
       
    else:
        raw_dem = mosaic_list[0]
    
    reproject_dem(raw_dem,
                   gt_out = gt_out,
                   prj_out = prj_out,
                   shape_out = shape_out,
                   outname = outname,
                   resampling = gdal.gdalconst.GRA_Bilinear)

    # Remove uncompressed files
    for file in mosaic_list:
        anc_files = glob.glob(pth.splitext(file)[0]+'*')            
        [os.remove(f) for f in anc_files if '.zip' not in f] 
    
    return outname

def slope_from_dem(dem_file_path, output=None, scale = 1):
    
    if not output:
        output = dem_file_path.replace('.tif','_slope.tif')

    gdal_command = 'gdaldem slope "%s" "%s" -s %s -compute_edges'%(dem_file_path, 
                                                    output, 
                                                    scale)
                                                    
    gdal_command = pth.join(GDAL_BIN, gdal_command)
    
    proc=sp.Popen(gdal_command,
                  shell=True,
                  stdout=sp.PIPE,
                  stdin=sp.PIPE,
                  stderr=sp.STDOUT,
                  universal_newlines=True)
    
    for line in iter(proc.stdout.readline, ''):
        print(line.rstrip('\r\n'))
    proc.stdout.close()
    proc.wait()

def aspect_from_dem(dem_file_path, output=None, scale = 1):
    
    if not output:
        output = pth.splitext(dem_file_path)[0]+'_slope.tif'

    gdal_command = 'gdaldem aspect "%s" "%s" -s %s -compute_edges'%(dem_file_path, 
                                                    output, 
                                                    scale) 
                                                    
    gdal_command = pth.join(GDAL_BIN, gdal_command)
    
    proc=sp.Popen(gdal_command,shell=True,
                  stdout=sp.PIPE,
                  stdin=sp.PIPE,
                  stderr=sp.STDOUT,
                  universal_newlines=True)
    
    for line in iter(proc.stdout.readline, ''):
        print(line.rstrip('\r\n'))
    proc.stdout.close()
    proc.wait()

def process_topography_for_s2_tile(s2_file_path, dem_folder = None):
    
    if not dem_folder:
        dem_folder = pth.dirname(s2_file_path)
        
    dem_file = process_dem_for_s2_tile(s2_file_path, dem_folder)

    basename = pth.basename(dem_file)    
    
    slope_file = pth.join(dem_folder, basename.replace('DEM', 'slope'))
    slope_from_dem(dem_file, output = slope_file, scale = 1)
    
    aspect_file = pth.join(dem_folder,basename.replace('DEM', 'aspect'))
    aspect_from_dem(dem_file, output = aspect_file, scale = 1)
    
    dem_lr_file = pth.join(dem_folder,basename.replace('HR', 'LR'))
    
    fid = gdal.Open(s2_file_path, gdal.GA_ReadOnly)
    gt = fid.GetGeoTransform()
    prj = fid.GetProjection()
    gt_new = [gt[0], 1000, 0, gt[3], 0, -1000]
    
    sen.resample_GDAL(dem_file, 
                      gt_new , 
                      prj, 
                      noDataValue=-99999,
                      outFile = dem_lr_file, 
                      band_id=1, 
                      resampling = gdal.gdalconst.GRA_Average)