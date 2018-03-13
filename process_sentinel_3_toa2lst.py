# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:20:01 2017

@author: hnieto
"""

from pySentinel.sentinel_download import sentinel_configuration_download
import pySentinel.sentinel_process as sen
import os
import os.path as pth
import gdal
import osr
import glob 
import shutil

geographic_EPSG = 4326
resolution = (1000,1000)
site = 'Borden'
emis_veg = [0.98, 0.975]   
emis_soil = [0.95, 0.945]
workdir=os.getcwd()

platformname = 'Sentinel-3'
footprint_template = '"Intersects(POLYGON((%s %s,%s %s,%s %s,%s %s,%s %s)))"'
server = 'https://scihub.copernicus.eu/s3'
user = 's3guest'
password = 's3guest'
instrumentshortname = 'SLSTR'
productlevel = 'L2'
timeliness = '%22Non Time Critical%22'
orbitdirection = 'descending'

query_file = pth.join(workdir, 'SENET_S3_SLSTR_L2.txt')

l2_band_list = ('LST', 'TCWV', 'NDVI', 'biome', 'fraction', 'cloud_in', 
           'sat_azimuth_tn', 'sat_zenith_tn', 'solar_azimuth_tn', 'solar_zenith_tn', 
           'dew_point_tx', 'temperature_tx', 
           'u_wind_tx_time_1_tx', 'u_wind_tx_time_2_tx', 'u_wind_tx_time_3_tx', 'u_wind_tx_time_4_tx', 'u_wind_tx_time_5_tx', 
           'v_wind_tx_time_1_tx', 'v_wind_tx_time_2_tx', 'v_wind_tx_time_3_tx', 'v_wind_tx_time_4_tx', 'v_wind_tx_time_5_tx')

l1_band_list = ('S8_BT_in','S9_BT_in','S8_BT_io','S9_BT_io','cloud_in','cloud_io',
                'sat_zenith_tn','sat_zenith_to','total_column_water_vapour_tx',
                'dew_point_tx', 'temperature_tx', 
                'u_wind_tx_time_1_tx', 'u_wind_tx_time_2_tx', 'u_wind_tx_time_3_tx', 'u_wind_tx_time_4_tx', 'u_wind_tx_time_5_tx', 
                'v_wind_tx_time_1_tx', 'v_wind_tx_time_2_tx', 'v_wind_tx_time_3_tx', 'v_wind_tx_time_4_tx', 'v_wind_tx_time_5_tx')


def convert_cooordinates(input_coordinate, input_EPSG, output_EPSG=4326): 

    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(input_EPSG)

     
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(output_EPSG)
    
    ct = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    
    # transform point
    X_0, Y_0, _ = ct.TransformPoint(input_coordinate[0], input_coordinate[1], 0)

    return X_0, Y_0

def get_map_coordinates(row, col, geoTransform):
    X = geoTransform[0] + geoTransform[1]*col + geoTransform[2]*row
    Y = geoTransform[3] + geoTransform[4]*col + geoTransform[5]*row 
                    
    return X, Y


def get_sentinel2_date_extent(setinel2_dim_folder):
    # Get the acquisition date from the filename
    filename = pth.basename(setinel2_dim_folder)
    year = filename[11:15]
    month = filename[15:17]
    day = filename[17:19]

    date_query = '%s-%s-%s'%(year, month, day)
    
    # Get the extent date from the lai.img gdal data structure
    fid = gdal.Open(pth.join(setinel2_dim_folder, 'lai.img'))
    geo = fid.GetGeoTransform()
    prj = fid.GetProjection()
    input_src = osr.SpatialReference()
    input_src.ImportFromWkt(prj)
    input_epsg = int(input_src.GetAttrValue('AUTHORITY',1))
    n_rows = fid.RasterXSize
    n_cols = fid.RasterYSize
    
    ul = get_map_coordinates(0, 0, geo)
    lr = get_map_coordinates(n_rows, n_cols, geo)
    
    extent=(ul[0], ul[1], lr[0], lr[1])
     
    return date_query,  input_epsg, extent

def sentinel_hub_footprint(extent, input_epsg):
    
    # Convert the projected coordinates into geographic coordinates
    
    ul = convert_cooordinates((extent[0],extent[1]), input_epsg, output_EPSG=4326)
    lr = convert_cooordinates((extent[2],extent[3]), input_epsg, output_EPSG=4326)
    
    footprint = (ul[0], ul[1], lr[0], ul[1], lr[0], lr[1], ul[0], lr[1], ul[0], ul[1])

    return footprint

def write_query_file(query_file, date_query, footprint):
    
    fid = open(query_file, 'w')
    
    footprint_str = footprint_template%footprint
    footprint_str = footprint_str.replace('"','%22')
    outstring = 'platformname=%s \n'%platformname
    outstring += 'footprint=%s \n'%footprint_str
    outstring += 'server=%s \n'%server
    outstring += 'user=%s \n'%user
    outstring += 'password=%s \n'%password
    outstring += 'instrumentshortname=%s \n'%instrumentshortname
    outstring += 'productlevel=%s \n'%productlevel
    #outstring+='producttype=%s \n'%producttype
    outstring += 'timeliness=%s \n'%timeliness
    outstring += 'orbitdirection=%s \n'%orbitdirection
    beginPosition = '[%sT00:00:00.000Z TO %sT23:59:59.999Z]'%(date_query,date_query)
    endPosition = beginPosition
    outstring += 'beginPosition=%s \n'%beginPosition
    outstring += 'endPosition=%s \n'%endPosition

    fid.write(outstring)
    fid.flush()
    fid.close()
    
if __name__=='__main__':

    file_list = glob.glob(pth.join(workdir, 'S2', site, '*.data'))
    
    for s2_file in file_list:
    
        date_query, input_epsg, extent = get_sentinel2_date_extent(s2_file)
        footprint = sentinel_hub_footprint(extent, input_epsg)
        print(date_query)
        write_query_file(query_file, date_query, footprint)
        test = sentinel_configuration_download(query_file,logfile=None)
        test.parse_input_config()
        test.get_query_command()
        wget_opts = {'user': test.config_data['user'],
                     'password': test.config_data['password']}
        
        output_dir = pth.join(workdir, 'S3', site)
        if not pth.isdir(output_dir):
            os.mkdir(output_dir)
        
#==============================================================================
#         test.download_query_products(output_dir,wget_opts,logfile=None)
#         
#         date_str = date_query.replace('-','')
#         download_list =  glob.glob(pth.join(output_dir, 'S3?_SL_2_LST____%sT*.zip'%date_str)) 
#         
#         for input_file in download_list:
#             input_file=input_file.replace('.zip','')
#             output_file = sen.reproject_S3_SNAP(input_file,
#                                           input_epsg,
#                                           [1000,1000],
#                                           extent=extent,
#                                           output_file=None,
#                                           paralellism=120)
#            
#             if pth.exists(output_file+'.data'):
#                 sen.delete_DIMAP_bands(output_file, l2_band_list)
#                 #shutil.rmtree(input_file)
#                 os.remove(input_file+'.zip')
#==============================================================================
        
        # Download L1 S3 data
        test.config_data['productlevel'] = 'L1'
        test.get_query_command()
        test.download_query_products(output_dir,wget_opts,logfile=None)
        
        date_str = date_query.replace('-','')
        download_list =  glob.glob(pth.join(output_dir, 'S3?_SL_1_RBT____%sT*.zip'%date_str)) 
        
        for input_file in download_list:
            input_file=input_file.replace('.zip','')
                  
            output_file = sen.reproject_S3_SNAP(input_file,
                                          input_epsg,
                                          resolution,
                                          extent=extent,
                                          output_file=None,
                                          paralellism=120)
           
            if pth.exists(output_file+'.data'):
                sen.delete_DIMAP_bands(output_file, l1_band_list)
                shutil.rmtree(input_file)
                os.remove(input_file+'.zip')
                
                LAI = sen.resample_S2_LAI(pth.join(s2_file,'lai.img'), 
                                      extent, 
                                      input_epsg, 
                                      resolution=resolution, 
                                      out_file=pth.join(output_file+'.data', 'LAI_S2.tif'))

                sen.sentinel3_LST_processor(output_file, 
                                        LAI,    
                                        emis_veg = emis_veg,    
                                        emis_soil = emis_soil, 
                                        VALID_PIXEL = 1, 
                                        out_dir = None)