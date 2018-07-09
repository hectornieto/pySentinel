# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:20:01 2017

@author: hnieto
"""

from pySentinel.sentinel_download import sentinel_configuration_download, parse_sentinel_hub_query, wget_download
import pySentinel.sentinel_process as sen
import pySentinel.gdal_utils as gu
import os
import os.path as pth
import gdal
import osr
import glob 
import shutil

geographic_EPSG = 4326
input_epsg = 32717
resolution = (1000,1000)
site = 'AgriOlmos'
emis_veg = [0.990, 0.990]   
emis_soil = [0.969, 0.977]
workdir=os.getcwd()
s2_bands_to_resample = ['lai', 'fcover', 'B4', 'B8A', 'fapar', 'lai_cab', 'lai_cw']

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
                'sat_zenith_tn','sat_zenith_to', 'solar_azimuth_tn', 'solar_zenith_tn',
                'total_column_water_vapour_tx','dew_point_tx', 'temperature_tx', 
                'u_wind_tx_time_1_tx', 'u_wind_tx_time_2_tx', 'u_wind_tx_time_3_tx', 'u_wind_tx_time_4_tx', 'u_wind_tx_time_5_tx', 
                'v_wind_tx_time_1_tx', 'v_wind_tx_time_2_tx', 'v_wind_tx_time_3_tx', 'v_wind_tx_time_4_tx', 'v_wind_tx_time_5_tx')

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
    n_cols = fid.RasterXSize
    n_rows = fid.RasterYSize
    
    ul = gu.get_map_coordinates(0, 0, geo)
    lr = gu.get_map_coordinates(n_rows, n_cols, geo)
    
    extent=(ul[0], ul[1], lr[0], lr[1])
     
    return date_query, input_src, extent

def sentinel_hub_footprint(extent, input_src):
    
    # Convert the projected coordinates into geographic coordinates
    
    ul_x, ul_y, _ = gu.convert_coordinate((extent[0],extent[1]), 
                                           input_src)

    lr_x, lr_y, _ = gu.convert_coordinate((extent[2],extent[3]), 
                                        input_src)
    
    footprint = (ul_x, ul_y, lr_x, ul_y, lr_x, lr_y, ul_x, lr_y, ul_x, ul_y)

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
    
    # Get the list of S2 files processed
    file_list = glob.glob(pth.join(workdir, 'S2', site, '*.data'))
    
    output_dir = pth.join(workdir, 'S3', site)
    if not pth.isdir(output_dir):
        os.mkdir(output_dir)
    
    # Get the list of already processed files
    processed_files=[]
    if not pth.isfile(pth.join(output_dir,'processed_%s.txt'%site)):
        logfid=open(pth.join(output_dir,'processed_%s.txt'%site),'w')
        logfid.close()
    else:
        logfid=open(pth.join(output_dir,'processed_%s.txt'%site),'r')
        for line in logfid.readlines():
            processed_files.append(line.rstrip('\n').rstrip('\r'))
        logfid.close()
        
    for s2_file in file_list:
        tile = gu.tile_from_file_name(s2_file)
        date_query, input_src, extent = get_sentinel2_date_extent(s2_file)
        
        footprint = sentinel_hub_footprint(extent, input_src)
        print(date_query)
        write_query_file(query_file, date_query, footprint)
        test = sentinel_configuration_download(query_file, 
                                               output_dir = output_dir,
                                               logfile=None)
        test.parse_input_config()
        test.get_query_command()
        wget_opts = {'user': test.config_data['user'],
                     'password': test.config_data['password']}

        # Download L1 S3 data
        test.config_data['productlevel'] = 'L1'
        test.get_query_command()
        links, filenames, orbits = parse_sentinel_hub_query(test.query_file)
        
        for link, filename, orbit in zip(links, filenames, orbits):
            
            input_file = pth.join(output_dir,filename)
            test_filename = filename.rstrip('.zip')[16:]
            
            if test_filename not in processed_files and orbit == orbitdirection:
                
                if not pth.isfile(input_file+'.zip'):
                    wget_download(link,input_file,options_dict=wget_opts)
        
                date_str = date_query.replace('-','')
                     
                l1_file = sen.reproject_S3_SNAP(input_file,
                                              input_epsg,
                                              resolution,
                                              extent = extent,
                                              output_file = input_file + '_%s'%tile,
                                              paralellism = 120)
                
                out_filename = list(filename.rstrip('.zip'))
                out_filename[7] = '2'
                out_filename = ''.join(out_filename)
                
                if pth.exists(l1_file + '.data'):
                    sen.delete_DIMAP_bands(l1_file, l1_band_list)
                    shutil.rmtree(input_file)
                    #os.remove(input_file+'.zip')
                    outfile_dir = pth.join(output_dir,out_filename+'_%s.data'%tile)
                    if not pth.isdir(outfile_dir):
                        os.makedirs(outfile_dir)
                    
                    resampled_images = {}
                    for band in s2_bands_to_resample:
                        resampled_images[band] = sen.resample_S2_LAI(
                                                        pth.join(s2_file,
                                                                 '%s.img'%band), 
                                                        extent, 
                                                        input_src, 
                                                        resolution = resolution, 
                                                        out_file = pth.join(outfile_dir, 
                                                                          '%s_S2.img'%band))

                    sen.sentinel3_LST_processor(l1_file, 
                                            resampled_images['lai'],    
                                            emis_veg = emis_veg,    
                                            emis_soil = emis_soil, 
                                            VALID_PIXEL = 0, 
                                            out_dir = outfile_dir)
                    
                    sen.copy_bands(l1_file+'.data',
                                   outfile_dir, 
                                   band_list = ['sat_zenith_tn',
                                                'solar_azimuth_tn', 
                                                'solar_zenith_tn',
                                                'total_column_water_vapour_tx',
                                                'cloud_in'])
        
                    logfid=open(pth.join(output_dir,'processed_%s.txt'%site),'a')
                    logfid.write('\n'+test_filename)
                    logfid.flush()
                    logfid.close()
