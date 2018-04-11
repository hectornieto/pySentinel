#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 07:51:42 2018

@author: hector
"""
import os
import numpy as np
import os.path as pth
import scipy.stats as st
import gdal, osr
import subprocess
from ecmwfapi import ECMWFDataServer
import pygrib
import datetime as dt
from pySentinel import gdal_utils as gu

TA_CODE = 167            #, '2 metre temperature')
TDEW_CODE = 168         #, '2 metre dewpoint temperature')
UWIND_CODE = 165         #, '10 metre U wind component')
VWIND_CODE = 166         #, '10 metre V wind component')
AOT550_CODE = 210207    #, 'Total Aerosol Optical Depth at 550nm')
TCWV_CODE = 137         #, 'Total column water vapour')

VAR_NAMES = {167: '2t',
             168: '2d',
             165: '10u',
             166: '10v',
             210207: 'aod550',
             137: 'tcwv'}

HOURS_ANALYSIS=(0,6,12,18)
HOURS_FORECAST=(3,9,15,21)

TIME_STEP = 3

class linear_interpolation():
    
    def __init__(self, t, t0, t1):
       
        self.coeffs=np.asarray([(t1-t)/(t1-t0),(t-t0)/(t1-t0)])
    
    def __call__(self, z_0, z_1):
        z_0, z_1 = map(np.asarray,(z_0, z_1))
        result = self.coeffs[0]*z_0 + self.coeffs[1]*z_1
        return result

def find_grib_data (gid, variable, date, time):
    # Loop all GRIB messages and find the variable and time needed,
    # Returns the data for the desired extent (lat0, lat1, lon0, lon1)
    print(variable, date, time)
    message=gid.select(paramId=variable,validityTime=int(time*100),validityDate=date)
    if len(message) == 1:
        
        scale=message[0]['scaleValuesBy']
        offset=message[0]['offsetValuesBy']
    
        data,lats,lons=message[0].data()
        data=offset+scale*data
        return data, lats, lons
    else:
        return None, None, None
    

def get_cams (date_time, 
              path = None,
              cams_type = None,
              variables = [TA_CODE, 
                           TDEW_CODE, 
                           UWIND_CODE, 
                           VWIND_CODE, 
                           AOT550_CODE, 
                           TCWV_CODE]):
    
    if type(path)==type(None):
        path=os.getcwd()
        
    params = [str(VAR_NAMES[param]) for param in variables if param != AOT550_CODE]
    params = '/'.join(params)
    
    date_obj =  date_time.date()
    date_0 = date_obj - dt.timedelta(1)
    date_1 = date_obj + dt.timedelta(1)
    
    retrievalPeriod = '%s/to/%s'%(date_0.strftime('%Y-%m-%d'), 
                                  date_1.strftime('%Y-%m-%d'))
    
    date_str=date_time.strftime('%Y%m%d')

    cams_fc = { "class": "mc",
                "dataset": "cams_nrealtime",
                "date": retrievalPeriod,
                "expver": "0001",
                "levtype": "sfc",
                "param": params,
                "step": "0/3/9",
                "stream": "oper",
                "time": "00:00:00/06:00:00/12:00:00/18:00:00",
                "type": "fc",
                "target": pth.join(path,'cams_fc_%s.grib'%date_str),
                'expect': 'any',
                }
                       
    cams_an = {"class": "mc",
               "dataset": "cams_nrealtime",
               "date": retrievalPeriod,
               "expver": "0001",
               "levtype": "sfc",
               "param": params,
               "step": "0",
               "stream": "oper",
               "time": "00:00:00/06:00:00/12:00:00/18:00:00",
               "type": "an",
               "target": pth.join(path,'cams_an_%s.grib'%date_str),
               'expect': 'any',
               }
    
    if AOT550_CODE in variables:
        cams_aot = {"class": "mc",
                   "dataset": "cams_nrealtime",
                   "date": retrievalPeriod,
                   "expver": "0001",
                   "levtype": "sfc",
                   "param": str(AOT550_CODE),
                   "step": "0/3/6/9",
                   "stream": "oper",
                   "time": "00:00:00/12:00:00",
                   "type": "fc",
                   "target": pth.join(path,'cams_fc_aot_%s.grib'%date_str),
                   'expect': 'any',
                   }
        
        if not os.path.exists(cams_aot['target']):
            download_ecmwf_data(cams_aot)
        
    if not cams_type:
        if not os.path.exists(cams_fc['target']):
            download_ecmwf_data(cams_fc)
     
        if not os.path.exists(cams_an['target']):
            download_ecmwf_data(cams_an)
        
        gid_fc = pygrib.open(cams_fc['target'])
        gid_an = pygrib.open(cams_an['target'])
        
        return gid_an, gid_fc
   
    else:
        if cams_type == 'an':
            if not os.path.exists(cams_an['target']):
                download_ecmwf_data(cams_an)
        
            gid_an = pygrib.open(cams_an['target'])
            return gid_an
                
        elif cams_type == 'fc':
            if not os.path.exists(cams_fc['target']):
                download_ecmwf_data(cams_fc)
        
            gid_fc = pygrib.open(cams_fc['target'])
            return gid_fc

        else:
            if not os.path.exists(cams_fc['target']):
                download_ecmwf_data(cams_fc)
            
            if not os.path.exists(cams_an['target']):
                download_ecmwf_data(cams_an)
            
            gid_fc = pygrib.open(cams_fc['target'])
            gid_an = pygrib.open(cams_an['target'])
            
            return gid_an, gid_fc

def get_cams_timeseries (date_time_ini,
                         date_time_end,
                         path = None,
                         cams_type = None,
                         variables = [TA_CODE, 
                                       TDEW_CODE, 
                                       UWIND_CODE, 
                                       VWIND_CODE]):
    
    if type(path)==type(None):
        path=os.getcwd()
        
    params = [str(VAR_NAMES[param]) for param in variables if param != AOT550_CODE]
    params = '/'.join(params)
    
    
    retrievalPeriod = '%s/to/%s'%(date_time_ini.strftime('%Y-%m-%d'), 
                                  date_time_end.strftime('%Y-%m-%d'))
    
    date_str = retrievalPeriod

    cams_fc = { "class": "mc",
                "dataset": "cams_nrealtime",
                "date": retrievalPeriod,
                "expver": "0001",
                "levtype": "sfc",
                "param": params,
                "step": "3",
                "stream": "oper",
                "time": "00:00:00/06:00:00/12:00:00/18:00:00",
                "type": "fc",
                "target": pth.join(path,'cams_fc_%s.grib'%date_str),
                'expect': 'any',
                }
                       
    cams_an = {"class": "mc",
               "dataset": "cams_nrealtime",
               "date": retrievalPeriod,
               "expver": "0001",
               "levtype": "sfc",
               "param": params,
               "step": "0",
               "stream": "oper",
               "time": "00:00:00/06:00:00/12:00:00/18:00:00",
               "type": "an",
               "target": pth.join(path,'cams_an_%s.grib'%date_str),
               'expect': 'any',
               }
    
        
    if not cams_type:
        if not os.path.exists(cams_fc['target']):
            download_ecmwf_data(cams_fc)
     
        if not os.path.exists(cams_an['target']):
            download_ecmwf_data(cams_an)
        
        return cams_an['target'], cams_fc['target']
   
    else:
        if cams_type == 'an':
            if not os.path.exists(cams_an['target']):
                download_ecmwf_data(cams_an)
        
            return cams_an['target']
                
        elif cams_type == 'fc':
            if not os.path.exists(cams_fc['target']):
                download_ecmwf_data(cams_fc)
        
            return cams_fc['target']

        else:
            if not os.path.exists(cams_fc['target']):
                download_ecmwf_data(cams_fc)
            
            if not os.path.exists(cams_an['target']):
                download_ecmwf_data(cams_an)
            
            
            return cams_an['target'], cams_fc['target']
   
def download_ecmwf_data (request):
    # To run this example, you need an API key 
    # available from https://api.ecmwf.int/v1/key/
 
    server = ECMWFDataServer()
      
    server.retrieve(request)
    
def get_time_boundaries(time):
    
    # round the time down and up to TIME_STEP
    time=float(time)
    time_0 = np.floor(time/TIME_STEP) * TIME_STEP
    time_1 = np.ceil(time/TIME_STEP) * TIME_STEP
    
    return time_0, time_1

def get_ecmwf_datetime_boundaries(date, time_step):
    
    # round the time down and up to TIME_STEP
    date_0 = date + dt.timedelta(hours = time_step)
    date_ecmwf = date_0.date()
    date_ecmwf = int(date_ecmwf.strftime('%Y%m%d'))
    

    return date_ecmwf




def get_ecmwf_resample_and_interpolate(date_time, 
                                       gt_out, 
                                       prj_out, 
                                       shape_out,
                                       outfile_basename, 
                                       variables =[TA_CODE, 
                                                   TDEW_CODE, 
                                                   UWIND_CODE, 
                                                   VWIND_CODE]):
    
    path = pth.dirname(outfile_basename)
    
    gid_an, gid_fc = get_cams (date_time, 
                              path = path,
                              cams_type = None,
                              variables = variables)
    
    
    time = float(date_time.hour) + float(date_time.minute)/60.  + float(date_time.second)/3600.
    time_0, time_1 = get_time_boundaries(time)
    time_0, time_1 = map(int, [time_0, time_1])
    date_ecmwf_0 = get_ecmwf_datetime_boundaries(date_time, time_0)
    date_ecmwf_1 = get_ecmwf_datetime_boundaries(date_time, time_1)
    
    time_interp = linear_interpolation(time, time_0, time_1)
    
    for var in variables:
        if time_0 in HOURS_ANALYSIS and var != AOT550_CODE:
            data_0, lats, lons = find_grib_data (gid_an, var, date_ecmwf_0, time_0)
        else:
            data_0, lats, lons = find_grib_data (gid_fc, var, date_ecmwf_0, time_0)
        
        if time_1 in HOURS_ANALYSIS and var != AOT550_CODE:
            data_1, lats, lons = find_grib_data (gid_an, var, date_ecmwf_1, time_1)
        else:
            data_1, lats, lons = find_grib_data (gid_fc, var, date_ecmwf_1, time_1)
        
        data = time_interp(data_0, data_1)
        
        res_lats = np.mean(lats[1:,:] - lats[:-1,:])
        res_lons = np.mean(lons[:,1:] - lons[:,:-1])
        
        gt_in = [lons[0,0], res_lons, 0,
                 lats[0,0], 0, res_lats]
        
        prj_in = osr.SpatialReference()
        prj_in.ImportFromEPSG(4326)
        prj_in = prj_in.ExportToWkt()
        outfile = pth.join(outfile_basename, '%s.img'%VAR_NAMES[var])
        
        gu.resample_array (data,
                   gt_in,
                   prj_in,
                   gt_out = gt_out,
                   prj_out = prj_out,
                   shape_out = shape_out,
                   outname = outfile)