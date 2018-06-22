# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 13:33:17 2018

@author: hnieto
"""

import os.path as pth
import os
import shutil
import datetime as dt
import glob
from pySentinel import sentinel_process as sen
from pySentinel import dem
from pySentinel import gdal_utils as gu
import subprocess as sp
import gdal
import osr

MAJA_BIN = 'maja'
DTM_FOLDER = '/mnt/archive/DTM'

L2BACKWARD = -1
L2INIT = 0
L2NOMINAL = 1

MODE_STR = {-1: 'L2BACKWARD',
            0: 'L2INIT',
            1: 'L2NOMINAL'}

L2_TEMPLATE = 'S2?_OPER_SSC_L2VALD_%s_*_%s.*'
L1_TEMPLATE = 'S2?_MSIL1C_%sT*_T%s_*.SAFE'

MAX_DAY_DELAY = 15

def s2filename_to_date(s2_filename, level):
    s2_filename = pth.basename(s2_filename)
    if level == 2:
        date_str = s2_filename[11:19]
    else:
        date_str = get_date_from_maja_image(s2_filename)
        
    date = dt.datetime.strptime(date_str, '%Y%m%d')
    return date

def get_date_from_maja_image(maja_file):
    date_str = pth.basename(maja_file)[29:37]
    return date_str

def get_tile_from_maja_image(maja_file):
    tile = pth.basename(maja_file)[20:25]
    return tile

    
def get_latest_maja_image(maja_files):
    dates = [int(get_date_from_maja_image(file)) for file in maja_files]
    dates = sorted(dates)
    latest_date = dates[-1]
    indices = [i for i, s in enumerate(maja_files) if str(latest_date) in s]
    latest_maja_file = maja_files[indices[0]]
    return latest_maja_file
    
def get_nearest_sentinel2_image(s2_filename, s2_folder, day_step, level = 2):
       
    # get image date
    date_new = s2filename_to_date(s2_filename, level)
    
    delay = 0
    if level==2:
        # get tile
        tile = sen.tile_from_file_name(s2_filename)
        tile = tile[1:]
        while delay <= MAX_DAY_DELAY:
            delay += abs(day_step)
            date_old = date_new + dt.timedelta(days=day_step)
            date_new = date_old
            date_old = date_old.strftime('%Y%m%d')
            print(pth.join(s2_folder, L2_TEMPLATE%(tile, date_old)))
            files = glob.glob(pth.join(s2_folder, L2_TEMPLATE%(tile, date_old)))
            if len(files) >= 1:
                return files
        return None

    else:
        # get tile
        tile = get_tile_from_maja_image(s2_filename)
        while delay <= MAX_DAY_DELAY:
            delay += abs(day_step)
            date_old = date_new + dt.timedelta(days=day_step)
            date_new = date_old
            date_old = date_old.strftime('%Y%m%d')
            print(pth.join(s2_folder, L1_TEMPLATE%(date_old, tile)))
            files = glob.glob(pth.join(s2_folder, L1_TEMPLATE%(date_old, tile)))
            if len(files) >= 1:
                return files


def prepare_maja_input(l1_input, l2_folder, workdir, mode = L2NOMINAL):
    
    l1_folder = pth.dirname(l1_input)
    symbolic_links = []	
    if mode == L2BACKWARD:
        # input data directory contains a list of l1 products
        l1_files = glob.glob(pth.join(l1_folder, '*.SAFE'))
        for l1_file in l1_files:
            print(l1_file)
            filename = pth.basename(l1_file)
            os.symlink(l1_file, pth.join(workdir, filename))
            symbolic_links.append(pth.join(workdir, filename))		
             
    elif mode == L2INIT:
        # input data director contains only one l1 product
        # copy the l1_input to the temp folder
        filename = pth.basename(l1_file)
        os.symlink(l1_file, pth.join(workdir, filename))
        symbolic_links.append(pth.join(workdir, filename))		
        
    elif mode == L2NOMINAL:
        # input data directory contains only one l1 product and only one l2 product
        # copy the l1_input to the temp folder
        filename = pth.basename(l1_input)
        os.symlink(l1_input, pth.join(workdir, filename)) 
        symbolic_links.append(pth.join(workdir, filename))
        # find nearest l2 image in time
        l2_input = get_nearest_sentinel2_image(l1_input, l2_folder, day_step=-1)
        for l2_file in l2_input:
            filename = pth.basename(l2_file)
            os.symlink(l2_file, pth.join(workdir, filename))
            symbolic_links.append(pth.join(workdir, filename))
      
    return symbolic_links
    
def maja(l1_input, l2_output, workdir, mode = L2NOMINAL):
    
    symbolic_links = prepare_maja_input(l1_input, l2_output, workdir, mode = mode)
    command = [MAJA_BIN, 
               '--mode', 
               MODE_STR[mode], 
               '--input', 
               workdir, 
               '--output',
               l2_output]
    
    print(" ".join(command))
    progress = sp.Popen(" ".join(command),
                                shell=True,
                                universal_newlines=True,
                                stdout=sp.PIPE,
                                stdin=open(os.devnull),
                                stderr=sp.STDOUT,).stdout
    for line in iter(progress.readline, ''):
        print(line)
    
    # delete input symbolic links
    for link in symbolic_links:
        os.unlink(link)
    
    if mode == L2NOMINAL:
        l2_folders = glob.glob(pth.join(workdir, 'S2A_OPER_SSC_L2????_*.DBL.DIR'))
        for l2_folder in l2_folders:
            shutil.rmtree(l2_folder)
