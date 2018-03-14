# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 12:51:34 2017

@author: hnieto
"""
from __future__ import print_function
import zipfile
import os.path as pth
import subprocess as sp
import os
import osr
import glob
import gdal
from pyPro4Sail import FourSAIL as sail 
import numpy as np
from pySentinel.gdal_merge import gdal_merge
from pySentinel.sen2cor_gipp_template import get_sen2cor_template
import re

C_SEN2COR_OUTPUT_DIR = "<!SEN2COR_OUTPUT_DIR!>" # Path to the default L2A folder
sen2cor_bin = 'L2A_Process' # Path to the sen2cor binary file
gptsnap_bin = 'gpt' # Path to the SNAP gpt binary file

def point2pix(coords , gt, upperBound = False):
    ''' Convert map coordinates into pixel coordinates
    
    Parameters
    ----------
    coords : tuple or array
        map coordinates (X,Y)
    gt : list or array
        GDAL geotransform array
    upperBound : bool
        flag to get the upper pixel, default=False
    
    Returns
    -------
    col : int or array
       image column
    row : int or array
       image row
    '''
       
    X = np.asarray(coords[0])
    Y = np.asarray(coords[1])
    if not upperBound:    
        px = np.floor((X - gt[0]) / gt[1]) #x pixel
        py = np.floor((Y - gt[3]) / gt[5]) #y pixel
    else:
        px = np.ceil((X - gt[0]) / gt[1]) #x pixel
        py = np.ceil((Y - gt[3]) / gt[5]) #y pixel
    return [px.astype(int), py.astype(int)]

def aggregate_image(input_image,extent,resolution,output_file=None):
    ''' Resamples an image to a lower resolution by pixel averaging
    
    Parameters
    ----------
    input_image : string
        GDAL image full path
    extent : list or array
        output image extent (xmin, ymin, xmax, ymax)
    resolution : tuple
        output image resolution (x_size,y_size)
    output_file : string
        GDAL output image full path, if None reuse the input file with suffix _<resolution>m.tif
    ''' 
    
    if not output_file:
        output_file=pth.splitext(input_image)[0]+'_%sm.tif'%(resolution[0])
    gdal_command='gdalwarp -te %s %s %s %s -tr %s %s -r average "%s" "%s"'%(
            extent[0],extent[3],extent[2],extent[1],resolution[0],resolution[1],input_image,output_file) 
    proc=sp.Popen(gdal_command,shell=True,stdout=sp.PIPE,
                               stdin=sp.PIPE,stderr=sp.STDOUT,universal_newlines=True)
    for line in iter(proc.stdout.readline, ''):
        print(line.rstrip('\r\n'))
    proc.stdout.close()
    proc.wait()

def biophysical_SNAP(input_file,
                     calcLAI=True,
                     calcFAPAR=True,
                     calcCab=False,
                     CalcCw=False,
                     CalcFVC=False,
                     output_file=None, 
                     paralellism=1):
    
    ''' Produces biophysical L2b product using SNAP biophysicalOp module
    Requires ESA-SNAP gpt binary installed 
    
    Parameters
    ----------
    input_image : string
        BEAM-DIMAP valid L2a full path without extension
    calcLAI : bool
        flag to compute LAI
    calcFAPAR : bool
        flag to compute fAPAR
    calcCab : bool
        flag to compute Ca+b
    CalcCw : bool
        flag to compute EWT
    CalcFVC= : bool
        flag to compute fcover
    output_file : string
        BEAM-DIMAP output image full path without extension, if None input file BEAM-DIMAP folder
    paralellism : int
        number of CPS to use in parallelization
        
    Returns
    -------
    output_file : string
        BEAM-DIMAP output image full path without extension, if None input file BEAM-DIMAP folder
    
    References
    ----------
    INFO: org.esa.snap.python.gpf.PyOperatorSpi: Python operator 'S2RutOp' registere
    d (Python module: 's2_rut', class: 'S2RutOp'
    INFO: org.esa.snap.core.gpf.operators.tooladapter.ToolAdapterIO: Initializing ex
    ternal tool adapters
    Usage:
      gpt biophysicalOp [options]
    
    Description:
      The 'Biophysical Processor' operator retrieves LAI from atmospherically correc
    ted Sentinel-2 products
    
    
    Source Options:
      -Ssource=<file>    The source product.
                         This is a mandatory source.
    
    Parameter Options:
      -PcomputeCab=<boolean>       Compute Cab (Chlorophyll content in the leaf)
                                   Default value is "true".
      -PcomputeCw=<boolean>        Compute Cw (Canopy Water Content)
                                   Default value is "true".
      -PcomputeFapar=<boolean>     Compute FAPAR (Fraction of Absorbed Photosyntheti
    cally Active Radiation)
                                   Default value is "true".
      -PcomputeFcover=<boolean>    Compute FVC (Fraction of Vegetation Cover)
                                   Default value is "true".
      -PcomputeLAI=<boolean>       Compute LAI (Leaf Area Index)
      '''
    
    if not output_file:
        output_file=pth.splitext(input_file)[0]+'_biophysical'
    biophysical='%s biophysicalOp -SsourceProduct="%s.dim" -PcomputeLAI=%s -PcomputeFapar=%s -PcomputeCab=%s -PcomputeCw=%s -PcomputeFcover=%s -t "%s" -q %s'%(gptsnap_bin,input_file,calcLAI,calcFAPAR,calcCab,CalcCw,CalcFVC,output_file,str(paralellism)) 
    print('Estimating Biophysical file '+input_file)
    proc=sp.Popen(biophysical,shell=True,stdout=sp.PIPE,
                               stdin=sp.PIPE,stderr=sp.STDOUT,universal_newlines=True)
    for line in iter(proc.stdout.readline, ''):
        print(line.rstrip('\r\n'))
    proc.stdout.close()
    proc.wait()
    
    if pth.isfile(output_file+'.dim'):
        return output_file
    else:
        return None  

def resample_S2_LAI(lai_file, 
                    extent, 
                    epsg_projection, 
                    resolution = (1000,1000), 
                    out_file='MEM'):
    ''' Resamples Sentinel 2 LAI product to Sentinel 3 scale by pixel aggregation
    
    Parameters
    ----------
    lai_file : string
        path to the input Sentinel lai image
    extent : tuple
        output extent in projected coordinates (ul_x, ul_y, lr_x, lr_y)
    epsg_projection : int
        EPSG code for the output projection system
    resolution : tuple, optional
        output resolution, default 1km
    out_file : string, optional
        Output path to save the resulting resampled image
        
    Returns
    -------
    LAI : array
        resampled Leaf Area Index    
    '''
    
    geo_transform = [extent[0], resolution[0], 0, 
                                 extent[1], 0, -resolution[1]]
    
    input_src = osr.SpatialReference()
    input_src.ImportFromEPSG(epsg_projection)
    projection = input_src.ExportToWkt()
    # Resample Sentinel2 LAI file
    
    file_LAI=resample_GDAL(lai_file, 
                               geo_transform, 
                               projection,
                               outFile = out_file, 
                               band_id = 1, 
                               noDataValue=np.nan,
                               resampling = gdal.gdalconst.GRA_Average)
                               
    LAI = file_LAI.GetRasterBand(1).ReadAsArray()
    
    return LAI


def sentinel3_LST_processor(toa_file, 
                            LAI,    
                            emis_veg = [0.98, 0.975],    
                            emis_soil = [0.95, 0.945], 
                            VALID_PIXEL = 1, 
                            out_dir = None):
    ''' Computes L2 Land Surface Temperature from Sentinel 3 L1c
    Top of the Atmosphere product
    
    Parameters
    ----------
    toa_file : string
        path to the L1c S3 TOA BEAM-Dimap product
    LAI : array
        LAI array
    emis_veg : float
        endmember value for pure vegetation emissivity
    emis_soil : float
        endmember value for bare soil emissivity
    VALID_PIXEL : int
        cloud-free flag value
    out_dir : string
        path to the output directory, if None the input folder will be used
        
    Returns
    -------
    success : bool
        boolean that flags a successfull process
    '''
    
    filename = pth.basename(toa_file)    
    
    if not out_dir:
        out_dir =  pth.join(pth.dirname(toa_file), 'L2')
    
    if not pth.isdir(out_dir):
        os.mkdir(out_dir)
        
    # Read total column water vapour
    toa_file+='.data'
    cloudFile=pth.join(toa_file,'cloud_in.img')
    fid=gdal.Open(cloudFile,gdal.GA_ReadOnly)
    cloud = fid.GetRasterBand(1).ReadAsArray()
    fid = None
    valid = cloud <= VALID_PIXEL
    if not valid.any():
        print('Not valid data found')
        return False
        
    # Read total column water vapour
    waterVapourFile=pth.join(toa_file,'total_column_water_vapour_tx.img')
    fid=gdal.Open(waterVapourFile,gdal.GA_ReadOnly)
    totalColumnWaterVapour = fid.GetRasterBand(1).ReadAsArray()
    # Convert from kg/m**2 to g/cm**2
    totalColumnWaterVapour[valid] = totalColumnWaterVapour[valid] * 0.1
    totalColumnWaterVapour[~valid]=np.nan
    fid = None
    
    # Read view zenith angles
    viewZenithAngleFile=pth.join(toa_file,'sat_zenith_tn.img')
    fid=gdal.Open(viewZenithAngleFile,gdal.GA_ReadOnly)
    viewZenithAngle = fid.GetRasterBand(1).ReadAsArray()
    viewZenithAngle[~valid]=np.nan
    fid = None    
    
    # Read brightness temperatures
    BT11_File=pth.join(toa_file,'S8_BT_in.img')
    fid=gdal.Open(BT11_File,gdal.GA_ReadOnly)
    bt11 = fid.GetRasterBand(1).ReadAsArray()
    bt11[~valid]=np.nan
    fid = None    
   
    # Read brightness temperatures
    BT12_File=pth.join(toa_file,'S9_BT_in.img')
    fid=gdal.Open(BT12_File,gdal.GA_ReadOnly)
    bt12 = fid.GetRasterBand(1).ReadAsArray()
    bt12[~valid]=np.nan
    
    LAI[~valid]=np.nan

    # We assume a spherical leaf inclination distribution function
    alpha=57.3
    emissivity=np.zeros((LAI.shape[0],LAI.shape[1],2))
    for i,emis_band in enumerate(zip(emis_veg, emis_soil)):

        emissivity[:,:,i] = calc_emissivity_4SAIL(LAI,
                                  viewZenithAngle,
                                  alpha,
                                  emis_veg=emis_band[0],
                                  emis_soil=emis_band[1])
    
    brightnessTemperature=np.dstack((bt11,bt12))
    LST=calc_LST_Sobrino(brightnessTemperature, emissivity, totalColumnWaterVapour, viewZenithAngle)
    LST[~valid]=np.nan
    
    outPath=pth.join(out_dir,'%s_LST_n.tif'%(filename.replace('SL_1_RBT','SL_2_RBT')))
    saveImg (np.dstack((LST,emissivity[:,:,0],emissivity[:,:,1])), fid.GetGeoTransform(), fid.GetProjection(), outPath, noDataValue =np.nan)
    fid = None
    
    return True

def calc_emissivity_4SAIL(lai,vza,alpha,emis_veg=0.98,emis_soil=0.95,tau=0.0):
    ''' Estimates surface directional emissivity using 4SAIL Radiative Transfer Model
    
    Parameters
    ----------
    lai : array_like
        Leaf Area Index
    vza : array_like
        View Zenith Angle
    alpha : array_like
        Average leaf angle in Campbell ellipsoidal leaf inclination distribution function
    emis_veg : float
        Leaf emissivity
    emis_soil : bool
        Bare soil emissivity
    tau : float
        Leaf thermal transmittance, default=0
    
    Returns
    -------
    emissivity : array_like
        surface directional emissivity
    '''
    
    # Get array dimensions and check for consistency in the dimensions
    lai=np.asarray(lai)
    vza=_check_default_parameter_size(vza, lai)
    alpha=_check_default_parameter_size(alpha, lai)
    dims=lai.shape
 
    # Vectorize the inputs
    lai=lai.reshape(-1)
    vza=vza.reshape(-1)
    alpha=alpha.reshape(-1)
    # Solar parameters (illumination angles and hotpost) in FourSAIL are not relevant:
    hotspot,sza,psi=np.zeros(lai.shape),np.zeros(lai.shape),np.zeros(lai.shape)
    
    lidf=sail.CalcLIDF_Campbell_vec(alpha,n_elements=18)
    [_,_,_,_,_,_,_,_,_,_,_,_,_,_,rdot,_,_,_,_,_,_] = sail.FourSAIL_vec(
                                    lai,
                                    hotspot,
                                    lidf,
                                    sza,
                                    vza,
                                    psi,
                                    np.ones(lai.shape)[np.newaxis,:]-emis_veg,
                                    np.zeros(lai.shape)[np.newaxis,:],
                                    np.ones(lai.shape)[np.newaxis,:]-emis_soil)
    
    # Kirchoff law
    emissivity=1.0-rdot
    # Convert output vector to original array
    emissivity=emissivity.reshape(dims)
    
    return emissivity

def calc_LST_Sobrino(BT, emissivity, totalColumnWaterVapour, viewZenithAngle):
    ''' Estimates split-window surface temperature based on AATSR coefficients
    
    Parameters
    ----------
    BT : array_like (rows,columms,2)
        Split-window TOA brightness temperature 
    emissivity : array_like (rows,colums,2)
        Split-window surface emissivity
    totalColumnWaterVapour : array_like (rows,columms)
        Total precipitable water vapour
    viewZenithAngle : array_like (rows,columms)
        View Zenith Angle 
    
    Returns
    -------
    emissivity : array_like
        surface directional emissivity
    
    References
    ----------
    .. [Sobrino2015] Sobrino JA, Jimenez-Munoz JC, Soria G, Brockmann C, Ruescas A, 
      Danne O, North P, Phillipe P, Berger M, Merchant C, Ghent D. 
      A Prototype Algorithm for Land Surface Temperature Retrieval from 
      Sentinel-3 Mission. In Sentinel-3 for Science Workshop 2015 Dec (Vol. 734, p. 38).
    '''
    
    W = totalColumnWaterVapour/np.cos(np.radians(viewZenithAngle))
    eps = (emissivity[:,:,0] + emissivity[:,:,1])/2
    deltaEps = emissivity[:,:,0] - emissivity[:,:,1]
    
    # constants from Table 1 [Sobrino2015]
    c0 = -0.268
    c1 = 1.084
    c2 = 0.2771
    c3 = 45.1
    c4 = -0.73
    c5 = -125.0
    c6 = 16.7
    
    # Equation (1) [Sobrino2015]
    LST = (BT[:,:,0] + c1*(BT[:,:,0]-BT[:,:,1]) + c2*(BT[:,:,0]-BT[:,:,1])**2 + c0 +
           (c3 + c4*W)*(1 - eps) + (c5 + c6*W)*deltaEps)
    LST[np.isnan(LST)] = 0.0    
    LST[LST<0] = 0.0
    return LST	

def coordinate_convert(X_in,Y_in,inputEPSG,outputEPSG,Z_in=0):
    ''' Coordinate conversion between two coordinate systems
    
    Parameters
    ----------
    X_in : float
        input X coordinate
    Y_in : float
        input Y coordinate
    inputEPSG : int
        EPSG coordinate code of input coordinates
    outputEPSG : int
       EPSG coordinate code of output coordinates
    Z_in : float
        input altitude, default=0
        
    Returns
    -------
    X_out : float
        output X coordinate    
    Y_out : float
        output X coordinate    
    Z_out : float
        output X coordinate   
    '''
    
    # create coordinate transformation
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(inputEPSG)

    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(outputEPSG)

    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # transform point
    X_out,Y_out,Z_out=coordTransform.TransformPoint(X_in,Y_in,Z_in)
    
    # print point in EPSG 4326
    return X_out, Y_out, Z_out

def create_binary_mask(input_array, valid_values=[1]):
    ''' Creates a binary array masking out a list of values
    
    Parameters
    ----------
    input_array : array_like
        Input array
    valid_values : list
        list of values to be masked
        
    Returns
    -------
    mask : array_like
        binary mask
    '''
    
    # Create empty array with zeros
    mask=np.zeros(input_array.shape).astype(bool)
    for value in valid_values:
        mask[input_array==value]=1
    return mask
    

def delete_DIMAP_bands(dimap_file,band_list):
    ''' Deletes a list of bands from a BEAM-DIMAP image
    
    Parameters
    ----------
    dimap_file : string
        Input BEAM-dimap file path, excluding the path
    band_list : list
        list of strings with the names of the bands to be removed
    '''
    
    dimap_bands=glob.glob(pth.join(dimap_file+'.data','*'))
    for band_file in dimap_bands:
        band_name=pth.splitext(pth.basename(band_file))[0]
        if band_name not in band_list:
            try:
                os.remove(band_file)
            except:
                print('Failed to remove %s'%(band_file))
                
def mosaic_SNAP(product_list, 
                outputfile, 
                variables=None, 
                resolution=10, 
                crs='EPSG:4326', 
                bounds=(0,40,3,43), 
                fmt='BEAM-DIMAP',
                paralellism=1):
    ''' Mosaics a list of Beam-DIMAP files. 
    Requires ESA-SNAP gpt binary installed. It is recommended to use mosaic_SAGA instead 
    
    Parameters
    ----------
    product_list : list
        list of string with full path to image files
    outputfile : string
        full path to output mosaic file
    variabes : None or list
        list of variable names to be included in the mosaic, if None all variables are included
    resolution : float
        output resolution
    crs : string
        output resolution in EPSG format, default='EPSG:4326'
    bounds : tuple
        output extent (xmin ymin xmax ymax)
    fmt : string
        output image format
    parallelism : float
        number of CPUs to be use in parallelization
    '''
    
    mosaic='%s Mosaic -Pcrs=%s -PpixelSizeX=%s -PpixelSizeY=%s -q %s'%(gptsnap_bin,crs,resolution,resolution,paralellism)
    mosaic+=' -PwestBound=%s -PsouthBound=%s -PeastBound=%s -PnorthBound=%s'%(bounds)
    if variables:
        mosaic+=' -Pvariables='
        for variable in variables:
             mosaic+='%s,'%(variable)
        mosaic=mosaic[:-1]
    mosaic+=' -t %s -f %s'%(outputfile,fmt)
    for product in product_list:
        mosaic+=' %s'%(product)
    
    print(mosaic)
    proc=sp.Popen(mosaic,shell=True,stdout=sp.PIPE,
                               stdin=sp.PIPE,stderr=sp.STDOUT,universal_newlines=True)
    for line in iter(proc.stdout.readline, ''):
        print(line.rstrip('\r\n'))
    proc.stdout.close()
    proc.wait()

def mosaic_sentinel_SAGA(product_list,
                         output_file, 
                         band_names, 
                         resolution, 
                         extent=None):
    ''' Mosaics a list Sentinel BEAM-DIMAP images. 
    Requires SAGA and GDAL installed
    
    Parameters
    ----------
    product_list : list
        list of string with full path to image files
    output_file : string
        full path to output mosaic file
    band_names : None or list
        list of variable names to be included in the mosaic, if None all variables are included
    resolution : float
        output resolution
    extent : None or tuple
        output extent (xmin ymin xmax ymax), if None use full extent of image set
        
    Returns
    -------
    output_file : string
        full path to output mosaic file in correctly processed, None othewise
    '''
    
    def mosaic_SAGA(filelist, mosaicfile,pixel_size,TYPE=7, RESAMPLING=0,OVERLAP=6,
            BLEND_DIST=10,MATCH=0):
        ''' Mosaics a list Sentinel BEAM-DIMAP images. 
        Requires SAGA and GDAL installed
        
        Parameters
        ----------
        filelist : list
            list of string with full path to image files
        mosaicfile : string
            full path to output mosaic file
        pixel_size : float
            output resolution
        extent : None or tuple
            output extent (xmin ymin xmax ymax), if None use full extent of image set
        TYPE : int
            output data type, default=7 4byte floating point
        RESAMPLING : int
            resampling method, default=0, nearest neighbour
        OVERLAP : int
            mosacking method in overlaping areas, default=6, feathering
        BLEND_DIST : int
            Blending Distance, default=10
        MATCH : int
            Default matching method, defautl=0, None
            
        Returns
        -------
        output_file : string
            full path to output mosaic file in correctly processed, None othewise
    
        References
        ----------
        -TYPE:<str>             	Preferred data storage type
        	[0] 1 bit
        	[1] 1 byte unsigned integer
        	[2] 1 byte signed integer
        	[3] 2 byte unsigned integer
        	[4] 2 byte signed integer
        	[5] 4 byte unsigned integer
        	[6] 4 byte signed integer
        	[7] 4 byte floating point
        	[8] 8 byte floating point
        	Default: 7
          -INTERPOL:<str>         	Interpolation
        	[0] Nearest Neighbor
        	[1] Bilinear Interpolation
        	[2] Inverse Distance Interpolation
        	[3] Bicubic Spline Interpolation
        	[4] B-Spline Interpolation
        	Default: 0
          -OVERLAP:<str>          	Overlapping Areas
        	[0] first
        	[1] last
        	[2] minimum
        	[3] maximum
        	[4] mean
        	[5] blend boundary
        	[6] feathering
        	Default: 1
          -BLEND_DIST:<str>       	Blending Distance
        	Minimum: 0.000000
        	Default: 10.000000
          -MATCH:<str>            	Match
        	[0] none
           Returns
        -------
        output_file : string
            full path to output mosaic file in correctly processed, None othewise
     	[1] regression
        	Default: 0
        '''  
        
        # Mosaic each band
        basedir=pth.dirname(mosaicfile)
        if not pth.isdir(basedir):
            os.mkdir(basedir)
        tmpdir=basedir+'/tmp/'
        if not pth.isdir(tmpdir):
            os.mkdir(tmpdir)
        
        gridlist=''
        for i,infile in enumerate(filelist):
            # Convert the band to a SAGA grid.
            workdir=pth.dirname(infile)
            filename=tmpdir+'/'+pth.splitext(pth.basename(infile))[0]+'_%s.sdat'%(i)
            if pth.isfile(infile):
                print('Converting file ' +filename +' to SAGA raster')
                command='gdal_translate -of SAGA "'+ infile +'" "'+ filename+'"'
                #command='saga_cmd io_gdal 0 -GRIDS="'+filename+'" -FILES="'+file+'"'
                proc=sp.Popen(command,cwd=workdir,shell=True,stdout=sp.PIPE,stdin=sp.PIPE,stderr=sp.STDOUT,universal_newlines=True)
                proc.communicate()
                filename=filename.replace('.sdat','.sgrd')
                if pth.exists(filename):
                    gridlist=gridlist+filename+';'
        
        gridlist=gridlist[:-1]
        # Mosaick each band
        # Run Saga: GRIDS, TYPE, INTERPOL,OVERLAP,BLEND_DIST,MATCH,OUTPUT_EXTENT,TARGET_USER_SIZE,TARGET_USER_FITS 
        saga_opt=(gridlist,TYPE,OVERLAP,BLEND_DIST,MATCH,RESAMPLING,pixel_size,mosaicfile)
        command='saga_cmd grid_tools 3 -GRIDS="%s" -TYPE=%s -OVERLAP=%s -BLEND_DIST=%s -MATCH=%s -RESAMPLING=%s -TARGET_USER_SIZE=%s -TARGET_OUT_GRID="%s"' %saga_opt
        proc=sp.Popen(command,cwd=tmpdir,shell=True,stdout=sp.PIPE,stdin=sp.PIPE,stderr=sp.STDOUT,universal_newlines=True)
        proc.communicate()
    
        if pth.isfile(mosaicfile+'.sdat'):
            print('Mosaic Successfully created  ' +mosaicfile)
        else:
            print('ERROR! creating  ' +mosaicfile)
            
        for infile in glob.glob(tmpdir+'/*'):os.remove(infile)
        os.rmdir(tmpdir)  
    
    mosaiclist=''
    mosaic_list=[]
    # Make temporal folder:
    workdir=pth.dirname(output_file)
    filename=pth.splitext(pth.basename(output_file))[0]
    tmp_dir=pth.join(workdir,'tmp')
    if not pth.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    # Mosaic each reflectance band
    for band_name in band_names:
        filelist=[]
        tmpmosaicfile=pth.join(tmp_dir,'%s_%s'%(filename,band_name))  
        
        for infile in product_list:
            filelist.append(pth.join(infile+'.data',band_name+'.img'))
        if len(filelist)>1:
            mosaic_SAGA(filelist, tmpmosaicfile,resolution,TYPE=7, 
                                      RESAMPLING=1,OVERLAP=4,BLEND_DIST=10,MATCH=0) 
            
        mosaic_list.append(tmpmosaicfile+'.sdat')
    mosaiclist=mosaiclist[:-1]
    # Run gdal merge to stack all single band mosaics (filelist,PCT,SEPARATE,RTYPE,OUTPUT
    print('Stacking all mosaic bands in a sigle file')
    if extent:
        gdal_merge( output_file, mosaic_list, ul_lr = extent, separate = True, nodata = np.nan, a_nodata = -99999)

    else:
        gdal_merge( output_file, mosaic_list, separate = True,nodata = np.nan, a_nodata = -99999)
   
    
    tmpfiles=glob.glob(pth.join(workdir,'tmp','*'))
    for file in tmpfiles:
        os.remove(file)
    if pth.exists(output_file):    
        return output_file
    else:
        return None
    

    
def reproject_S3_SNAP(input_file,
                      epsg,
                      resolution,
                      extent,
                      output_file=None,
                      paralellism=1):
    ''' Reprojects Sentinel 3 product using ESA-SNAP
    
    Parameters
    ----------
    input_file : string
        full path of input file
    epsg : int
        EPSG code of output coordinate system
    resolution : tuple
        output resolution in coordinate system units (x_size,y_size)
    extent : tuple
        output image extent (UL_x, UL_y, LR_x, LR_y)  
    output_file : string
        output file path
    paralellism : int
        number of CPUs to be used in parallelization
    
    Returns
    -------
    output_file : string
        full path to output mosaic file in correctly processed, None othewise
            
    References
    ----------
    INFO: org.esa.snap.python.gpf.PyOperatorSpi: Python operator 'S2RutOp' registere
    d (Python module: 's2_rut', class: 'S2RutOp')
    INFO: org.esa.snap.core.gpf.operators.tooladapter.ToolAdapterIO: Initializing ex
    ternal tool adapters
    Usage:
      gpt reproject [options]
    
    Description:
      Reprojection of a source product to a target Coordinate Reference System.
     
    Source Options:
      -cd =<file>    The source product will be collocated with this product.
      -Ssource=<file>           The product which will be reprojected. This is a mandatory source.
    
    Parameter Options:
      -PaddDeltaBands=<boolean>           Whether to add delta longitude and latitude bands.
                                          Default value is 'false'.
      -Pcrs=<string>                      A text specifying the target Coordinate Reference System
                                          Examples: EPSG:4326, AUTO:42001
      -Peasting=<double>                  The easting of the reference pixel.
      -PelevationModelName=<string>       The name of the elevation model for the or
    thorectification. If not given tie-point data is used.
      -Pheight=<integer>                  The height of the target product.
      -PincludeTiePointGrids=<boolean>    Whether tie-point grids should be includedin the output product.
                                          Default value is 'true'.
      -PnoDataValue=<double>              The value used to indicate no-data.
      -Pnorthing=<double>                 The northing of the reference pixel.
      -Porientation=<double>              The orientation of the output product (in degree).
                                          Valid interval is [-360,360].
                                          Default value is '0'.
      -Porthorectify=<boolean>            Whether the source product should be orthorectified. (Not applicable to all products)
                                          Default value is 'false'.
      -PpixelSizeX=<double>               The pixel size in X direction given in CRS
     units.
      -PpixelSizeY=<double>               The pixel size in Y direction given in CRS
     units.
      -PreferencePixelX=<double>          The X-position of the reference pixel.
      -PreferencePixelY=<double>          The Y-position of the reference pixel.
      -Presampling=<string>               The method used for resampling of floating-point raster data.
                                          Value must be one of 'Nearest', 'Bilinear', 'Bicubic'.
                                          Default value is 'Nearest'.
      -PtileSizeX=<integer>               The tile size in X direction.
      -PtileSizeY=<integer>               The tile size in Y direction.
      -Pwidth=<integer>                   The width of the target product.
      -PwktFile=<file>                    A file which contains the target Coordinate Reference System in WKT format.
    '''
    
    easting=extent[0]
    northing=extent[1]
    width=int(np.ceil((extent[2]-extent[0])/resolution[0]))
    height=int(np.ceil((extent[1]-extent[3])/resolution[1]))

    # SNAP-GPT graph file to enforce reprojecting the 1km dataset
    xml_text='''<graph id="Graph">
                  <version>1.0</version>
                  <node id="ReadNode">
                      <operator>Read</operator>
                      <sources/>
                          <parameters>
                              <file>${slstrFile}</file>
                              <formatName>Sen3_SLSTRL1B_1km</formatName>
                          </parameters>
                  </node>
                  <node id="ReprojNode">
                      <operator>Reproject</operator>
                      <sources>
                          <source>ReadNode</source>
                      </sources>
                      <parameters>
                          <crs>EPSG:%s</crs>
                          <pixelSizeX>%s</pixelSizeX>
                          <pixelSizeY>%s</pixelSizeY>
                          <resampling>Bilinear</resampling>
                          <easting>%s</easting>
                          <northing>%s</northing>
                          <width>%s</width>
                          <height>%s</height>
                          <referencePixelX>0</referencePixelX>
                          <referencePixelY>0</referencePixelY>
                      </parameters>  
                  </node>
              </graph>'''%(epsg,resolution[0],resolution[1],easting,northing,width,height)
    # Uncompress data
    basedir=pth.dirname(input_file)
    graph_file=pth.join(basedir,'graph.xml')
    fid=open(graph_file,'w')
    fid.write(xml_text)
    fid.flush()
    fid.close()
    
    if not pth.isdir(input_file):
        print('Uncompressing '+input_file)
        safe_unzip(input_file+'.zip',extractpath=basedir)
        
    input_xml_file=pth.join(input_file,'xfdumanifest.xml')
    
    reproject='%s %s -PslstrFile="%s" -q %s'%(gptsnap_bin,graph_file, input_xml_file,paralellism)
    
    if not output_file:
        output_file=input_file+'_reproject'
    
    reproject+=' -t "%s.dim"'%(output_file)
     
    print('Reprojecting file '+input_file)
    proc=sp.Popen(reproject,shell=True,stdout=sp.PIPE,
                               stdin=sp.PIPE,stderr=sp.STDOUT,universal_newlines=True)
    for line in iter(proc.stdout.readline, ''):
        print(line.rstrip('\r\n'))
    proc.stdout.close()
    proc.wait()
    
    if pth.isfile(output_file+'.dim'):
        return output_file
    else:
        return None
    
def resample_SNAP(input_file,
                  resolution,
                  upsampling='Bilinear',
                  downsampling='Mean',
                  output_file=None,
                  paralellism=8):
    
    ''' Resamples Sentinel image using ESA-SNAP module
    
    Parameters
    ----------
    input_file : string
        full path of input file
    resolution : tuple
        output resolution in coordinate system units (x_size,y_size)
    upsampling : string
        resampling method to a finer resolution, default Bilinear 
    downsampling : string
        resampling method to a coarser resolution, default Mean 
    output_file : string
        output file path
    paralellism : int
        number of CPUs to be used in parallelization

    Returns
    -------
    output_file : string
        full path to output mosaic file in correctly processed, None othewise
            
    References
    ----------
    INFO: org.esa.snap.python.gpf.PyOperatorSpi: Python operator 'S2RutOp' registere
    d (Python module: 's2_rut', class: 'S2RutOp', 
    INFO: org.esa.snap.core.gpf.operators.tooladapter.ToolAdapterIO: Initializing ex
    ternal tool adapters
    Usage:
      gpt resample [options]
    
    Description:
      Resampling of a multi-size source product to a single-size target product.
    
    
    Source Options:
      -SsourceProduct=<file>    The source product which is to be resampled.
                                This is a mandatory source.
    
    Parameter Options:
      -Pdownsampling=<string>                The method used for aggregation (downsampling to a coarser resolution).
                                             Value must be one of 'First', 'Min', 'Max', 'Mean', 'Median'.
                                             Default value is 'First'.
      -PflagDownsampling=<string>            The method used for aggregation (downsampling to a coarser resolution) of flags.
                                             Value must be one of 'First', 'FlagAnd', 'FlagOr', 'FlagMedianAnd', 'FlagMedianOr'.
                                             Default value is 'First'.
      -PreferenceBand=<string>               The name of the reference band. All other bands will be re-sampled to match its size and resolution. 
                                              Either this or targetResolutionor targetWidth and targetHeight must be set.
      -PresampleOnPyramidLevels=<boolean>    This setting will increase performance when viewing the image, 
                                              but accurate resamplings are only retrieved when zooming in on a pixel.
                                             Default value is 'true'.
      -PtargetHeight=<integer>               The height that all bands of the target product shall have. 
                                              If this is set, targetWidth must be set, too. Either this and targetWidth or referenceBand or targetResolution must be set.
      -PtargetResolution=<integer>           The resolution that all bands of the target product shall have. 
                                              The same value will be applied to scale image widths and heights. 
                                              Either this or referenceBand or targetwidth and targetHeight must be
      -PtargetWidth=<integer>                The width that all bands of the target product shall have. 
                                              If this is set, targetHeight must be set, too. 
                                              Either this and targetHeight or referenceBand or targetResolution must be set.
      -Pupsampling=<string>                  The method used for interpolation (upsampling to a finer resolution).
                                             Value must be one of 'Nearest', 'Bilinear', 'Bicubic'.
                                             Default value is "Nearest"
    '''
    
    if not output_file:
        output_file=pth.splitext(input_file)[0]+'_resample_%sm'%(resolution)
    resample='%s resample -SsourceProduct="%s" -PtargetResolution=%s -Pupsampling=%s -Pdownsampling=%s -t "%s" -q %s'%(gptsnap_bin,input_file,resolution,upsampling,downsampling,output_file,paralellism) 
    print('Resampling file '+input_file)
    proc=sp.Popen(resample,shell=True,stdout=sp.PIPE,
                               stdin=sp.PIPE,stderr=sp.STDOUT,universal_newlines=True)
    for line in iter(proc.stdout.readline, ''):
        print(line.rstrip('\r\n'))
    proc.stdout.close()
    proc.wait()
    
    if pth.isfile(output_file+'.dim'):
        return output_file
    else:
        return None    

def resample_GDAL(inFile, 
                  gtNew , 
                  projInfoNew, 
                  noDataValue=-99999,
                  outFile = "MEM", 
                  band_id=1, 
                  resampling = gdal.gdalconst.GRA_NearestNeighbour):
    
    fid=gdal.Open(inFile, gdal.GA_ReadOnly)
    gtOrig=fid.GetGeoTransform()
    projInfoOrig=fid.GetProjection()
    data=fid.GetRasterBand(band_id).ReadAsArray()
    fid = None
    fileOrig = saveImg(data, gtOrig, projInfoOrig, "MEM")
    # Set Input noData value
    fileOrig.GetRasterBand(1).SetNoDataValue(noDataValue)
    shapeNew = (round(data.shape[0]*gtOrig[1]/gtNew[1]), round(data.shape[1]*gtOrig[5]/gtNew[5]))   
    fileNew = saveImg(np.empty(shapeNew)*np.nan, gtNew, projInfoNew, outFile)
    gdal.ReprojectImage(fileOrig, fileNew, projInfoOrig, projInfoNew, resampling)
    fileOrig = None    
    
    return fileNew
    
def safe_unzip(zip_file, extractpath='.'):
    with zipfile.ZipFile(zip_file, 'r') as zf:
        for member in zf.infolist():
            abspath = os.path.abspath(os.path.join(extractpath, member.filename))
            if abspath.startswith(os.path.abspath(extractpath)):
                zf.extract(member, extractpath)
                
def prepareSen2CorGippFile(gippFile, outdir):
    # Get the template xml text
    template_xml = get_sen2cor_template()
    
    with open(gippFile, 'w') as gf:
        gipp = template_xml%(outdir)
        gf.write(gipp)    # Prepare Sen2Cor GIPP file for each date
    
def sen2cor(args):
    if len(args) != 3:
        raise TypeError("runSen2Cor() takes exactly 3 arguments ("+str(len(args))+" given)")
    l1c_file = args[0]
    sen2CorGippFile = args[1]
    resolution = args[2]
    
    print("Calling Sen2Cor for %s" %l1c_file)
    command = [sen2cor_bin, "--resolution", str(resolution), "--GIP_L2A", sen2CorGippFile, l1c_file]
    print(" ".join(command))
    progress = sp.Popen(" ".join(command),
                                shell=True,
                                universal_newlines=True,
                                stdout=sp.PIPE,
                                stdin=open(os.devnull),
                                stderr=sp.STDOUT,).stdout
    for line in iter(progress.readline, ''):
        print(line)

def sen2lai_worker(args):
    sen2lai(*args)
    
def sen2lai(l2a_file, 
            resolution=10, 
            calcLAI=True,
            calcFAPAR=True, 
            calcCab=True, 
            CalcCw=True, 
            CalcFVC=False, 
            outdir=None,
            paralellism=8):
    
    if not outdir:
       outdir=pth.dirname(l2a_file)
      
    l2b_file=l2a_file.replace('.SAFE','_L2B_%sm'%(resolution))

    if l2a_file:
        l2a_resample=resample_SNAP(l2a_file,resolution,upsampling='Bilinear',downsampling='Mean',output_file=l2b_file)
        if l2a_resample:
            l2b_file=biophysical_SNAP(l2a_resample,calcLAI=calcLAI,calcFAPAR=calcFAPAR,calcCab=calcCab,CalcCw=CalcCw,CalcFVC=CalcFVC,output_file=l2b_file,paralellism=paralellism)
        else:
            print('Resampled L2A file %s already exists'%(l2a_resample))
            l2b_file=biophysical_SNAP(l2a_resample,calcLAI=calcLAI,calcFAPAR=calcFAPAR,calcCab=calcCab,CalcCw=CalcCw,CalcFVC=CalcFVC,output_file=l2b_file,paralellism=paralellism)
        
    return l2b_file
	
def saveImg (data, geotransform, proj, outPath, noDataValue = np.nan, dtype=gdal.GDT_Float32):
    
    # Start the gdal driver for GeoTIFF
    if outPath == "MEM":
        driver = gdal.GetDriverByName("MEM")
        driverOpt = []
    else:
        driver = gdal.GetDriverByName("GTiff")
        driverOpt = ['COMPRESS=DEFLATE', 'PREDICTOR=1', 'BIGTIFF=IF_SAFER']  
    
    shape=data.shape
    if len(shape) > 2:
        ds = driver.Create(outPath, shape[1], shape[0], shape[2], dtype, driverOpt)
        ds.SetProjection(proj)
        ds.SetGeoTransform(geotransform)
        for i in range(shape[2]):
            ds.GetRasterBand(i+1).WriteArray(data[:,:,i])  
            ds.GetRasterBand(i+1).SetNoDataValue(noDataValue)
    else:
        ds = driver.Create(outPath, shape[1], shape[0], 1, dtype)
        ds.SetProjection(proj)
        ds.SetGeoTransform(geotransform)
        ds.GetRasterBand(1).WriteArray(data)
        ds.GetRasterBand(1).SetNoDataValue(noDataValue)
        
    print('Saved ' +outPath )

    return ds

def _check_default_parameter_size(parameter, input_array):

    parameter = np.asarray(parameter)
    if parameter.size == 1:
        parameter = np.ones(input_array.shape) * parameter
        return np.asarray(parameter)
    elif parameter.shape != input_array.shape:
        raise ValueError(
            'dimension mismatch between parameter array and input array with shapes %s and %s' %
            (parameter.shape, input_array.shape))
    else:
        return np.asarray(parameter)
    
def _setGraphPlaceholder(graph, placeholder, replacement):
    while re.search(placeholder, graph):
        rep = _handleVaryingPlaceholder(graph, placeholder, replacement)
        graph = re.sub(placeholder, rep, graph, count=1)
    return graph

# Handle placeholders which might have a varying ending (e.g. C_OUTPUT_FILE)
def _handleVaryingPlaceholder(graph, placeholder, replacement):
    match = re.search(placeholder, graph)
    if match and match.lastindex == 1:
        replacement = os.path.splitext(replacement)[0] + match.group(2) + os.path.splitext(replacement)[1]
    return replacement
