# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 10:29:31 2018

@author: Hector Nieto, hector.nieto.solana@gmail.com 
adapted from gdal_merge.py by Frank Warmerdam, warmerdam@pobox.com
"""

import math
import sys
import time
import gdal
import numpy as np
import osr
import os.path as pth

try:
    progress = gdal.TermProgress_nocb
except:
    progress = gdal.TermProgress

__version__ = '$id$'[5:-1]
verbose = 0
quiet = 0

def prj_to_epsg(prj):
    src = osr.SpatialReference()
    src.ImportFromWkt(prj)
    epsg = int(src.GetAttrValue('AUTHORITY',1))
    return epsg

def epsg_to_prj(epsg):
    src = osr.SpatialReference()
    src.ImportFromEPSG(epsg)
    prj = src.ExportFromWkt()
    return prj

def resample_array (data,
                   gt_in,
                   prj_in,
                   gt_out = None,
                   prj_out = None,
                   shape_out = None,
                   outname = 'MEM',
                   resampling = gdal.gdalconst.GRA_Bilinear):
    
    infid = save_img (data, 
                      gt_in, 
                      prj_in, 
                      'MEM', 
                      noDataValue = np.nan, 
                      dtype=gdal.GDT_Float32)
    if not gt_out:
        gt_out = gt_in
        
    if not prj_out:
        prj_out = prj_in
    
    if not shape_out:
        shape_out = data.shape
    
    
    outfid = save_img(np.empty(shape_out)*np.nan, gt_out, prj_out, outname)
    gdal.ReprojectImage(infid, outfid, prj_in, prj_out, resampling)
    infid = None    

    data = outfid.GetRasterBand(1).ReadAsArray()
    return data
    
def tile_from_file_name(s2_file_path):

    file_base_name = pth.basename(s2_file_path)
    tile = file_base_name[38:44]
    return tile 
    
def resample_file(inFile, 
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
    fileOrig = save_img(data, gtOrig, projInfoOrig, "MEM")
    # Set Input noData value
    fileOrig.GetRasterBand(1).SetNoDataValue(noDataValue)
    shapeNew = (int(round(data.shape[0]*gtOrig[1]/gtNew[1])), 
                int(round(data.shape[1]*gtOrig[5]/gtNew[5])))

    fileNew = save_img(np.empty(shapeNew)*np.nan, gtNew, projInfoNew, outFile)
    gdal.ReprojectImage(fileOrig, fileNew, projInfoOrig, projInfoNew, resampling)
    fileOrig = None    
    
    return fileNew

def reproject_file(input_file_path,
                   gt_out = None,
                   prj_out = None,
                   shape_out = None,
                   outname = 'MEM',
                   resampling = gdal.gdalconst.GRA_NearestNeighbour):


    infid =  gdal.Open(input_file_path, gdal.GA_ReadOnly)
    prj_in = infid.GetProjection()
    gt_in =  infid. GetGeoTransform()
    data = infid.GetRasterBand(1).ReadAsArray().astype(np.uint8)

    if not gt_out:
        gt_out = gt_in
        
    if not prj_out:
        prj_out = prj_in
    
    if not shape_out:
        shape_out = data.shape
    
    
    outfid = save_img(np.empty(shape_out)*np.nan, 
                         gt_out, 
                         prj_out, 
                         outname,
                         dtype = gdal.GDT_UInt16)
                         
    gdal.ReprojectImage(infid, outfid, prj_in, prj_out, resampling)
    del infid
    del outfid    


    
def get_coordinates_image(shape, gt_in):
    image_coords = np.indices(shape)
    # Convert the image extent into geographic coordinates
    x_coords, y_coords = get_map_coordinates(image_coords[0], 
                                             image_coords[1], 
                                             gt_in)
    
    return x_coords, y_coords    

def get_map_coordinates(row,col,geoTransform):
    X=geoTransform[0]+geoTransform[1]*col+geoTransform[2]*row
    Y=geoTransform[3]+geoTransform[4]*col+geoTransform[5]*row
    return X,Y
    
def get_pixel_coordinates(X, Y, geoTransform):
    row = (Y - geoTransform[3]) / geoTransform[5]
    col =( X - geoTransform[0]) / geoTransform[1]
    return int(row), int(col)

def convert_coordinate_array(input_coordinate, input_EPSG, output_EPSG=4326): 
    ''' Coordinate conversion between two coordinate systems
    
    Parameters
    ----------
    X_in : array
        input X coordinates
    Y_in : array
        input Y coordinates
    inputEPSG : int
        EPSG coordinate code of input coordinates
    outputEPSG : int
       EPSG coordinate code of output coordinates
    Z_in : float
        input altitude, default=0
        
    Returns
    -------
    X_0 : array
        output X coordinates 
    Y_0 : array
        output X coordinates   
    '''
    from pyproj import Proj, transform
    
    inProj = Proj(init='epsg:%s'%input_EPSG)
    outProj = Proj(init='epsg:%s'%output_EPSG)
    
    X_0, Y_0 = transform(inProj, outProj, input_coordinate[0], input_coordinate[1])

    return X_0, Y_0

def convert_coordinate(input_coordinate,
                       inputEPSG,
                       outputEPSG = 4326,
                       Z_in = 0):
    ''' Coordinate conversion between two coordinate systems
    
    Parameters
    ----------
    input_coordinate : tuple
        input coordinate (x,y)
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
    X_out,Y_out,Z_out=coordTransform.TransformPoint(input_coordinate[0],
                                                    input_coordinate[1],
                                                    Z_in)
    
    # print point in EPSG 4326
    return X_out, Y_out, Z_out

def read_single_band_image(image_path):
    fid = gdal.Open(image_path, gdal.GA_ReadOnly)
    geo = fid.GetGeoTransform()
    prj = fid.GetProjection()
    data = fid.GetRasterBand(1).ReadAsArray()
    
    return data, geo, prj

def save_img (data, 
              geotransform, 
              proj, 
              outPath, 
              noDataValue = np.nan, 
              dtype=gdal.GDT_Float32):
    
    # Start the gdal driver for GeoTIFF
    if outPath == "MEM":
        driver = gdal.GetDriverByName("MEM")
        driverOpt = []
    elif pth.splitext(outPath)[1] == '.img':
        driver = gdal.GetDriverByName("ENVI")
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

###############################################################################
# $Id$
#
# Project:  InSAR Peppers
# Purpose:  Module to extract data from many rasters into one output.
# Author:   Frank Warmerdam, warmerdam@pobox.com
#
###############################################################################
# Copyright (c) 2000, Atlantis Scientific Inc. (www.atlsci.com)
# Copyright (c) 2009-2011, Even Rouault <even dot rouault at mines-paris dot org>
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
###############################################################################
# changes 29Apr2011
# If the input image is a multi-band one, use all the channels in
# building the stack.
# anssi.pekkarinen@fao.org
# =============================================================================
def raster_copy( s_fh, s_xoff, s_yoff, s_xsize, s_ysize, s_band_n,
                 t_fh, t_xoff, t_yoff, t_xsize, t_ysize, t_band_n,
                 nodata=None ):

    if verbose != 0:
        print('Copy %d,%d,%d,%d to %d,%d,%d,%d.'
              % (s_xoff, s_yoff, s_xsize, s_ysize,
             t_xoff, t_yoff, t_xsize, t_ysize ))

    if nodata is not None:
        return raster_copy_with_nodata(
            s_fh, s_xoff, s_yoff, s_xsize, s_ysize, s_band_n,
            t_fh, t_xoff, t_yoff, t_xsize, t_ysize, t_band_n,
            nodata )

    s_band = s_fh.GetRasterBand( s_band_n )
    m_band = None
    # Works only in binary mode and doesn't take into account
    # intermediate transparency values for compositing.
    if s_band.GetMaskFlags() != gdal.GMF_ALL_VALID:
        m_band = s_band.GetMaskBand()
    elif s_band.GetColorInterpretation() == gdal.GCI_AlphaBand:
        m_band = s_band
    if m_band is not None:
        return raster_copy_with_mask(
            s_fh, s_xoff, s_yoff, s_xsize, s_ysize, s_band_n,
            t_fh, t_xoff, t_yoff, t_xsize, t_ysize, t_band_n,
            m_band )

    s_band = s_fh.GetRasterBand( s_band_n )
    t_band = t_fh.GetRasterBand( t_band_n )

    data = s_band.ReadRaster( s_xoff, s_yoff, s_xsize, s_ysize,
                             t_xsize, t_ysize, t_band.DataType )
    t_band.WriteRaster( t_xoff, t_yoff, t_xsize, t_ysize,
                        data, t_xsize, t_ysize, t_band.DataType )

    return 0

# =============================================================================
def raster_copy_with_nodata( s_fh, s_xoff, s_yoff, s_xsize, s_ysize, s_band_n,
                             t_fh, t_xoff, t_yoff, t_xsize, t_ysize, t_band_n,
                             nodata ):
    try:
        import numpy as Numeric
    except ImportError:
        import Numeric

    s_band = s_fh.GetRasterBand( s_band_n )
    t_band = t_fh.GetRasterBand( t_band_n )

    data_src = s_band.ReadAsArray( s_xoff, s_yoff, s_xsize, s_ysize,
                                   t_xsize, t_ysize )
    data_dst = t_band.ReadAsArray( t_xoff, t_yoff, t_xsize, t_ysize )

    nodata_test = Numeric.equal(data_src,nodata)
    to_write = Numeric.choose( nodata_test, (data_src, data_dst) )

    t_band.WriteArray( to_write, t_xoff, t_yoff )

    return 0

# =============================================================================
def raster_copy_with_mask( s_fh, s_xoff, s_yoff, s_xsize, s_ysize, s_band_n,
                           t_fh, t_xoff, t_yoff, t_xsize, t_ysize, t_band_n,
                           m_band ):
    try:
        import numpy as Numeric
    except ImportError:
        import Numeric

    s_band = s_fh.GetRasterBand( s_band_n )
    t_band = t_fh.GetRasterBand( t_band_n )

    data_src = s_band.ReadAsArray( s_xoff, s_yoff, s_xsize, s_ysize,
                                   t_xsize, t_ysize )
    data_mask = m_band.ReadAsArray( s_xoff, s_yoff, s_xsize, s_ysize,
                                    t_xsize, t_ysize )
    data_dst = t_band.ReadAsArray( t_xoff, t_yoff, t_xsize, t_ysize )

    mask_test = Numeric.equal(data_mask, 0)
    to_write = Numeric.choose( mask_test, (data_src, data_dst) )

    t_band.WriteArray( to_write, t_xoff, t_yoff )

    return 0

# =============================================================================
def names_to_fileinfos( names ):
    """
    Translate a list of GDAL filenames, into file_info objects.

    names -- list of valid GDAL dataset names.

    Returns a list of file_info objects.  There may be less file_info objects
    than names if some of the names could not be opened as GDAL files.
    """

    file_infos = []
    for name in names:
        fi = file_info()
        if fi.init_from_name( name ) == 1:
            file_infos.append( fi )

    return file_infos

# *****************************************************************************
class file_info:
    """A class holding information about a GDAL file."""

    def init_from_name(self, filename):
        """
        Initialize file_info from filename

        filename -- Name of file to read.

        Returns 1 on success or 0 if the file can't be opened.
        """
        fh = gdal.Open( filename )
        if fh is None:
            return 0

        self.filename = filename
        self.bands = fh.RasterCount
        self.xsize = fh.RasterXSize
        self.ysize = fh.RasterYSize
        self.band_type = fh.GetRasterBand(1).DataType
        self.projection = fh.GetProjection()
        self.geotransform = fh.GetGeoTransform()
        self.ulx = self.geotransform[0]
        self.uly = self.geotransform[3]
        self.lrx = self.ulx + self.geotransform[1] * self.xsize
        self.lry = self.uly + self.geotransform[5] * self.ysize

        ct = fh.GetRasterBand(1).GetRasterColorTable()
        if ct is not None:
            self.ct = ct.Clone()
        else:
            self.ct = None

        return 1

    def report( self ):
        print('Filename: '+ self.filename)
        print('File Size: %dx%dx%d'
              % (self.xsize, self.ysize, self.bands))
        print('Pixel Size: %f x %f'
              % (self.geotransform[1],self.geotransform[5]))
        print('UL:(%f,%f)   LR:(%f,%f)'
              % (self.ulx,self.uly,self.lrx,self.lry))

    def copy_into( self, t_fh, s_band = 1, t_band = 1, nodata_arg=None ):
        """
        Copy this files image into target file.

        This method will compute the overlap area of the file_info objects
        file, and the target gdal.Dataset object, and copy the image data
        for the common window area.  It is assumed that the files are in
        a compatible projection ... no checking or warping is done.  However,
        if the destination file is a different resolution, or different
        image pixel type, the appropriate resampling and conversions will
        be done (using normal GDAL promotion/demotion rules).

        t_fh -- gdal.Dataset object for the file into which some or all
        of this file may be copied.

        Returns 1 on success (or if nothing needs to be copied), and zero one
        failure.
        """
        t_geotransform = t_fh.GetGeoTransform()
        t_ulx = t_geotransform[0]
        t_uly = t_geotransform[3]
        t_lrx = t_geotransform[0] + t_fh.RasterXSize * t_geotransform[1]
        t_lry = t_geotransform[3] + t_fh.RasterYSize * t_geotransform[5]

        # figure out intersection region
        tgw_ulx = max(t_ulx,self.ulx)
        tgw_lrx = min(t_lrx,self.lrx)
        if t_geotransform[5] < 0:
            tgw_uly = min(t_uly,self.uly)
            tgw_lry = max(t_lry,self.lry)
        else:
            tgw_uly = max(t_uly,self.uly)
            tgw_lry = min(t_lry,self.lry)

        # do they even intersect?
        if tgw_ulx >= tgw_lrx:
            return 1
        if t_geotransform[5] < 0 and tgw_uly <= tgw_lry:
            return 1
        if t_geotransform[5] > 0 and tgw_uly >= tgw_lry:
            return 1

        # compute target window in pixel coordinates.
        tw_xoff = int((tgw_ulx - t_geotransform[0]) / t_geotransform[1] + 0.1)
        tw_yoff = int((tgw_uly - t_geotransform[3]) / t_geotransform[5] + 0.1)
        tw_xsize = int((tgw_lrx - t_geotransform[0])/t_geotransform[1] + 0.5) \
                   - tw_xoff
        tw_ysize = int((tgw_lry - t_geotransform[3])/t_geotransform[5] + 0.5) \
                   - tw_yoff

        if tw_xsize < 1 or tw_ysize < 1:
            return 1

        # Compute source window in pixel coordinates.
        sw_xoff = int((tgw_ulx - self.geotransform[0]) / self.geotransform[1])
        sw_yoff = int((tgw_uly - self.geotransform[3]) / self.geotransform[5])
        sw_xsize = int((tgw_lrx - self.geotransform[0]) \
                       / self.geotransform[1] + 0.5) - sw_xoff
        sw_ysize = int((tgw_lry - self.geotransform[3]) \
                       / self.geotransform[5] + 0.5) - sw_yoff

        if sw_xsize < 1 or sw_ysize < 1:
            return 1

        # Open the source file, and copy the selected region.
        s_fh = gdal.Open( self.filename )

        return raster_copy( s_fh, sw_xoff, sw_yoff, sw_xsize, sw_ysize, s_band,
                            t_fh, tw_xoff, tw_yoff, tw_xsize, tw_ysize, t_band,
                            nodata_arg )

def gdal_merge( out_file, input_files,
            format = 'GTiff',
            ul_lr = None,
            psize = None,
            separate = False,
            copy_pct = False,
            nodata = None,
            a_nodata = None,
            create_options = [],
            pre_init=[],
            band_type = 'Float32',
            createonly = False,
            bTargetAlignedPixels = False):

    start_time = time.time()

    gdal.AllRegister()
    
    band_type = gdal.GetDataTypeByName( band_type )
    if band_type == gdal.GDT_Unknown:
        print('Unknown GDAL data type: %s' %band_type)
        return False

    if ul_lr:
        ulx = float(ul_lr[0])
        uly = float(ul_lr[1])
        lrx = float(ul_lr[2])
        lry = float(ul_lr[3])

    if len(input_files) == 0:
        print('No input files selected.')
        return False

    Driver = gdal.GetDriverByName(format)
    if Driver is None:
        print('Format driver %s not found, pick a supported driver.' % format)
        return False

    DriverMD = Driver.GetMetadata()
    if 'DCAP_CREATE' not in DriverMD:
        print('Format driver %s does not support creation and piecewise writing.\nPlease select a format that does, such as GTiff (the default) or HFA (Erdas Imagine).' % format)
        return False

    # Collect information on all the source files.
    file_infos = names_to_fileinfos( input_files )

    if ul_lr is None:
        ulx = file_infos[0].ulx
        uly = file_infos[0].uly
        lrx = file_infos[0].lrx
        lry = file_infos[0].lry

        for fi in file_infos:
            ulx = min(ulx, fi.ulx)
            uly = max(uly, fi.uly)
            lrx = max(lrx, fi.lrx)
            lry = min(lry, fi.lry)

    if psize is None:
        psize = [file_infos[0].geotransform[1],file_infos[0].geotransform[5]]

    if band_type is None:
        band_type = file_infos[0].band_type

    # Try opening as an existing file.
    gdal.PushErrorHandler( 'CPLQuietErrorHandler' )
    t_fh = gdal.Open( out_file, gdal.GA_Update )
    gdal.PopErrorHandler()

    # Create output file if it does not already exist.
    if t_fh is None:

        if bTargetAlignedPixels:
            ulx = math.floor(ulx / psize[0]) * psize[0]
            lrx = math.ceil(lrx / psize[0]) * psize[0]
            lry = math.floor(lry / -psize[1]) * -psize[1]
            uly = math.ceil(uly / -psize[1]) * -psize[1]

        geotransform = [ulx, psize[0], 0, uly, 0, psize[1]]

        xsize = int((lrx - ulx) / geotransform[1] + 0.5)
        ysize = int((lry - uly) / geotransform[5] + 0.5)


        if separate != 0:
            bands=0

            for fi in file_infos:
                bands=bands + fi.bands
        else:
            bands = file_infos[0].bands


        t_fh = Driver.Create( out_file, xsize, ysize, bands,
                              band_type, create_options )
        if t_fh is None:
            print('Creation failed, terminating gdal_merge.')
            return False

        t_fh.SetGeoTransform( geotransform )
        t_fh.SetProjection( file_infos[0].projection )

        if copy_pct:
            t_fh.GetRasterBand(1).SetRasterColorTable(file_infos[0].ct)
    else:
        if separate != 0:
            bands=0
            for fi in file_infos:
                bands=bands + fi.bands
            if t_fh.RasterCount < bands :
                print('Existing output file has less bands than the input files. You should delete it before. Terminating gdal_merge.')
                return False
        else:
            bands = min(file_infos[0].bands,t_fh.RasterCount)

    # Do we need to set nodata value ?
    if a_nodata is not None:
        for i in range(t_fh.RasterCount):
            t_fh.GetRasterBand(i+1).SetNoDataValue(a_nodata)

    # Do we need to pre-initialize the whole mosaic file to some value?
    if pre_init is not None:
        if t_fh.RasterCount <= len(pre_init):
            for i in range(t_fh.RasterCount):
                t_fh.GetRasterBand(i+1).Fill( pre_init[i] )
        elif len(pre_init) == 1:
            for i in range(t_fh.RasterCount):
                t_fh.GetRasterBand(i+1).Fill( pre_init[0] )

    # Copy data from source files into output file.
    t_band = 1

    if quiet == 0 and verbose == 0:
        progress( 0.0 )
    fi_processed = 0

    for fi in file_infos:
        if createonly != 0:
            continue

        if verbose != 0:
            print("")
            print("Processing file %5d of %5d, %6.3f%% completed in %d minutes."
                  % (fi_processed+1,len(file_infos),
                     fi_processed * 100.0 / len(file_infos),
                     int(round((time.time() - start_time)/60.0)) ))
            fi.report()

        if separate == 0 :
            for band in range(1, bands+1):
                fi.copy_into( t_fh, band, band, nodata )
        else:
            for band in range(1, fi.bands+1):
                fi.copy_into( t_fh, band, t_band, nodata )
                t_band = t_band+1

        fi_processed = fi_processed+1
        if quiet == 0 and verbose == 0:
            progress( fi_processed / float(len(file_infos))  )

    # Force file to be closed.
    t_fh = None
    return True