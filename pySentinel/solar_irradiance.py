# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 14:01:44 2018

@author: hector
"""
import numpy as np
import os.path as pth
import gdal
import osr
import time
from pySentinel import gdal_utils as gu

S2_AOT_scale = 0.001
S2_TCWV_scale = 0.01 #mm

def declination_angle(doy):
    ''' Calculates the Earth declination angle
    
    Parameters
    ----------
    doy : float or int
        day of the year
    
    Returns
    -------
    declination : float
        Declination angle (radians)
    '''    
    declination = np.radians(23.45) * np.sin((2.0 * np.pi * doy / 365.0) - 1.39)
    
    return declination

def hour_angle (ftime, declination, lon, stdlon = 0): 
    '''Calculates the hour angle
    
    Parameters
    ----------
    ftime : float
        Time of the day (decimal hours)
    declination : float
        Declination angle (radians)
    lon : float
        longitude of the site (degrees).
    stdlon : float
        Longitude of the standard meridian that represent the ftime time zone    
    
    Returns
    w : float
        hour angle (radians)
    '''    

    EOT = 0.258 * np.cos(declination) - 7.416 * np.sin(declination) - \
    3.648 * np.cos(2.0 * declination) - 9.228 * np.sin(2.0 * declination)
    LC = (stdlon - lon) / 15.
    time_corr = (-EOT / 60.) + LC
    solar_time = ftime - time_corr
    # Get the hour angle
    w = np.radians((12.0 - solar_time) * 15.)
    
    return w
    
def calc_sun_angles(lat, lon, stdlon, doy, ftime):
    '''Calculates the Sun Zenith and Azimuth Angles (SZA & SAA).

    Parameters
    ----------
    lat : float
        latitude of the site (degrees).
    lon : float
        longitude of the site (degrees).
    stdlon : float
        central longitude of the time zone of the site (degrees).
    doy : float
        day of year of measurement (1-366).
    ftime : float
        time of measurement (decimal hours).

    Returns
    -------
    sza : float
        Sun Zenith Angle (degrees).
    saa : float
        Sun Azimuth Angle (degrees).

    '''

    lat, lon, stdlon, doy, ftime = map(
        np.asarray, (lat, lon, stdlon, doy, ftime))
    # Calculate declination
    declination = declination_angle(doy)
    # Get the hour angle
    w = hour_angle (ftime, declination, lon, stdlon = stdlon)
    w = np.degrees(w)
    # Get solar elevation angle
    sin_thetha = np.cos(np.radians(w)) * np.cos(declination) * np.cos(np.radians(lat)) \
                 + np.sin(declination) * np.sin(np.radians(lat))
    
    sun_elev = np.arcsin(sin_thetha)
    # Get solar zenith angle
    sza = np.pi / 2.0 - sun_elev
    sza = np.asarray(np.degrees(sza))
    # Get solar azimuth angle
    cos_phi = np.asarray((np.sin(declination) * np.cos(np.radians(lat))
                        - np.cos(np.radians(w)) * np.cos(declination) * np.sin(np.radians(lat)))
                        / np.cos(sun_elev))
    
    saa = np.zeros(sza.shape)
    saa[w <= 0.0] = 360. - np.degrees(np.arccos(cos_phi[w <= 0.0]))
    saa[w > 0.0] = np.degrees(np.arccos(cos_phi[w > 0.0]))
    return np.asarray(sza), np.asarray(saa)
    
def incidence_angle_tilted(lat, 
                           lon, 
                           doy, 
                           ftime, 
                           stdlon = 0, 
                           A_ZS = 0, 
                           slope = 0):
    ''' Calculates the incidence solar angle over a tilted flat surface
    
    Parameters
    ----------
    lat :  float or array
        latitude (degrees)
    lon :  float or array
        longitude (degrees)
    doy : int
        day of the year
    ftime : float
        Time of the day (decimal hours)
    stdlon : float
        Longitude of the standard meridian that represent the ftime time zone    
    A_ZS : float or array
        surface azimuth angle, measured clockwise from north (degrees)
    slope : float or array
        slope angle (degrees)
    
    Returns
    -------
    cos_theta_i : float or array
        cosine of the incidence angle
    '''
     
    # Get the dclination and hour angle
    delta = declination_angle(doy)
    omega = hour_angle(ftime, delta, lon, stdlon = stdlon)

    # Convert remaining angles into radians    
    lat, A_ZS, slope = map(np.radians, [lat, A_ZS, slope])    
    
    cos_theta_i = (np.sin(delta) * np.sin(lat) * np.cos(slope) 
                   + np.sin(delta) * np.cos(lat) * np.sin(slope) * np.cos(A_ZS)
                   + np.cos(delta) * np.cos(lat) * np.cos(slope) * np.cos(omega)
                   - np.cos(delta) * np.sin(lat) * np.sin(slope) * np.cos(A_ZS) * np.cos(omega)
                   - np.cos(delta) * np.sin(slope) * np.sin(A_ZS) * np.sin(omega))
    
    return cos_theta_i

def calc_global_horizontal_radiance_clear_sky(DOY,
                                              sza,
                                              aot_550, 
                                              pw, 
                                              press, 
                                              altitude=0, 
                                              solar_constant=1367.7):
                                                  
    # Calculate extraterrestrial solar irradiance
    I_0=calc_extraterrestrial_radiation(DOY, solar_constant=1367.7)
    # Calculate the air mass
    air_mass=calc_air_mass_planeparallel(sza)
    
    # Calculate the optical thickness of a water and aerosol free atmosphere
    #delta_cda=calc_delta_cda(air_mass)
    # Calculate Linke turbidity index
    tl=calc_Linke_turbidity_Ineichen(aot_550,press,pw,altitude=altitude)
    
    # Calculate Kasten coefficients
    f_h1=np.exp(-altitude/8000.0)
    f_h2=np.exp(-altitude/1250.0)    
    
    # Calculate Ineichen and Perez coefficients
    c_g1=5.09e-5*altitude+0.868
    c_g2=3.92e-5*altitude+0.0387
    
    I=c_g1*I_0*np.cos(np.radians(sza))*np.exp(-c_g2*air_mass*(f_h1+f_h2*(tl-1.0)))*np.exp(0.01*air_mass**1.8)
    
    return I

def calc_air_mass_planeparallel(sza):
    
    air_mass=1./np.cos(np.radians(sza))
    return air_mass

def calc_delta_cda(air_mass):
    delta_cda=0.128-0.054*np.log10(air_mass)
    return delta_cda
    
def calc_Linke_turbidity_Ineichen(aot_550,press,pw,altitude=0):
    
    p0_p=calc_pressure(altitude)/press
    
    tl=3.91*np.exp(0.689*p0_p)*aot_550+0.376*np.log(pw)+2+0.54*p0_p-0.5*p0_p**2+0.15*p0_p**3
    return tl
    

def calc_extraterrestrial_radiation(DOY, solar_constant=1367.7):
    I_0=solar_constant*(1.0+0.033*np.cos(2.*np.pi*DOY/365.))
    return I_0

    
def calc_pressure(z):
    ''' Calculates the barometric pressure above sea level.

    Parameters
    ----------
    z: float
        height above sea level (m).

    Returns
    -------
    p: float
        air pressure (mb).'''

    p = 1013.25 * (1.0 - 2.225577e-5 * z)**5.25588
    return np.asarray(p)

def get_solar_irradiance_and_incidence(s3_image, 
                                       s2_image_folder,
                                       dem_folder,
                                       Sdn_out_file,
                                       incidence_out_file):
    
    # Get the date and time from the image filename
    s3_base_file = pth.basename(s3_image)
    date_str = s3_base_file[16:24]
    time_s3 = float(s3_base_file[25:27]) + float(s3_base_file[27:29])/60.\
                + float(s3_base_file[29:31])/3600.
    
    date = time.strptime(date_str, '%Y%m%d')
    DOY = date.tm_yday
    
    # Get projection and extent information            
    fid = gdal.Open(s3_image, gdal.GA_ReadOnly)
    geo_s3 = fid.GetGeoTransform()
    prj_s3 = fid.GetProjection()
    s3_src = osr.SpatialReference()
    s3_src.ImportFromWkt(prj_s3)
    s3_epsg = int(s3_src.GetAttrValue('AUTHORITY',1))
    shape_s3 = fid.RasterYSize, fid.RasterXSize
    del fid
    
    # Obtain pixel coordinates
    x_s3, y_s3 = gu.get_coordinates_image(shape_s3, geo_s3)    
    lon_s3, lat_s3 = gu.convert_coordinate_array((x_s3, y_s3), 
                                               s3_epsg, 
                                               output_EPSG=4326)
    
    # Get solar angles                                   
    sza_s3, saa_s3 =calc_sun_angles(lat_s3, lon_s3, 0, DOY, time_s3)
    
    
    # Get atmospheric information from Sentinel2 tile
    fid = gdal.Open(pth.join(s2_image_folder,'mask.img'), gdal.GA_ReadOnly)
    geo_s2 = fid.GetGeoTransform()
    prj_s2 = fid.GetProjection()
    s2_src = osr.SpatialReference()
    s2_src.ImportFromWkt(prj_s2)
    s2_epsg = int(s2_src.GetAttrValue('AUTHORITY',1))
    mask = fid.GetRasterBand(1).ReadAsArray()

    fid = gdal.Open(pth.join(s2_image_folder,'quality_wvp.img'), gdal.GA_ReadOnly)
    tcwv = S2_TCWV_scale * fid.GetRasterBand(1).ReadAsArray()

    fid = gdal.Open(pth.join(s2_image_folder,'quality_aot.img'), gdal.GA_ReadOnly)
    aot_550 = S2_AOT_scale * fid.GetRasterBand(1).ReadAsArray()
    

    # Retrieve topographic data
    s2_base_file = pth.basename(s2_image_folder)
    tile = s2_base_file[38:44]
    
    fid = gdal.Open(pth.join(dem_folder,'%s_DEM_LR.tif'%tile), gdal.GA_ReadOnly)
    dem_lr = fid.GetRasterBand(1).ReadAsArray()
     
    fid = gdal.Open(pth.join(dem_folder,'%s_slope_HR.tif'%tile), gdal.GA_ReadOnly)
    slope_hr = fid.GetRasterBand(1).ReadAsArray()

    fid = gdal.Open(pth.join(dem_folder,'%s_aspect_HR.tif'%tile), gdal.GA_ReadOnly)
    aspect_hr = fid.GetRasterBand(1).ReadAsArray()
    
    # Calculate irradiance
    press = calc_pressure(dem_lr)
    Sdn = calc_global_horizontal_radiance_clear_sky(DOY,
                                              sza_s3,
                                              np.nanmean(aot_550[mask>0]), 
                                              np.nanmean(tcwv[mask>0]), 
                                              press, 
                                              altitude = dem_lr, 
                                              solar_constant=1367.7)
    
    
    gu.save_img(Sdn, 
                geo_s3, 
                prj_s3, 
                Sdn_out_file, 
                noDataValue = np.nan, 
                dtype=gdal.GDT_Float32)    
    

    press = calc_pressure(dem_lr)
    
    x_s2, y_s2 = gu.get_coordinates_image(mask.shape, geo_s2)    
    lon_s2, lat_s2 = gu.convert_coordinate_array((x_s2, y_s2), 
                                               s2_epsg, 
                                               output_EPSG=4326)
    
    cos_theta_i = incidence_angle_tilted(lat_s2, 
                                     lon_s2, 
                                     DOY, 
                                     time_s3, 
                                     stdlon = 0, 
                                     A_ZS = aspect_hr, 
                                     slope = slope_hr)

    gu.save_img(cos_theta_i, 
                geo_s2, 
                prj_s2, 
                incidence_out_file, 
                noDataValue = np.nan, 
                dtype=gdal.GDT_Float32)

    return Sdn, sza_s3, saa_s3, press, cos_theta_i