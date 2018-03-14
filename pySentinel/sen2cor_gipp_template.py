# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 09:06:25 2018

@author: hnieto
"""

def get_sen2cor_template():
    
    sen2Cor_template_xml = '''<?xml version="1.0" encoding="UTF-8"?>
    <Level-2A_Ground_Image_Processing_Parameter xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="L2A_GIPP.xsd">
      <Common_Section>
        <Log_Level>INFO</Log_Level>
        <!-- can be: NOTSET, DEBUG, INFO, WARNING, ERROR, CRITICAL -->
        <Nr_Processes>7</Nr_Processes>
        <!-- can be an unsigned integer value specifying the number or processes you intend to operate in parallel or: AUTO. If AUTO is chosen, the processor determines the number of processes automatically, using cpu_count() -->
        <Target_Directory>%s</Target_Directory>
        <!-- should be either a directory or 'DEFAULT'. If default, target will be created at root of L1C product -->
        <DEM_Directory>~/DEM</DEM_Directory>
        <!-- should be either a directory in the sen2cor home folder or 'NONE'. If NONE, no DEM will be used -->
        <DEM_Reference>http://data_public:GDdci@data.cgiar-csi.org/srtm/tiles/GeoTIFF/</DEM_Reference>
        <!-- will be ignored if DEM is NONE. A SRTM DEM will be downloaded from this reference, if no local DEM is available -->
        <PSD_Scheme PSD_Version="13" PSD_Reference="S2-PDGS-TAS-DI-PSD-V13.1_Schema">
    		<UP_Scheme_1C>S2_User_Product_Level-1C_Metadata.xsd</UP_Scheme_1C>
    		<UP_Scheme_2A>S2_User_Product_Level-2A_Metadata.xsd</UP_Scheme_2A>
    		<Tile_Scheme_1C>S2_PDI_Level-1C_Tile_Metadata.xsd</Tile_Scheme_1C>
    		<Tile_Scheme_2A>S2_PDI_Level-2A_Tile_Metadata.xsd</Tile_Scheme_2A>
    		<DS_Scheme_1C>S2_PDI_Level-1C_Datastrip_Metadata.xsd</DS_Scheme_1C>
    		<DS_Scheme_2A>S2_PDI_Level-2A_Datastrip_Metadata.xsd</DS_Scheme_2A>
    	</PSD_Scheme>
        <PSD_Scheme PSD_Version="14" PSD_Reference="S2-PDGS-TAS-DI-PSD-V14.2_Schema">
    		<UP_Scheme_1C>S2_User_Product_Level-1C_Metadata.xsd</UP_Scheme_1C>
    		<UP_Scheme_2A>S2_User_Product_Level-2A_Metadata.xsd</UP_Scheme_2A>
    		<Tile_Scheme_1C>S2_PDI_Level-1C_Tile_Metadata.xsd</Tile_Scheme_1C>
    		<Tile_Scheme_2A>S2_PDI_Level-2A_Tile_Metadata.xsd</Tile_Scheme_2A>
    		<DS_Scheme_1C>S2_PDI_Level-1C_Datastrip_Metadata.xsd</DS_Scheme_1C>
    		<DS_Scheme_2A>S2_PDI_Level-2A_Datastrip_Metadata.xsd</DS_Scheme_2A>
    	</PSD_Scheme>
        <GIPP_Scheme>L2A_GIPP.xsd</GIPP_Scheme>
        <SC_Scheme>L2A_CAL_SC_GIPP.xsd</SC_Scheme>
        <AC_Scheme>L2A_CAL_AC_GIPP.xsd</AC_Scheme>
      </Common_Section>
      <Scene_Classification>
        <Filters>
          <Median_Filter>0</Median_Filter>
        </Filters>
      </Scene_Classification>
      <Atmospheric_Correction>
        <Look_Up_Tables>
          <Aerosol_Type>AUTO</Aerosol_Type>
          <!-- RURAL, MARITIME, AUTO -->
          <Mid_Latitude>AUTO</Mid_Latitude>
          <!-- SUMMER, WINTER, AUTO -->
          <Ozone_Content>0</Ozone_Content>
          <!-- 0, f-k, t-y -->
          <!-- The atmospheric temperature profile and ozone content:
          	"0" means: get best approximation from metadata (this is smallest difference between metadata and column DU)
          	
            For midlatitude summer atmosphere:
            "f" 250 DU
            "g" 290 DU
            "h" 331 DU (standard MS)
            "i" 370 DU
            "j" 410 DU
            "k" 450 DU
            
            For midlatitude winter atmosphere:
            "t" 250 DU
            "u" 290 DU
            "v" 330 DU
            "w" 377 DU (standard MW)
            "x" 420 DU
            "y" 460 DU
           -->
        </Look_Up_Tables>
        <Flags>
          <WV_Correction>1</WV_Correction>
          <!-- 0: No WV correction, 1: only 940 nm bands, 2: only 1130 nm bands , 3: both regions used during wv retrieval, 4: Thermal region -->
          <VIS_Update_Mode>1</VIS_Update_Mode>
          <!-- 0: constant, 1: variable visibility -->
          <WV_Watermask>1</WV_Watermask>
          <!-- 0: not replaced, 1: land-average, 2: line-average -->
          <Cirrus_Correction>1</Cirrus_Correction>
          <!-- 0: no, 1: yes -->
          <BRDF_Correction>0</BRDF_Correction>
          <!-- 0: no BRDF correction, 1: , 2: ,11, 12, 22, 21: -->
          <BRDF_Lower_Bound>0.22</BRDF_Lower_Bound>
          <!-- In most cases, g=0.2 to 0.25 is adequate, in extreme cases of overcorrection g=0.1 should be applied -->
        </Flags>
        <Calibration>
          <DEM_Unit>0</DEM_Unit>
          <!-- (0=[m], 1=[dm], 2=[cm]) -->
          <Adj_Km>1.000</Adj_Km>
          <!-- Adjancency Range [km] -->
          <Visibility>40.0</Visibility>
          <!-- visibility (5 <= visib <= 120 km) -->
          <Altitude>0.100</Altitude>
          <!-- [km] -->
          <Smooth_WV_Map>100.0</Smooth_WV_Map>
          <!-- length of square box, [meters] -->
          <WV_Threshold_Cirrus>0.25</WV_Threshold_Cirrus>
          <!-- water vapor threshold to switch off cirrus algorithm [cm]Range: 0.1-1.0 -->
        </Calibration>
      </Atmospheric_Correction>
    </Level-2A_Ground_Image_Processing_Parameter>
    '''
    
    return sen2Cor_template_xml