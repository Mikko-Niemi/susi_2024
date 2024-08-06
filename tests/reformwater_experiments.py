# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 17:26:29 2023

@author: alauren
"""

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def FIN(vege=True):
    intervals = np.array([  0,  38,  73, 109, 162])
    cols = ['Column_no', 'Month','Code', 'Country', 'Site',
           'Sample_ID', 'Treat', 'WT',
           'WT_cm_from_top_of_the_column', 'Site','bd_a','bd_b','bd_c','Peat_surface_from_top_of_tube_cm',
           'T_C_mean', 'CO_flux_mg_s1_m2_mean','Vegetation','Total_dry_g', 'Vascular_photosynth_dry_g','TotalMoss_dry_g',
           'nmonth', 'pvm', 'T_soil_C',
           'TC', 'SUVA254', 'E2E3']
    df = pd.read_excel('C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/egu2023/reformwater_data.xlsx', sheet_name='data', usecols =cols)
    df['pvm'][0].date()
    #samples = np.unique(df['Column_no'].values)
    
    #-------PRE guessed values-------------------
    water_content_low = 63.0/100                         # m3 m-3 from Laurén et al. 2019
    water_content_high = 71.0/100
    
    #--------------- FOR Paroninkorpi only---------------------
    if vege:
        dfFIN = df[(df['Country']=='FIN')&(df['Treat']!='FN')&(df['Vegetation']=='Yes')].copy()
    else: 
        dfFIN = df[(df['Country']=='FIN')&(df['Treat']!='FN')&(df['Vegetation']!='Yes')].copy()
    cyls = np.unique(dfFIN['Column_no'].values)
    dfFIN['sample_height'] = 50.0 - dfFIN['Peat_surface_from_top_of_tube_cm']
    dfFIN['WT_cm_from_bottom'] = 50.0 + dfFIN['WT_cm_from_top_of_the_column']       # water table up from the bottom cm
    dfFIN['WT_cm'] = dfFIN['sample_height']- dfFIN['WT_cm_from_bottom']             # water table cm down from the sample surface
    #dfFIN['peat_volume_m3'] = dfFIN['sample_height']/100.0 * np.pi * 0.1**2         # peat volume in m3 
    dfFIN['peat_volume_m3'] = dfFIN['sample_height']/100.0 * np.pi * 0.08**2         # peat volume in m3 
    
    dfFIN['water_content'] = water_content_low
    dfFIN['water_content'][dfFIN['WT']=='H'] = water_content_high
    dfFIN['water_m3'] = dfFIN['water_content'] * dfFIN['peat_volume_m3'] 
    dfFIN['DOC_in_water_mg'] = dfFIN['TC'] * dfFIN['water_m3'] * 1000.0
    dfFIN['DOC_in_water_mg'].mean() 
    
    mean_a, mean_b, mean_c =  dfFIN[['bd_a','bd_b','bd_c']].mean().values
    
    bds_a = np.where(np.isfinite(dfFIN['bd_a'].values) , dfFIN['bd_a'].values, mean_a)
    bds_b = np.where(np.isfinite(dfFIN['bd_b'].values) , dfFIN['bd_b'].values, mean_b)
    bds_c = np.where(np.isfinite(dfFIN['bd_c'].values) , dfFIN['bd_c'].values, mean_c)
    dfFIN['bds_a'] = bds_a
    dfFIN['bds_b'] = bds_b
    dfFIN['bds_c'] = bds_c
    bds = dfFIN[['bds_a', 'bds_b','bds_c']].mean(axis = 1)
    dfFIN['bulk_density'] = bds                                                     #kg/m3
    dfFIN['peat_mass'] = dfFIN['peat_volume_m3'] * dfFIN['bulk_density']
    
    columns = np.unique(dfFIN['Column_no'].values)
    #dfFIN['plant_kg_ha-1'] = dfFIN['Total_dry_g']/(0.1**2*np.pi)*10000/1000 *1.5
    #dfFIN['plant_kg_ha-1'] = dfFIN['Vascular_photosynth_dry_g']/(0.1**2*np.pi)*10000/1000 *  5.0 \
    #    + dfFIN['TotalMoss_dry_g']/(0.1**2*np.pi)*10000/1000 
    dfFIN['plant_kg_ha-1'] = dfFIN['Vascular_photosynth_dry_g']/(0.1**2*np.pi)*10000/1000 *  5.0 \
        + dfFIN['TotalMoss_dry_g']/(0.08**2*np.pi)*10000/1000 
    
    
    dfFIN['plant_respiration_CO2'] = dfFIN['plant_kg_ha-1'] *0.5* 44/12 
    
    co2out = np.zeros(len(columns))
    WTs = np.zeros(len(columns))
    bd_as = np.zeros(len(columns))
    bd_bs = np.zeros(len(columns))
    bd_cs = np.zeros(len(columns))
    delta_DOC = np.zeros(len(columns))
    
    for n, col in enumerate(columns):
      dfa = dfFIN[dfFIN['Column_no']==col].copy()                                   
      dfa = dfa.sort_values('nmonth')                                               # sort data to chronological order   
      m_peat = dfa['peat_mass'].mean()*1000                                         # peat mass in g
      meas = dfa['CO_flux_mg_s1_m2_mean'].values*0.879               # measured CO2 cflux mg s-1, corrected from 0.15 diameter to 0.16
      respi = dfa['plant_respiration_CO2'].mean()
      WTs[n] = dfa['WT_cm'].mean()/100
      bd_as[n] = dfa['bds_a'].mean()
      bd_bs[n] = dfa['bds_b'].mean()
      bd_cs[n] = dfa['bds_c'].mean()
      delta_DOC[n] = dfa['TC'].values[-1] - dfa['TC'].values[0] 
      co2s = np.zeros(162)                                                          # output array, lenght in days
      for k in range(162):
        fco2 = interp1d(intervals, meas)   
        co2s[k] = fco2(k) * 86400 /  1000000 *  10000                                 #  CO2 kg/ha/day
      co2out[n]=sum(co2s) - respi
      #print (col, respi, co2out[n])
    print ('Finland, vege:', vege)
    print (columns, WTs, co2out)
    print ('xxxxxxxxxxxxxxxxxx')
    print ('WT', np.mean(WTs), 'WTsd',np.std(WTs), 'CO2', np.nanmean(co2out))
    return columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs

def otherlands(country, vege=True):
    intervals = np.array([  0,  38,  73, 109, 162])
    cols = ['Column_no', 'Month','Code', 'Country', 'Site',
           'Sample_ID', 'Treat', 'WT',
           'WT_cm_from_top_of_the_column', 'Site','bd_a','bd_b','bd_c','Peat_surface_from_top_of_tube_cm',
           'T_C_mean', 'CO_flux_mg_s1_m2_mean','Vegetation','Total_dry_g', 'Vascular_photosynth_dry_g','TotalMoss_dry_g',
           'nmonth', 'pvm', 'T_soil_C',
           'TC', 'SUVA254', 'E2E3']
    df = pd.read_excel('C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/egu2023/reformwater_data.xlsx', sheet_name='data', usecols =cols)
    df['pvm'][0].date()
    #samples = np.unique(df['Column_no'].values)
    
    #-------PRE guessed values-------------------
    water_content_low = 63.0/100                         # m3 m-3 potential -8 kPa -80 cm h2o Laurén et al. 2019
    water_content_high = 71.0/100                        # m3 m-3 potential -3 kPa -30 cm h2o
    
    if vege:
        dfFIN = df[(df['Country']==country)&(df['Vegetation']=='Yes')].copy()
    else: 
        dfFIN = df[(df['Country']==country)&(df['Vegetation']!='Yes')].copy()

    #dfFIN = df[(df['Country']==country)].copy()
    cyls = np.unique(dfFIN['Column_no'].values)
    dfFIN['sample_height'] = 50.0 - dfFIN['Peat_surface_from_top_of_tube_cm']
    dfFIN['WT_cm_from_bottom'] = 50.0 + dfFIN['WT_cm_from_top_of_the_column']       # water table up from the bottom cm
    dfFIN['WT_cm'] = dfFIN['sample_height']- dfFIN['WT_cm_from_bottom']             # water table cm down from the sample surface
    #dfFIN['peat_volume_m3'] = dfFIN['sample_height']/100.0 * np.pi * 0.1**2         # peat volume in m3 
    dfFIN['peat_volume_m3'] = dfFIN['sample_height']/100.0 * np.pi * 0.08**2         # peat volume in m3 
    dfFIN['water_content'] = water_content_low
    dfFIN['water_content'][dfFIN['WT']=='H'] = water_content_high
    dfFIN['water_m3'] = dfFIN['water_content'] * dfFIN['peat_volume_m3'] 
    dfFIN=dfFIN[dfFIN['WT_cm']>0]
    #dfFIN['DOC_in_water_mg'] = dfFIN['TC'] * dfFIN['water_m3'] * 1000.0
    #dfFIN['DOC_in_water_mg'].mean() 
    
    mean_a, mean_b, mean_c =  dfFIN[['bd_a','bd_b','bd_c']].mean().values
    
    bds_a = np.where(np.isfinite(dfFIN['bd_a'].values) , dfFIN['bd_a'].values, mean_a)
    bds_b = np.where(np.isfinite(dfFIN['bd_b'].values) , dfFIN['bd_b'].values, mean_b)
    bds_c = np.where(np.isfinite(dfFIN['bd_c'].values) , dfFIN['bd_c'].values, mean_c)
    dfFIN['bds_a'] = bds_a
    dfFIN['bds_b'] = bds_b
    dfFIN['bds_c'] = bds_c
    bds = dfFIN[['bds_a', 'bds_b','bds_c']].mean(axis = 1)
    dfFIN['bulk_density'] = bds                                                     #kg/m3
    dfFIN['peat_mass'] = dfFIN['peat_volume_m3'] * dfFIN['bulk_density']
    
    columns = np.unique(dfFIN['Column_no'].values)
    #dfFIN['plant_kg_ha-1'] = dfFIN['Total_dry_g']/(0.1**2*np.pi)*10000/1000 *1.5
    #dfFIN['plant_kg_ha-1'] = dfFIN['Vascular_photosynth_dry_g']/(0.1**2*np.pi)*10000/1000 *  5.0 \
    #    + dfFIN['TotalMoss_dry_g']/(0.1**2*np.pi)*10000/1000 
    dfFIN['plant_kg_ha-1'] = dfFIN['Vascular_photosynth_dry_g']/(0.1**2*np.pi)*10000/1000 *  5.0 \
        + dfFIN['TotalMoss_dry_g']/(0.08**2*np.pi)*10000/1000 
    
    
    dfFIN['plant_respiration_CO2'] = dfFIN['plant_kg_ha-1'] *0.5* 44/12 
    
    co2out = np.zeros(len(columns))
    WTs = np.zeros(len(columns))
    bd_as = np.zeros(len(columns))
    bd_bs = np.zeros(len(columns))
    bd_cs = np.zeros(len(columns))
    delta_DOC = np.zeros(len(columns))
    
    for n, col in enumerate(columns):
      dfa = dfFIN[dfFIN['Column_no']==col].copy()                                   
      dfa = dfa.sort_values('nmonth')                                               # sort data to chronological order   
      m_peat = dfa['peat_mass'].mean()*1000                                         # peat mass in g
      meas = dfa['CO_flux_mg_s1_m2_mean'].values*0.879                                    # measured CO2 cflux mg s-1
      respi = dfa['plant_respiration_CO2'].mean()
      WTs[n] = dfa['WT_cm'].mean()/100
      bd_as[n] = dfa['bds_a'].mean()
      bd_bs[n] = dfa['bds_b'].mean()
      bd_cs[n] = dfa['bds_c'].mean()
      delta_DOC[n] = dfa['TC'].values[-1] - dfa['TC'].values[0] 
      co2s = np.zeros(162)                                                          # output array, lenght in days
      for k in range(162):
        lenarr = min([len(meas),len(intervals)])
        fco2 = interp1d(intervals[:lenarr], meas[:lenarr], bounds_error=False, fill_value=(meas[0],meas[-1]))   
        co2s[k] = fco2(k) * 86400 /  1000000 *  10000                                 #  CO2 kg/ha/day
      co2out[n]=sum(co2s) - respi
      #print (col, respi, co2out[n])
    
    print (columns, WTs, co2out)
    #print (np.mean(WTs), np.nanmean(co2out))
    print (country, vege)
    print ('xxxxxxxxxxxxxxxxxx')
    print ('WT', np.mean(WTs), 'WTsd',np.std(WTs), 'CO2', np.nanmean(co2out))
    
    return columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs

#FIN(vege = True)
#otherlands('ES', vege = False)
#otherlands('IR', vege = False)
#otherlands('SWE', vege = False)
