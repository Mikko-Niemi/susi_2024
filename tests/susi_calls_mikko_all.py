# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:10:42 2020

@author: alauren
"""
import numpy as np
import datetime
from susi.susi_utils import read_FMI_weather
from inputs.susi_para import get_susi_para
from susi.susi_main import Susi

stands = [
          'SF_22.xlsx','SF_31.xlsx','SF_32.xlsx', 'SF_41.xlsx','SF_51.xlsx',
          'CF_22.xlsx','CF_31.xlsx','CF_32.xlsx', 'CF_41.xlsx','CF_51.xlsx',
          'NOBK_22.xlsx','NOBK_31.xlsx','NOBK_32.xlsx', 'NOBK_41.xlsx','NOBK_51.xlsx',
          'Lap_22.xlsx','Lap_31.xlsx','Lap_32.xlsx', 'Lap_41.xlsx','Lap_51.xlsx'
          ]

weather_inputs = ['SFw.csv','SFw.csv','SFw.csv','SFw.csv','SFw.csv',
                  'CFw.csv','CFw.csv','CFw.csv','CFw.csv','CFw.csv',
                  'NOBKw.csv','NOBKw.csv','NOBKw.csv','NOBKw.csv','NOBKw.csv',
                  'Lapw.csv','Lapw.csv','Lapw.csv','Lapw.csv','Lapw.csv'
                  ]
outnames = ['SF_22.nc','SF_31.nc','SF_32.nc','SF_41.nc','SF_51.nc',
            'CF_22.nc','CF_31.nc','CF_32.nc','CF_41.nc','CF_51.nc',
            'NOBK_22.nc','NOBK_31.nc','NOBK_32.nc','NOBK_41.nc','NOBK_51.nc',
            'Lap_22.nc','Lap_31.nc','Lap_32.nc','Lap_41.nc','Lap_51.nc',
            ]
sfcs = [2,3,3,4,5,
        2,3,3,4,5,
        2,3,3,4,5,
        2,3,3,4,5
        ]

humus =[6.9, 5.8, 7.0, 5.5, 5.2,                 #20yrs spinoff
        6.7, 5.6, 6.7, 5.3, 4.9,
        7.0, 5.4, 7.0, 5.5, 4.9,
        6.9, 5.9, 7.2, 6.0, 5.2
        ]
for stand, weather_input,outname,sitetype,hum in zip(stands, weather_inputs,outnames, sfcs, humus):
    #***************** local call for SUSI*****************************************************
    folderName=r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi//outputs/' #_shallow/' 
    wpath = r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi/weather_data/'
    
    
    
    mottifile = {'path': r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi/motti_files/',
                  'dominant':{1: stand},
                  'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
                  'under':{0:'susi_motti_input_lyr_2.xlsx'}} 
    
    wdata=weather_input
    
    start_date = datetime.datetime(2004,1,1)
    end_date=datetime.datetime(2023,12,31)
    start_yr = start_date.year 
    end_yr = end_date.year
    yrs = (end_date - start_date).days/365.25
    
    sarkaSim = 40.                                                                  # strip width, ie distance between ditches, m
    n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip
    
    ageSim = {'dominant': 60.*np.ones(n),
              'subdominant': 0*np.ones(n),
              'under': 0*np.ones(n)}                                                         # age of the stand in each node
    
    sfc =  np.ones(n, dtype=int)*sitetype                                                                        # site fertility class
    
    #ageSim['dominant'][int(n/2):] = 2.
    #ageSim[4:-4] = 2.
    
    site = 'develop_scens'
    
    forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
                
    wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                              folderName=folderName, hdomSim=None,  
                                                                              ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                              n=n)
    #spara['canopylayers']['dominant'][int(n/2):] = 2                                                                        
    #spara['canopylayers']['subdominant'][:int(n/2)] = 1                                                                        
    #spara['cutting_yr'] = 3001
    spara['drain_age'] =  100.
    mass_mor = hum # 1.616*np.log(spara['drain_age'])-1.409     #Pitkänen et al. 2012 Forest Ecology and Management 284 (2012) 100–106
    
    if np.median(sfc) > 4:
        spara['peat type']=['S','S','S','S','S','S','S','S']
        spara['peat type bottom']=['A']
        spara['vonP top'] =  [2,5,5,5,6,6,7,7] 
        spara['anisotropy'] = 10
        spara['rho_mor'] = 80.0
    else:
        spara['vonP top'] =  [2,5,5,5,6,6,7,7] 
        spara['anisotropy'] = 10
        spara['rho_mor'] = 90.0
        
    
    spara['h_mor'] = mass_mor/ spara['rho_mor'] 
    
    spara['ditch depth west'] = [-0.3, -0.6, -0.9]   #nLyrs kerrosten lkm, dzLyr kerroksen paksuus m, saran levys m, n laskentasolmulen lukumäärä, ditch depth pjan syvyys simuloinnin alussa m  
    spara['ditch depth east'] = [-0.3, -0.6, -0.9]
    #spara['ditch depth 20y west'] = [-0.2, -0.4, -0.6]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
    #para['ditch depth 20y east'] = [-0.2, -0.4, -0.6]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
    spara['ditch depth 20y west'] = [-0.3, -0.6, -0.9]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
    spara['ditch depth 20y east'] = [-0.3, -0.6, -0.9]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
    
    spara['scenario name'] =  ['D30', 'D60','D90']                                #kasvunlisaykset
    
    #spara['enable_peatmiddle'] = False,
    #spara['enable_peatbottom'] = False
    #print (spara)
    outpara['netcdf'] = outname
    susi = Susi()
     
    susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                    mottifile=mottifile, peat= 'other', photosite='All data', 
                                    folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc)
    #if stand == 'SF_32.xlsx': break
        
#%%          
             
from netCDF4 import Dataset 
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

# @title Valitse muuttuja { run: "auto" }
muuttuja = 'kangashumus' #@param ['kangashumus', 'turve']
aine = 'Mass' #@param ['Mass', 'N', 'P', 'K']
scen = 0

outnames = ['SF_22.nc','SF_31.nc','SF_32.nc','SF_41.nc','SF_51.nc',
            'CF_22.nc','CF_31.nc','CF_32.nc','CF_41.nc','CF_51.nc',
            'NOBK_22.nc','NOBK_31.nc','NOBK_32.nc','NOBK_41.nc','NOBK_51.nc',
            'Lap_22.nc','Lap_31.nc','Lap_32.nc','Lap_41.nc','Lap_51.nc',
            ]

#folderName=r'C:/Users/laurenan/OneDrive - University of Helsinki/codes/susi_11/outputs/' #'sensitivity/'
folderName=r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi//outputs/'

for outn in outnames:
    ff = folderName + outn
    
    ncf=Dataset(ff, mode='r')
    peat = ncf['esom'][aine]['P1'][scen,:, :]/10000. + ncf['esom'][aine]['P2'][scen,:, :]/10000. + ncf['esom'][aine]['P3'][scen,:, :]/10000.
    dfpeat = pd.DataFrame(peat*10000)
    mor = ncf['esom'][aine]['LL'][scen,:, :]/10000. + ncf['esom'][aine]['LW'][scen,:, :]/10000. + ncf['esom'][aine]['FL'][scen,:, :]/10000.\
     + ncf['esom'][aine]['FW'][scen,:, :]/10000. + ncf['esom'][aine]['H'][scen,:, :]/10000.
    dfmor = pd.DataFrame(mor*10000)
    
    litterin = ncf['esom'][aine]['L0L'][scen, :, :] + ncf['esom'][aine]['L0W'][scen, :, :] 
    
    #vars = {'kangashumus': dfmor, 'turve': dfpeat}
    
    #data_table.DataTable(vars[muuttuja], include_index=False, num_rows_per_page=25, max_columns = 50)
    ncf.close()
    print(outn)
    #print (np.mean(litterin, axis=1))
    print (np.mean(dfmor, axis=1))
    print('*********************')
    #print (np.std(dfmor, axis = 1))
    plt.figure(outn)
    plt.plot(np.mean(dfmor, axis=1))
    plt.plot(np.cumsum(np.mean(litterin, axis=1)))
    plt.title(outn)
    plt.ylim([0, 150000])
    
