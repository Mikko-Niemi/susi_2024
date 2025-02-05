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

stands = ['HelsinkiP22.xls','HelsinkiP32.xls','HelsinkiP42.xls',
          'JyvaskylaP22.xls', 'JyvaskylaP32.xls', 'JyvaskylaP42.xls',
          'KemijarviP22.xls', 'KemijarviP32.xls', 'KemijarviP42.xls']

weather_inputs = ['katila_weather.csv','katila_weather.csv','katila_weather.csv',
                  'jaakkoinsuo_weather.csv', 'jaakkoinsuo_weather.csv', 'jaakkoinsuo_weather.csv',
                  'koirasuo_weather.csv', 'koirasuo_weather.csv', 'koirasuo_weather.csv']
outnames = ['H22.nc', 'H32.nc', 'H42.nc',
            'J22.nc','J32.nc','J42.nc', 
            'K22.nc','K32.nc','K42.nc']
sfcs = [2,3,4,
        2,3,4,
        2,3,4]
for stand, weather_input,outname,sitetype in zip(stands, weather_inputs,outnames, sfcs):
    #***************** local call for SUSI*****************************************************
    folderName=r'C:/Users/laurenan/OneDrive - University of Helsinki/codes/susi_11/outputs/' #'sensitivity/'
    wpath = r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi/spinup/'
    
    
    
    mottifile = {'path': r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi/spinup/',
                  'dominant':{1: stand},
                  'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
                  'under':{0:'susi_motti_input_lyr_2.xlsx'}} 
    
    wdata=weather_input
    
    start_date = datetime.datetime(1961,1,1)
    end_date=datetime.datetime(2010,12,31)
    start_yr = start_date.year 
    end_yr = end_date.year
    yrs = (end_date - start_date).days/365.25
    
    sarkaSim = 40.                                                                  # strip width, ie distance between ditches, m
    n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip
    
    ageSim = {'dominant': 20.*np.ones(n),
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
    mass_mor = 0.5 # 1.616*np.log(spara['drain_age'])-1.409     #Pitkänen et al. 2012 Forest Ecology and Management 284 (2012) 100–106
    
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
    
    spara['ditch depth west'] = [-0.5]   #nLyrs kerrosten lkm, dzLyr kerroksen paksuus m, saran levys m, n laskentasolmulen lukumäärä, ditch depth pjan syvyys simuloinnin alussa m  
    spara['ditch depth east'] = [-0.5]
    spara['ditch depth 20y west'] = [-0.5]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
    spara['ditch depth 20y east'] = [-0.5]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
    spara['scenario name'] =  ['D50']                                #kasvunlisaykset
    #spara['enable_peatmiddle'] = False,
    #spara['enable_peatbottom'] = False
    #print (spara)
    outpara['netcdf'] = outname
    susi = Susi()
     
    susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                    mottifile=mottifile, peat= 'other', photosite='All data', 
                                    folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc)
    
    break    
#%%          
             
from netCDF4 import Dataset 
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

# @title Valitse muuttuja { run: "auto" }
muuttuja = 'kangashumus' #@param ['kangashumus', 'turve']
aine = 'K' #@param ['Mass', 'N', 'P', 'K']
scen = 0
outnames = ['H22.nc', 'H32.nc', 'H42.nc',
            'J22.nc','J32.nc','J42.nc', 
            'K22.nc','K32.nc','K42.nc']
folderName=r'C:/Users/laurenan/OneDrive - University of Helsinki/codes/susi_11/outputs/' #'sensitivity/'

for outn in outnames:
    ff = folderName + outn
    
    ncf=Dataset(ff, mode='r')
    peat = ncf['esom'][aine]['P1'][scen,:, :]/10000. + ncf['esom'][aine]['P2'][scen,:, :]/10000. + ncf['esom'][aine]['P3'][scen,:, :]/10000.
    dfpeat = pd.DataFrame(peat*10000)
    mor = ncf['esom'][aine]['LL'][scen,:, :]/10000. + ncf['esom'][aine]['LW'][scen,:, :]/10000. + ncf['esom'][aine]['FL'][scen,:, :]/10000.\
     + ncf['esom'][aine]['FW'][scen,:, :]/10000. + ncf['esom'][aine]['H'][scen,:, :]/10000.
    dfmor = pd.DataFrame(mor*10000)
    
    litterin = ncf['esom'][aine]['L0L'][scen, :, :] + ncf['esom'][aine]['L0W'][scen, :, :] 
    
    #release = ncf['balance'][aine]['decomposition_tot'][scen, :, :]    
    vars = {'kangashumus': dfmor, 'turve': dfpeat}
    
    #data_table.DataTable(vars[muuttuja], include_index=False, num_rows_per_page=25, max_columns = 50)
    ncf.close()
    print (np.mean(litterin, axis=1))
    print (np.mean(dfmor, axis=1))
    print (np.std(dfmor, axis = 1))
    plt.figure(outn)
    plt.plot(np.mean(dfmor, axis=1))
    plt.plot(np.cumsum(np.mean(litterin, axis=1)))
    plt.title(outn)
    #plt.ylim([0, 300000])

#%%          
             
from netCDF4 import Dataset 
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

# @title Valitse muuttuja { run: "auto" }
muuttuja = 'kangashumus' #@param ['kangashumus', 'turve']
aine = 'K' #@param ['Mass', 'N', 'P', 'K']
scen = 0
outnames = ['H22.nc', 'H32.nc', 'H42.nc',
            'J22.nc','J32.nc','J42.nc', 
            'K22.nc','K32.nc','K42.nc']
folderName=r'C:/Users/laurenan/OneDrive - University of Helsinki/codes/susi_11/outputs/' #'sensitivity/'

for outn in outnames:
    ff = folderName + outn
    
    ncf=Dataset(ff, mode='r')
    peat = ncf['esom'][aine]['P1'][scen,:, :]/10000. + ncf['esom'][aine]['P2'][scen,:, :]/10000. + ncf['esom'][aine]['P3'][scen,:, :]/10000.
    dfpeat = pd.DataFrame(peat*10000)
    mor = ncf['esom'][aine]['LL'][scen,:, :]/10000. + ncf['esom'][aine]['LW'][scen,:, :]/10000. + ncf['esom'][aine]['FL'][scen,:, :]/10000.\
     + ncf['esom'][aine]['FW'][scen,:, :]/10000. + ncf['esom'][aine]['H'][scen,:, :]/10000.
    dfmor = pd.DataFrame(mor*10000)
    
    litterin = ncf['esom'][aine]['L0L'][scen, :, :] + ncf['esom'][aine]['L0W'][scen, :, :] 
    
    release = ncf['balance'][aine]['decomposition_tot'][scen, 1:, :]    
    dfrelease = np.mean(pd.DataFrame(release), axis=1) 
    morrelease = (dfrelease + np.diff(np.mean(dfpeat, axis = 1)))
    morfraction = (morrelease / dfrelease)
    #data_table.DataTable(vars[muuttuja], include_index=False, num_rows_per_page=25, max_columns = 50)
    ncf.close()
    #print (dfrelease)
    
    plt.figure(outn)
    plt.plot(morfraction)
    #plt.plot(np.cumsum(np.mean(litterin, axis=1)))
    plt.title(outn)
    #plt.ylim([0, 300000])
