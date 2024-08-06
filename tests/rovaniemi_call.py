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


#***************** local call for SUSI*****************************************************
#folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/' #'sensitivity/'
folderName=r'C:/Users/alauren/OneDrive - University of Eastern Finland/Opetus/Climate smart peatland mgmt/ccf/'
wpath = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/susi_experim/weather/'

mottifile = {'path': r'C:/Users/alauren/OneDrive - University of Eastern Finland/Henkilöt/Stenberg/motti_Arille/',
              'dominant':{1: 'koira11_A.xls'},
              'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
              'under':{0:'susi_motti_input_lyr_2.xlsx'}} 

wdata='koirasuo_weather.csv'

start_date = datetime.datetime(2008,1,1)
end_date=datetime.datetime(2014,12,31)
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25

sarkaSim = 37.                                                                  # strip width, ie distance between ditches, m
n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip

ageSim = {'dominant': 76.*np.ones(n),
          #'dominant': 2.*np.ones(n),
          'subdominant': 0*np.ones(n),
          'under': 0*np.ones(n)}                                                         # age of the stand in each node

sfc =  np.ones(n, dtype=int)*3                                                                        # site fertility class

#ageSim['dominant'][int(n/2):] = 2.
#ageSim[4:-4] = 2.

site = 'develop_scens'

forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
            
wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          n=n)

spara['vonP top']=[4,4,5,6,7,7]
spara['depoN']= 4.0, 
spara['depoP']= 0.1,
spara['depoK']= 0.65

spara['ditch depth west']= [-0.3, -0.5, -0.7, -0.9]   #nLyrs kerrosten lkm, dzLyr kerroksen paksuus m, saran levys m, n laskentasolmulen lukumäärä, ditch depth pjan syvyys simuloinnin alussa m  
spara['ditch depth east']= [-0.3, -0.5, -0.7, -0.9]
spara['ditch depth 20y west']= [-0.3, -0.5, -0.7, -0.9]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
spara['ditch depth 20y east']= [-0.3, -0.5, -0.7, -0.9]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
spara['scenario name']= ['DNM30', 'DNM50', 'DNM70', 'DNM90'] #kasvunlisaykset
#spara['fertilization']['application year']=2008

susi = Susi()
 
susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc)
    
          
             
#%%
from susi.figures import *
ff = folderName + 'susi.nc' #r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/susi.nc'
scen = 0
hydrology(ff, scen)
stand(ff, scen)
mass(ff, scen)
carbon(ff, scen)
nutrient_balance(ff, 'N', scen)
nutrient_balance(ff, 'P', scen)
nutrient_balance(ff, 'K', scen)
compare_1(ff, [0,1])
compare_scens(ff)

#%%
from susi.figures import compare_runs
scen = 1
folderName=r'C:/Users/alauren/OneDrive - University of Eastern Finland/Opetus/Climate smart peatland mgmt/ccf/'
ff_0 = folderName + 'susi_nocut.nc'
ff_1 = folderName + 'susi_ba02.nc'
compare_runs(ff_0, ff_1, scen)
