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
folderName = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/mikko_niemi/'

wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/mikko_niemi/'

mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/mikko_niemi/',
              'dominant':{1: 'susi_stand_development_ix_6.xlsx'},
              'subdominant':{0:'susi_stand_development_ix_6.xlsx'},
              'under':{0:'susi_stand_development_ix_6.xlsx'}} 
"""
mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/krycklan/',
              'dominant':{1: 'susi_krycklan_input_lyr_0.xlsx'},
              'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
              'under':{0:'susi_motti_input_lyr_2.xlsx'}} 
"""
wdata='ks_weather.csv'

start_date = datetime.datetime(2004,1,1)
end_date=datetime.datetime(2013,12,31)
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25

sarkaSim = 40.                                                                  # strip width, ie distance between ditches, m
n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip

ageSim = {'dominant': 100.*np.ones(n),
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
#spara['canopylayers']['dominant'][int(n/2):] = 2                                                                        
#spara['canopylayers']['subdominant'][:int(n/2)] = 1                                                                        
spara['cutting_yr'] = 2001
spara['drain_age'] =  100.
mass_mor = 1.616*np.log(spara['drain_age'])-1.409     #Pitkänen et al. 2012 Forest Ecology and Management 284 (2012) 100–106

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
spara['ditch depth 20y west'] = [-0.3, -0.6, -0.9]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
spara['ditch depth 20y east'] = [-0.3, -0.6, -0.9]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
spara['scenario name'] =  ['D30','D60', 'D90']                                #kasvunlisaykset
#spara['enable_peatmiddle'] = False,
#spara['enable_peatbottom'] = False
#print (spara)

susi = Susi()
 
susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc)
    
          
#%%             
from netCDF4 import Dataset 
import numpy as np
ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/mikko_niemi/susi.nc'
scen = 0
ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode


##TODO: balances only for 1:-1, exclude dtic nodes

soil_co2 = np.mean(ncf['balance']['C']['soil_c_balance_co2eq'][scen, :, 1:-1])
wt = np.mean((ncf['strip']['dwtyr_growingseason'][scen, :, 1:-1]))
stand_litter = np.mean(ncf['balance']['C']['stand_litter_in'] [scen, :, 1:-1])
gv_litter = np.mean(ncf['balance']['C']['gv_litter_in'] [scen, :, 1:-1])
esom_out = np.mean( ncf['balance']['C']['co2c_release'][scen, :, 1:-1] )
gv_change = np.mean(ncf['balance']['C']['gv_change'] [scen, :, 1:-1])
stand_change = np.mean(ncf['balance']['C']['stand_change'] [scen, :, 1:-1])


if np.median(sfc)< 4:
    ojanen2019 = -1*(-115 + 12*wt*-100)
else:
    ojanen2019 = -1*(-259 + 6*wt*-100)
print ('*************************************')
print ('soil co2 balance kg CO2/ha/yr', soil_co2) 
print ('mean growing season WT m', wt)
print ('stand litter kg C / yr', stand_litter)

print ('gv litter kg C / yr', gv_litter)
print ('esom out kg C', esom_out)
print ('esom out g CO2 m-2 yr-1', esom_out/10*44/12)

print ('litter in - esom out kg C/ha/yr', ((stand_litter+ gv_litter) - esom_out)*44/12/10)
print ('Ojanen 2019', ojanen2019)
print ('******* In compoents *******')
print ('gv mass change kg C/ha/yr', gv_change)
#print (ncf['balance']['C']['gv_change'] [scen, :, :])
#print (ncf['balance']['C']['soil_c_balance_co2eq'][scen, :, :])
#print (ncf['strip']['dwtyr_growingseason'][scen, :, :])
print ('stand_change', stand_change)
ncf.close()
#%%
ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/mikko_niemi/susi.nc'
scen = 0
ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
wts = ncf['strip']['dwt'][:, :, :]      #scen, day, column
ncf.close()
import matplotlib.pylab as plt
#print (np.min(wts[0,:,1:-1], axis=1))
plt.plot(np.min(wts[0,:,1:-1], axis=1))
plt.plot(np.min(wts[1,:,1:-1], axis=1))
plt.plot(np.min(wts[2,:,1:-1], axis=1))

