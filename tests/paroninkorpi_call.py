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
from netCDF4 import Dataset 

#***************** local call for SUSI*****************************************************
#folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/' #'sensitivity/'

wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/inputs/'
folderName = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/paroninkorpi_doc/'

mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/paroninkorpi_doc/',
              'dominant':{1: 'susi_lyr_0_paroni_30_0.xlsx'},
              'subdominant':{0:'susi_lyr_0_paroni_30_0.xlsx'},
              'under':{0:'susi_lyr_0_paroni_30_0.xlsx'}} 

wdata='parkano_weather.csv'

start_date = datetime.datetime(2005,1,1)
end_date=datetime.datetime(2014,12,31)
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25

sarkaSim = 70.                                                                  # strip width, ie distance between ditches, m
n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip

ageSim = {'dominant': 62.*np.ones(n),
          'subdominant': 0*np.ones(n),
          'under': 0*np.ones(n)}                                                         # age of the stand in each node

sfc =  np.ones(n, dtype=int)*2                                                                        # site fertility class

#ageSim['dominant'][int(n/2):] = 2.
#ageSim[4:-4] = 2.

site = 'develop_scens'

forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
            
wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          n=n)
    
spara['bd top'] = [0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.16, 0.16]
spara['bd bottom'] = [0.16]

outpara['netcdf'] = 'susi_paroninkorpi_70m.nc'

susi = Susi()
 
susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc)
    
          
             
#%%
from susi.figures import *
ff=r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/paroninkorpi_doc/susi_paroninkorpi_60m.nc'
scen = 1

#hydrology(ff, scen)

#stand(ff, scen)

#mass(ff, scen)
#carbon(ff, scen)
#nutrient_balance(ff, 'N', scen)
# nutrient_balance(ff, 'P', scen)
#nutrient_balance(ff, 'K', scen)
#compare_1(ff, [0,1])
compare_scens(ff)
#%%

folderName = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/paroninkorpi_doc/'

ncf=Dataset(folderName+'susi_paroninkorpi_80m.nc', mode='r')
for scen in [0,1,2]:

    roff = ncf['strip']['roff'][scen, :]*1000
    sflow = np.mean(ncf['strip']['surfacerunoff'][scen, :, :], axis=1)*1000
    totr = roff + sflow
    nofdays = len(roff)
    docexport = ncf['export']['hmwtoditch'][scen,:, :] + ncf['export']['lmwtoditch'][scen,:, :]
    
    
    import pandas as pd
    dates = pd.date_range(str(start_date.date()),periods=nofdays)
    dfrunoff = pd.DataFrame(data={'date':dates, 'runoff':totr})
    dfrunoff = dfrunoff.set_index('date')
    water = dfrunoff['runoff'].resample('y').sum()
    docs = np.mean(docexport, axis=1)[1:]
    conc = (docs*1e6 )/( water.values*1e4)
    print ('scenario',scen)
    print ('exp',np.mean(docs), np.std(docs))
    print ('conc', np.mean(conc), np.std(conc))
ncf.close()

#%%
folderName = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/paroninkorpi_doc/'

ncf=Dataset(folderName+'susi_paroninkorpi_60m.nc', mode='r')
hmwtoditch0 = np.mean(ncf['export']['hmwtoditch'][0,:, :] + ncf['export']['lmwtoditch'][0,:, :])
hmwtoditch1 = np.mean(ncf['export']['hmwtoditch'][1,:, :] + ncf['export']['lmwtoditch'][1,:, :])
hmwtoditch2 = np.mean(ncf['export']['hmwtoditch'][2,:, :] + ncf['export']['lmwtoditch'][2,:, :])
print (hmwtoditch0, hmwtoditch1,hmwtoditch2)
print (np.sum(ncf['strip']['roff'][0, :]))
print (np.sum(ncf['strip']['roff'][1, :]))
print (np.sum(ncf['strip']['roff'][2, :]))
print('-------------------')
print (np.sum(ncf['strip']['surfacerunoff'][0, :]))
print (np.sum(ncf['strip']['surfacerunoff'][1, :]))
print (np.sum(ncf['strip']['surfacerunoff'][2, :]))



#print ((hmwtoditch0*1e6)/(np.sum(ncf['strip']['roff'][0, :])*1000/10*1000000))
#print ((hmwtoditch1*1e6)/(np.sum(ncf['strip']['roff'][1, :])*1000/10*1000000))
#print ((hmwtoditch2*1e6)/(np.sum(ncf['strip']['roff'][2, :])*1000/10*1000000))
print ('---------------')
#print (ncf['export']['hmwtoditch'][2,:, :])
ncf.close()                                        # water netCDF, open in reading mode

ncf=Dataset(folderName+'susi_paroninkorpi_40m.nc', mode='r')
hmwtoditch0 = np.mean(ncf['export']['hmwtoditch'][0,:, :] + ncf['export']['lmwtoditch'][0,:, :])
hmwtoditch1 = np.mean(ncf['export']['hmwtoditch'][1,:, :] + ncf['export']['lmwtoditch'][1,:, :])
hmwtoditch2 = np.mean(ncf['export']['hmwtoditch'][2,:, :] + ncf['export']['lmwtoditch'][2,:, :])
print ('-----40m----------')
print ('exports',hmwtoditch0, hmwtoditch1,hmwtoditch2)
wt0 =  np.mean(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1 =  np.mean(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2 =  np.mean(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

wt0sd =  np.std(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1sd =  np.std(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2sd =  np.std(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

print ('wts', wt0,wt1,wt2)
print ('wtssd', wt0sd,wt1sd,wt2sd)

ncf.close()                                 

ncf=Dataset(folderName+'susi_paroninkorpi_60m.nc', mode='r')
hmwtoditch0 = np.mean(ncf['export']['hmwtoditch'][0,:, :] + ncf['export']['lmwtoditch'][0,:, :])
hmwtoditch1 = np.mean(ncf['export']['hmwtoditch'][1,:, :] + ncf['export']['lmwtoditch'][1,:, :])
hmwtoditch2 = np.mean(ncf['export']['hmwtoditch'][2,:, :] + ncf['export']['lmwtoditch'][2,:, :])
print ('-----60m----------')
print ('exports', hmwtoditch0, hmwtoditch1,hmwtoditch2)
wt0 =  np.mean(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1 =  np.mean(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2 =  np.mean(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

wt0sd =  np.std(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1sd =  np.std(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2sd =  np.std(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

print ('wts', wt0,wt1,wt2)
print ('wtssd', wt0sd,wt1sd,wt2sd)
ncf.close()                                 


ncf=Dataset(folderName+'susi_paroninkorpi_80m.nc', mode='r')
hmwtoditch0 = np.mean(ncf['export']['hmwtoditch'][0,:, :] + ncf['export']['lmwtoditch'][0,:, :])
hmwtoditch1 = np.mean(ncf['export']['hmwtoditch'][1,:, :] + ncf['export']['lmwtoditch'][1,:, :])
hmwtoditch2 = np.mean(ncf['export']['hmwtoditch'][2,:, :] + ncf['export']['lmwtoditch'][2,:, :])
print ('-----80m----------')
print ('exports',hmwtoditch0, hmwtoditch1,hmwtoditch2)
wt0 =  np.mean(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1 =  np.mean(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2 =  np.mean(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

wt0sd =  np.std(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1sd =  np.std(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2sd =  np.std(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

print ('wts', wt0,wt1,wt2)
print ('wtssd', wt0sd,wt1sd,wt2sd)






ncf.close()       
#%%

wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/inputs/'
wdata='parkano_weather.csv'

start_date = datetime.datetime(2005,1,1)
end_date=datetime.datetime(2014,12,31)

forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)    
print (np.mean(forc['T'].resample('Y').mean()))
print (np.std(forc['T'].resample('Y').mean()))

#print(forc['Prec'].resample('Y').sum())

print (np.mean(forc['Prec'].resample('Y').sum()))
print (np.std(forc['Prec'].resample('Y').sum()))
                          