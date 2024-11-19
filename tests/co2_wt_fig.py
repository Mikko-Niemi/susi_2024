# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 09:53:54 2024

@author: laurenan
"""

import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec

from scipy import stats
from netCDF4 import Dataset  

#from sklearn.metrics import r2_score
import seaborn as sns
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
sns.set()

def difference_plot(ax, x, y, cols, ylabel, label, fs, facecolor, colorin, hidex=False, hidey=False, elevation=None):
    
    return ax

facecolor = '#f2f5eb'
fs = 15
fig = plt.figure(num='hydro', figsize=(15,18))   #width, height
gs = gridspec.GridSpec(ncols=12, nrows=12, figure=fig, wspace=0.25, hspace=0.25)
ax = fig.add_subplot(gs[0:4, 0:4])


#%%
#REFERENCE data to Ojanen & Minkkinen 2019
fig = plt.figure(num='co2', figsize=(8,5))

fdata = r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi/CO2_vs_WT.xlsx'
ojanen_data = pd.read_excel(fdata,'CO2_vs_WT', 
                            usecols=['Site', 'N', 'E', 'Region', 'Site type', 'Trophic level', 'WT', 'CO2', 'Source'],
                            nrows=76)
sns.scatterplot(data=ojanen_data, x='WT', y='CO2', hue='Region', style='Site type', s = 150)
plt.ylabel('CO$_2$ g m$^{-2}$ yr$^{-1}$')
plt.xlabel('WT, cm')

# OWN simulations
scen = 0

outnames = ['SF_22.nc','SF_31.nc','SF_32.nc','SF_41.nc','SF_51.nc',
            'CF_22.nc','CF_31.nc','CF_32.nc','CF_41.nc','CF_51.nc',
            'NOBK_22.nc','NOBK_31.nc','NOBK_32.nc','NOBK_41.nc','NOBK_51.nc',
            'Lap_22.nc','Lap_31.nc','Lap_32.nc','Lap_41.nc','Lap_51.nc',
            ]

colors = ['blue','blue','blue','blue','blue',
           'orange', 'orange','orange','orange','orange',
           'green','green','green','green','green',
           'red', 'red','red','red','red' ]

mks =['o', 'X', 'X', 's', '+',
      'o', 'X', 'X', 's', '+',
      'o', 'X', 'X', 's', '+',
      'o', 'X', 'X', 's', '+']

#folderName=r'C:/Users/laurenan/OneDrive - University of Helsinki/codes/susi_11/outputs/' #'sensitivity/'
folderName=r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi//outputs/'

for outn, colos, mk in zip(outnames, colors, mks):
    ff = folderName + outn
    ncf=Dataset(ff, mode='r')
    
    bal = (ncf['balance']['C']['soil_c_balance_co2eq'][scen, 1:, 1:-1]
              + ncf['balance']['C']['LMWdoc_to_water'][scen, 1:, 1:-1]*44/12
              + ncf['balance']['C']['HMW_to_water'][scen, 1:, 1:-1]*44/12)/10*-1
    
    dfbal = pd.DataFrame(bal)
    wt =  ncf['strip']['dwtyr_growingseason'][scen, 1:, 1:-1]
    dfwt = pd.DataFrame(wt)
    ncf.close()
    #for c in range(17):
    #    wts = dfwt[c].values*-100    
    #    bals = dfbal[c].values    
    #    plt.plot(wts, bals, 'o', color = 'grey', alpha = 0.20)
    wts = dfwt.mean(axis=1).values*-100    
    bals = dfbal.mean(axis=1).values
    plt.plot(wts, bals, mk, color = colos, alpha = 0.5,markerfacecolor='none', markersize=4)
    plt.xlim([0,150])
    plt.ylim([-600, 2000])
    #print (dfwt)
    #dfwt.to_clipboard()
   

#%%
#CAI figure
fig = plt.figure(num='growth')

for outn, colos, mk in zip(outnames, colors, mks):
    ff = folderName + outn
    ncf=Dataset(ff, mode='r')
    
    growth0 = np.mean(ncf['stand']['volumegrowth'][0, 1:, 1:-1])
    growth1 = np.mean(ncf['stand']['volumegrowth'][1, 1:, 1:-1])
    growth2 = np.mean(ncf['stand']['volumegrowth'][2, 1:, 1:-1])
    
    wt0 =  np.mean(ncf['strip']['dwtyr_growingseason'][0, 1:, 1:-1])*-100
    wt1 =  np.mean(ncf['strip']['dwtyr_growingseason'][1, 1:, 1:-1])*-100
    wt2 =  np.mean(ncf['strip']['dwtyr_growingseason'][2, 1:, 1:-1])*-100
     
    wts = np.array([wt0, wt1,wt2]) 
    gr = np.array([growth0, growth1, growth2])
    plt.plot(wts,gr, mk, color = colos, alpha = 0.5,markerfacecolor='none', markersize=5, linestyle='-', label=outn[:-3])
    #plt.xlim([0,150])
    #plt.ylim([-600, 2000])
    plt.xlabel('WT, cm')
    plt.ylabel('CAI m$^3$ ha$^{-1}$ yr$^{-1}$')
    ncf.close()
#plt.legend(loc='best', ncols=5)

#plt.legend(['o', 'X', 's', '+'],['RhTkg', 'MTkg', 'PTkg', 'VaTkg'], loc='lower right')

    
#%%
# Volume change figure
fig = plt.figure(num='vol change')
t = np.arange(1,21,1)

for outn, colos, mk in zip(outnames[15:20], colors[15:20], mks[15:20]):
    ff = folderName + outn
    ffs = r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi//outputs_shallow/' + outn
 
    ncf=Dataset(ff, mode='r')
    ncfs=Dataset(ffs, mode = 'r') 
    
    growth = ncf['stand']['volumegrowth'][:, 1:, 1:-1]
    growths = ncfs['stand']['volumegrowth'][:, 1:, 1:-1]
    deltagr = growths - growth 
    wts = ncf['strip']['dwtyr_growingseason'][:, 1:, 1:-1]
    wtss = ncfs['strip']['dwtyr_growingseason'][:, 1:, 1:-1]
    deltawts = wtss - wts
    linestyles = ['dotted','dashed', 'solid' ]
    scens = [0, 1, 2]
    labs = ['0.3 m', '0.6 m', '0.9 m']
    for scen, ls, la in zip(scens, linestyles, labs):
        dg=np.cumsum(np.mean(deltagr[scen,:,:],axis=1))
        dw=np.mean(deltawts[scen,:,:],axis=1)
        plt.plot(t, dg, mk, color = colos, alpha = 0.5,markerfacecolor='none', markersize=5, linestyle=ls, label = la) 
    plt.xlabel('Time, years')
    plt.ylabel('$\Delta$ volume growth, m$^{3}$ ha$^{-1}$ ')
    #plt.legend(loc='lower left')
    plt.ylim([-30, 30])
    ncf.close()
    ncfs.close()
#%%
# Water table change
fig = plt.figure(num='wt change')
t = np.arange(1,21,1)
for outn, colos, mk in zip(outnames, colors, mks):
    ff = folderName + outn
    ffs = r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi//outputs_shallow/' + outn
 
    ncf=Dataset(ff, mode='r')
    ncfs=Dataset(ffs, mode = 'r') 
    
    wts = ncf['strip']['dwtyr_growingseason'][:, 1:, 1:-1]*100
    wtss = ncfs['strip']['dwtyr_growingseason'][:, 1:, 1:-1]*100
    deltawts = wtss - wts
 
    for scen in range(3):
        dw=np.mean(deltawts[scen,:,:],axis=1)
        plt.plot(t, dw, mk, color = colos, alpha = 0.5,markerfacecolor='none', markersize=5, linestyle='-') 
    plt.xlabel('Time, years')
    plt.ylabel('$\Delta$ WT, m')
    ncf.close()
    ncfs.close()
#%%
#Ecosystem C balance change
fig = plt.figure(num='Cbal change')
t = np.arange(1,21,1)
for outn, colos, mk in zip(outnames, colors, mks):
    ff = folderName + outn
    ffs = r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi//outputs_shallow/' + outn
 
    ncf=Dataset(ff, mode='r')
    ncfs=Dataset(ffs, mode = 'r') 
    
    standc = ncf['balance']['C']['stand_c_balance_c'][:, 1:, 1:-1]
    standcs = ncfs['balance']['C']['stand_c_balance_c'][:, 1:, 1:-1]
    deltastandc = standcs - standc 
   

    for scen in range(3):
        ds=np.cumsum(np.mean(deltastandc[scen,:,:],axis=1))
        plt.plot(t, ds, mk, color = colos, alpha = 0.5,markerfacecolor='none', markersize=5, linestyle='-') 
    plt.xlabel('Time, years')
    plt.ylabel('$\Delta$ ecosystem C balance, kg C ha$^{-1}$')
 
    ncf.close()
    ncfs.close()
#%%
# soil c balance change    
fig = plt.figure(num='soil Cbal change')
t = np.arange(1,21,1)
for outn, colos, mk in zip(outnames, colors, mks):
    ff = folderName + outn
    ffs = r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi//outputs_shallow/' + outn
 
    ncf=Dataset(ff, mode='r')
    ncfs=Dataset(ffs, mode = 'r') 
    
    soilc = ncf['balance']['C']['soil_c_balance_c'][:, 1:, 1:-1]
    soilcs = ncfs['balance']['C']['soil_c_balance_c'][:, 1:, 1:-1]
    deltasoilc = soilcs - soilc        

    for scen in range(3):
        ds=np.cumsum(np.mean(deltasoilc[scen,:,:],axis=1))
        plt.plot(t, ds, mk, color = colos, alpha = 0.5,markerfacecolor='none', markersize=5, linestyle='-') 
    plt.xlabel('Time, years')
    plt.ylabel('$\Delta$ soil C balance, kg C ha$^{-1}$')

    ncf.close()
    ncfs.close()
#%%
# N export change
fig = plt.figure(num='N export change')
t = np.arange(1,21,1)
for outn, colos, mk in zip(outnames, colors, mks):
    ff = folderName + outn
    ffs = r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi//outputs_shallow/' + outn
 
    ncf=Dataset(ff, mode='r')
    ncfs=Dataset(ffs, mode = 'r') 
    
    n = ncf['balance']['N']['to_water'][:, 1:, 1:-1]
    ns = ncfs['balance']['N']['to_water'][:, 1:, 1:-1]
    deltans = ns - n 
   

    for scen in range(3):
        ds=np.cumsum(np.mean(deltans[scen,:,:],axis=1))
        plt.plot(t, ds, mk, color = colos, alpha = 0.5,markerfacecolor='none', markersize=5, linestyle='-') 
    plt.xlabel('Time, years')
    plt.ylabel('$\Delta$ N export, kg ha$^{-1}$')

    ncf.close()
    ncfs.close()
#%%
# P export change
fig = plt.figure(num='P export change')
t = np.arange(1,21,1)
for outn, colos, mk in zip(outnames, colors, mks):
    ff = folderName + outn
    ffs = r'C:/Users/laurenan/OneDrive - University of Helsinki/SUSI/mikko_niemi//outputs_shallow/' + outn
 
    ncf=Dataset(ff, mode='r')
    ncfs=Dataset(ffs, mode = 'r') 
    
    n = ncf['balance']['P']['to_water'][:, 1:, 1:-1]
    ns = ncfs['balance']['P']['to_water'][:, 1:, 1:-1]
    deltans = ns - n 
   

    for scen in range(3):
        ds=np.cumsum(np.mean(deltans[scen,:,:],axis=1))
        plt.plot(t, ds, mk, color = colos, alpha = 0.5,markerfacecolor='none', markersize=5, linestyle='-') 
    plt.xlabel('Time, years')
    plt.ylabel('$\Delta$ P export, kg ha$^{-1}$')

    ncf.close()
    ncfs.close()