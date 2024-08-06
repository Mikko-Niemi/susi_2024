# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 00:54:09 2023

@author: alauren
"""

import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from susi.temperature import PeatTemperature
#from docclass import DocModel
from susi.esom import Esom

from susi.susi_io import print_site_description
from susi.susi_utils import  rew_drylimit
from susi.susi_utils import  get_temp_sum
from susi.strip import StripHydrology, drain_depth_development
from reformwater_experiments import FIN, otherlands


class Susi():
    def __init(self):
        pass

    def run_lab(self, WT, bd0, bd1, bd2): 
        from scipy.interpolate import interp1d
        ins = np.array([0,38,35,36,53])                                  # number of days between the measurements 
        intervals = np.cumsum(ins)                                       # number of day in the experiment  
        temps = np.array([24.9, 28.7, 26.1, 26.2,21.3 ])                 # temperatures in the measured days
        ft = interp1d(intervals, temps)                                  # interpolation function between the days and temperatures
        days = np.arange(0,162,1)                                        # contains 162 days
        temperatures = ft(days)-0                                         # interpolation   
       
        def ojanen_2010(sfc, v, t_peat, wt, bd):
          B = 350. ; T5zero = -46.02; T5ref = 10.
          Rref = (0.0695+3.7*10**(-4)*v +5.4*10**(-4)*bd + 1.2*10**(-3)*wt)*24.    #unit g/m2/h CO2-> g/m2/day
          Rhet = [rref*np.exp(B*((1./(T5ref-T5zero))-(1./(t_peat-T5zero)))) for rref in Rref]          
          Rhet = np.array(Rhet).T
          return np.sum(Rhet, axis=0)*10.

        spara ={ 

            'sfc':[3], 'sfc_specification': 1,
            'nLyrs':10, 'dzLyr': 0.05, 'L': 1, 'n':1, 

            'initial h': -0.2, 'slope': 0.0, 
            'peat type':['A','A','A','A','A','A','A','A'], 
            'peat type bottom':['A'],'anisotropy':10.,
            'vonP': False,
            'vonP top':  [2,2,3,4,4,5,6,6], 
            'vonP bottom': 8,
            'bd top':[bd0, bd0, bd1, bd1, bd1, bd1], 'bd bottom': bd2,
            'enable_peattop': True, 'enable_peatmiddle': True,
            'enable_peatbottom': True,
            
                }
        n = spara['n']                                                             # number of columns along the strip        
        
        sfc = spara['sfc']
        esmass = Esom(spara, sfc, 162, substance='Mass')                       # initializing organic matter decomposition instace for mass
        #WT = -0.35
        days = 162
        #******** Soil and strip parameterization *************************
        stp = StripHydrology(spara)                                                # initialize soil hydrology model
        
        pt = PeatTemperature(spara,24.0)                              # initialize peat temperature model       
        peat_temperatures = pt.create_outarrays(1, days, spara['nLyrs'])    
        df_peat_temperatures = pd.DataFrame(peat_temperatures[0,0:days,:],index=pd.date_range(datetime.date(2022,4,24),periods=days)) # all peat temperatures 
        #print (np.shape(df_peat_temperatures))
        for d,ta in enumerate(temperatures):
            df_peat_temperatures.iloc[d,:] =ta
        dwts = np.ones((days,1))*(WT)
        dfwt = pd.DataFrame(dwts,index=pd.date_range(datetime.date(2022,4,24),periods=days))                              # daily water table data frame
        afp = np.ones((days,10))*stp.dwtToAfp(WT)
        dfafp = pd.DataFrame(afp,index=pd.date_range(datetime.date(2022,4,24),periods=days))                              # air filled porosity
        
        data = {'T':temperatures}
        weather = pd.DataFrame(data,index=pd.date_range(datetime.date(2022,4,24),periods=days))
        
        nonwoodylitter  = 0
        woodylitter = 0
    # **************  Biogeochemistry ***********************************
        esmass.run_yr(weather, df_peat_temperatures, dfwt, nonwoodylitter, woodylitter)
        massout =  (esmass.mass[0,0,10,:])                #  kg m-2 day-1
        #import matplotlib.pylab as plt
        #plt.plot(np.arange(0,162,1),massout*10000*0.5*44/12)
       
        #esmass.compose_export(stp)
        docshare = 0.05 #0.05
        lmwtohmwshare = 0.04
        mass_to_c = 0.5
        hmw_prod = (1-lmwtohmwshare)*esmass.out * 1/(1+docshare) * docshare * mass_to_c 
        hmw_out = hmw_prod*np.exp(-0.0004*162) 
        lmw_prod = esmass.out * 1/(1+docshare) * docshare * mass_to_c 
        lmw_out = lmw_prod* lmwtohmwshare* np.exp(-0.15*162) 
        print ('----------------------------')
        #print ('mean T',np.mean(temperatures))
        co2emiss =massout[-1]*10000*0.5*44/12
        #print ('CO2 emission kg/ha', co2emiss)
        jauhis = (71.40*-WT+ 23.15)*1000/365*162
        #print ('Jauhis',jauhis)
        #print ('bulk density', bd0*1000)
        paavo = ojanen_2010(3, 10, temperatures, np.array([WT])*-100, bd0*1000)
        #print ('Ojanen 2010', paavo)
        #print ('Water table',WT)
        #print ('hmw out kg/ha',hmw_out)
        #print ('hmw prod', hmw_prod)
        #print ('lmw out kg/ha',lmw_out)
        #print ('water_content m', stp.dwtToSto(WT))
        delta_DOC = np.ravel(((hmw_out+lmw_out)*1000000)/(stp.dwtToSto(WT)*10000*1000))
        #print ('DOC change mg/l', delta_DOC)
        #print ('************************************')
        return co2emiss, delta_DOC[0]

MAXBD = 0.22
#%%
columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs = FIN(vege=True)
co2s = []; delta_docs =[]
for c, WT, co2meas, doc, bd0,bd1,bd2 in zip(columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs ):        
    s= Susi()
    a,b=s.run_lab(-WT, min(MAXBD,bd0/1000), min(MAXBD, bd1/1000), min(MAXBD,bd2/1000))
    co2s.append(a)
    delta_docs.append(b)

fin_co2_meas_v1 = co2out.copy()
fin_co2_calc_v1 = np.array(co2s)
fin_doc_meas_v1 = delta_DOC.copy()
fin_doc_calc_v1 = np.array(delta_docs)

print ('**********FINLAND*****************')
print ('Measured CO2 {:.2f} Calculated CO2 {:.2f} '. format(np.mean(fin_co2_meas_v1), np.mean(fin_co2_calc_v1)))
print ('Measured DOC {:.2f} Calculated DOC {:.2f} '. format(np.mean(fin_doc_meas_v1), np.mean(fin_doc_calc_v1)))

#%%
columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs = FIN(vege=False)
co2s = []; delta_docs =[]
for c, WT, co2meas, doc, bd0,bd1,bd2 in zip(columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs ):        
    
    s= Susi()
    a,b=s.run_lab(-WT, min(MAXBD,bd0/1000), min(MAXBD, bd1/1000), min(MAXBD,bd2/1000))
    co2s.append(a)
    delta_docs.append(b)

fin_co2_meas_v0 = co2out.copy()
fin_co2_calc_v0 = np.array(co2s)
fin_doc_meas_v0 = delta_DOC.copy()
fin_doc_calc_v0 = np.array(delta_docs)
print ('**********FINLAND*****************')
print ('Measured CO2 {:.2f} Calculated CO2 {:.2f} '. format(np.mean(fin_co2_meas_v0), np.mean(fin_co2_calc_v0)))
print ('Measured DOC {:.2f} Calculated DOC {:.2f} '. format(np.mean(fin_doc_meas_v0), np.mean(fin_doc_calc_v0)))

#%%
#print (fin_co2_meas, fin_co2_calc)
#print (len(fin_co2_meas), len(fin_co2_calc))
#%%
#s= Susi()
#bd =320 
#a,b=s.run_lab(-0.2, bd/1000, bd/1000, bd/1000)

#%%
columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs = otherlands('ES', vege=True)
co2s = []; delta_docs =[]
for c, WT, co2meas, doc, bd0,bd1,bd2 in zip(columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs ):        
    #print ('column',c, WT, bd0/1000, bd1/1000, bd2/1000)
    
    s= Susi()
    a,b=s.run_lab(-WT, min(MAXBD,bd0/1000), min(MAXBD, bd1/1000), min(MAXBD,bd2/1000))
    co2s.append(a)
    delta_docs.append(b)

es_co2_meas_v1 = co2out.copy()
es_co2_calc_v1 = np.array(co2s)
es_doc_meas_v1 = delta_DOC.copy()
es_doc_calc_v1 = np.array(delta_docs)

print ('**********ESTONIA*****************')
print ('Measured CO2 {:.2f} Calculated CO2 {:.2f} '. format(np.mean(es_co2_meas_v1), np.mean(es_co2_calc_v1)))
print ('Measured DOC {:.2f} Calculated DOC increase {:.2f}'. format(np.mean(es_doc_meas_v1), np.mean(es_doc_calc_v1)))
print (delta_DOC)
#%%
columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs = otherlands('ES', vege=False)
co2s = []; delta_docs =[]
for c, WT, co2meas, doc, bd0,bd1,bd2 in zip(columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs ):        
    #print ('column',c, WT, bd0/1000, bd1/1000, bd2/1000)
    
    s= Susi()
    a,b=s.run_lab(-WT, min(MAXBD,bd0/1000), min(MAXBD, bd1/1000), min(MAXBD,bd2/1000))
    co2s.append(a)
    delta_docs.append(b)

es_co2_meas_v0 = co2out.copy()
es_co2_calc_v0 = np.array(co2s)
es_doc_meas_v0 = delta_DOC.copy()
es_doc_calc_v0 = np.array(delta_docs)
print ('**********ESTONIA*****************')
print ('Measured CO2 {:.2f} Calculated CO2 {:.2f} '. format(np.mean(es_co2_meas_v0), np.mean(es_co2_calc_v0)))
print ('Measured DOC {:.2f} Calculated DOC increase {:.2f}'. format(np.mean(es_doc_meas_v0), np.mean(es_doc_calc_v0)))

#%%
columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs = otherlands('IR', vege=True)
co2s = []; delta_docs =[]
for c, WT, co2meas, doc, bd0,bd1,bd2 in zip(columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs ):        
    #print ('column',c, WT, bd0/1000, bd1/1000, bd2/1000)
    
    s= Susi()
    a,b=s.run_lab(-WT,min(MAXBD,bd0/1000), min(MAXBD, bd1/1000), min(MAXBD,bd2/1000))
    co2s.append(a)
    delta_docs.append(b)

ir_co2_meas_v1 = co2out.copy()
ir_co2_calc_v1 = np.array(co2s)
ir_doc_meas_v1 = delta_DOC.copy()
ir_doc_calc_v1 = np.array(delta_docs)
print ('**********IRELAND*****************')
print ('Measured CO2 {:.2f} Calculated CO2 {:.2f} '. format(np.mean(ir_co2_meas_v1), np.mean(ir_co2_calc_v1)))
print ('Measured DOC {:.2f} Calculated DOC increase {:.2f}'. format(np.mean(ir_doc_meas_v1), np.mean(ir_doc_calc_v1)))
#%%
columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs = otherlands('IR', vege=False)
co2s = []; delta_docs =[]

for c, WT, co2meas, doc, bd0,bd1,bd2 in zip(columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs ):        
    #print ('column',c, WT, bd0/1000, bd1/1000, bd2/1000)
    
    s= Susi()
    a,b=s.run_lab(-WT, min(MAXBD,bd0/1000), min(MAXBD, bd1/1000), min(MAXBD,bd2/1000))
    co2s.append(a)
    delta_docs.append(b)

ir_co2_meas_v0 = co2out.copy()
ir_co2_calc_v0 = np.array(co2s)
ir_doc_meas_v0 = delta_DOC.copy()
ir_doc_calc_v0 = np.array(delta_docs)
print ('**********IRELAND*****************')
print ('Measured CO2 {:.2f} Calculated CO2 {:.2f} '. format(np.mean(ir_co2_meas_v0), np.mean(ir_co2_calc_v0)))
print ('Measured DOC {:.2f} Calculated DOC increase {:.2f}'. format(np.mean(ir_doc_meas_v0), np.mean(ir_doc_calc_v0)))

#%%

columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs = otherlands('SWE', vege=True)
co2s = []; delta_docs =[]
for c, WT, co2meas, doc, bd0,bd1,bd2 in zip(columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs ):        
    #print ('column',c, WT, bd0/1000, bd1/1000, bd2/1000)
    
    s= Susi()
    a,b=s.run_lab(-WT, min(MAXBD,bd0/1000), min(MAXBD, bd1/1000), min(MAXBD,bd2/1000))
    co2s.append(a)
    delta_docs.append(b)

swe_co2_meas_v1 = co2out.copy()
swe_co2_calc_v1 = np.array(co2s)
swe_doc_meas_v1 = delta_DOC.copy()
swe_doc_calc_v1 = np.array(delta_docs)

print ('**********SWEDEN*****************')
print ('Measured CO2 {:.2f} Calculated CO2 {:.2f} '. format(np.mean(swe_co2_meas_v1), np.mean(swe_co2_calc_v1)))
print ('Measured DOC {:.2f} Calculated DOC increase {:.2f}'. format(np.mean(swe_doc_meas_v1), np.mean(swe_doc_calc_v1)))

#%%

columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs = otherlands('SWE', vege=False)
co2s = []; delta_docs =[]
for c, WT, co2meas, doc, bd0,bd1,bd2 in zip(columns, WTs, co2out, delta_DOC, bd_as, bd_bs, bd_cs ):        
    #print ('column',c, WT, bd0/1000, bd1/1000, bd2/1000)
    
    s= Susi()
    a,b=s.run_lab(-WT,min(MAXBD,bd0/1000), min(MAXBD, bd1/1000), min(MAXBD,bd2/1000))
    co2s.append(a)
    delta_docs.append(b)

swe_co2_meas_v0 = co2out.copy()
swe_co2_calc_v0 = np.array(co2s)
swe_doc_meas_v0 = delta_DOC.copy()
swe_doc_calc_v0 = np.array(delta_docs)

print ('**********SWEDEN*****************')
print ('Measured CO2 {:.2f} Calculated CO2 {:.2f} '. format(np.mean(swe_co2_meas_v0), np.mean(swe_co2_calc_v0)))
print ('Calculated DOC increase {:.2f}'.format(np.mean(swe_doc_calc_v0)))
print ('Measured DOC {:.2f} Calculated DOC increase {:.2f}'. format(np.mean(swe_doc_meas_v0), np.mean(swe_doc_calc_v0)))


#%%
#plt.plot(co2out,np.array(co2s), 'go')

data ={'measured': fin_co2_meas_v0, 'calculated':fin_co2_calc_v0}
dfout = pd.DataFrame(data=data)
fig = plt.figure(num='compare_scens', figsize=(12,6)) 
gs = gridspec.GridSpec(ncols=4, nrows=2, figure=fig, wspace=0.35, hspace=0.35)

ax = fig.add_subplot(gs[0,0])
colorin='orange'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
ax.set_ylabel('$CO_2$ emission $kg \ ha^{-1} \ yr^{-1} $')
ax.set_ylim([0,80000])
plt.title('FINLAND, No vegetation')
ax.text(0.03, 1, 'a', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)


ax = fig.add_subplot(gs[0,1])
data ={'measured': fin_co2_meas_v1, 'calculated':fin_co2_calc_v1}
dfout = pd.DataFrame(data=data)

colorin='green'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('CO2 emission kg/ha/yr')
ax.set_ylim([0,80000])
plt.title('FINLAND, Vegetation')
ax.text(0.03, 1, 'b', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)

#-------------------------------------------------------------


data ={'measured': es_co2_meas_v0, 'calculated':es_co2_calc_v0}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[0,2])
colorin='orange'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('CO2 emission kg/ha/yr')
ax.set_ylim([0,80000])
plt.title('ESTONIA, No vegetation')
ax.text(0.03, 1, 'c', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)


data ={'measured': es_co2_meas_v1, 'calculated':es_co2_calc_v1}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[0,3])
colorin='green'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('CO2 emission kg/ha/yr')
ax.set_ylim([0,80000])
plt.title('ESTONIA, Vegetation')
ax.text(0.03, 1, 'd', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)


#----------------------------------------------------------------------------

data ={'measured': ir_co2_meas_v0, 'calculated':ir_co2_calc_v0}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[1,0])
colorin='orange'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
ax.set_ylabel('$CO_2$ emission $kg \ ha^{-1} \ yr^{-1} $')
ax.set_ylim([0,80000])
plt.title('IRELAND, No Vegetation')
ax.text(0.03, 1, 'e', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)


data ={'measured': ir_co2_meas_v1, 'calculated':ir_co2_calc_v1}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[1,1])
colorin='green'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('co2 emission kg/ha/yr')
ax.set_ylim([0,80000])
plt.title('IRELAND, Vegetation')
ax.text(0.03, 1, 'f', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)


#---------------------------------------------------------------------------------


data ={'measured': swe_co2_meas_v0, 'calculated':swe_co2_calc_v0}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[1,2])
colorin='orange'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('co2 emission kg/ha/yr')
ax.set_ylim([0,80000])
plt.title('SWEDEN, No vegetation')
ax.text(0.03, 1, 'g', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)



data ={'measured': swe_co2_meas_v1, 'calculated':swe_co2_calc_v1}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[1,3])
colorin='green'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('co2 emission kg/ha/yr')
ax.set_ylim([0,80000])
plt.title('SWEDEN, Vegetation')
ax.text(0.03, 1, 'h', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)




#%%
#%%
#----------DOC----------------------------
fig = plt.figure(num='compare_doc', figsize=(12,6)) 
#gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, wspace=0.35, hspace=0.35)
gs = gridspec.GridSpec(ncols=4, nrows=2, figure=fig, wspace=0.35, hspace=0.35)

#------------FIN--------------------------
data ={'measured': fin_doc_meas_v0, 'calculated':fin_doc_calc_v0}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[0,0])
colorin='orange'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
ax.set_ylabel('$\Delta$ DOC $mg \ l^{-1}$')
ax.set_ylim([0,400])
plt.title('FINLAND, No Vegetation')
ax.text(0.03, 1, 'i', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)


data ={'measured': fin_doc_meas_v1, 'calculated':fin_doc_calc_v1}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[0,1])
colorin='green'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('delta DOC mg/l')
ax.set_ylim([0,400])
plt.title('FINLAND, Vegetation')
ax.text(0.03, 1, 'j', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)


#------------ES--------------------------
data ={'measured': es_doc_meas_v0, 'calculated':es_doc_calc_v0}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[0,2])
colorin='orange'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('$\Delta$ DOC $mg \ l^{-1}$')
ax.set_ylim([0,400])
plt.title('ESTONIA, No Vegetation')
ax.text(0.03, 1, 'k', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)

data ={'measured': es_doc_meas_v1, 'calculated':es_doc_calc_v1}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[0,3])
colorin='green'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('delta DOC mg/l')
ax.set_ylim([0,400])
plt.title('ESTONIA, Vegetation')
ax.text(0.03, 1, 'l', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)

#------------IR--------------------------
data ={'measured': ir_doc_meas_v0, 'calculated':ir_doc_calc_v0}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[1,0])
colorin='orange'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
ax.set_ylabel('$\Delta$ DOC $mg \ l^{-1}$')
ax.set_ylim([0,400])
plt.title('IRELAND, No Vegetation')
ax.text(0.03, 1, 'm', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)

data ={'measured': ir_doc_meas_v1, 'calculated':ir_doc_calc_v1}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[1,1])
colorin='green'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('delta DOC mg/l')
ax.set_ylim([0,400])
plt.title('IRELAND, Vegetation')
ax.text(0.03, 1, 'n', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)


#------------SWE--------------------------
data ={'measured': swe_doc_meas_v0, 'calculated':swe_doc_calc_v0}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[1,2])
colorin='orange'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('$\Delta$ DOC $mg \ l^{-1}$')
ax.set_ylim([0,400])
plt.title('SWEDEN, No Vegetation')
ax.text(0.03, 1, 'o', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)

data ={'measured': swe_doc_meas_v1, 'calculated':swe_doc_calc_v1}
dfout = pd.DataFrame(data=data)

ax = fig.add_subplot(gs[1,3])
colorin='green'
ax = dfout.boxplot(ax = ax,
              color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
              boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              flierprops=dict(linestyle='-', linewidth=1.5),
              medianprops=dict(linestyle='-', linewidth=1.5),
              whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
              showfliers=False, grid=False, rot=0)
#ax.set_ylabel('delta DOC mg/l')
ax.set_ylim([0,400])
plt.title('SWEDEN, Vegetation')
ax.text(0.03, 1, 'p', horizontalalignment='left', 
            verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)
#%% RMSE obs-pred figures
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec

# variables
fin_co2_meas_v1
fin_co2_calc_v1

fin_co2_meas_v0
fin_co2_calc_v0

es_co2_meas_v1
es_co2_calc_v1

es_co2_meas_v0
es_co2_calc_v0

ir_co2_meas_v1
ir_co2_calc_v1

ir_co2_meas_v0
ir_co2_calc_v0

swe_co2_meas_v1 
swe_co2_calc_v1 

swe_co2_meas_v0 
swe_co2_calc_v0 



fin_doc_meas_v1
fin_doc_calc_v1

fin_doc_meas_v0
fin_doc_calc_v0

es_doc_meas_v1
es_doc_calc_v1

es_doc_meas_v0
es_doc_calc_v0

ir_doc_meas_v1
ir_doc_calc_v1

ir_doc_meas_v0
ir_doc_calc_v0

swe_doc_meas_v1 
swe_doc_calc_v1 

swe_doc_meas_v0 
swe_doc_calc_v0 


fs=16
nsites = 11
figsi = (10,5)
fig= plt.figure(num = 'growth', figsize=figsi)   #Figsize(w,h), tuple inches 
col1 = 'yellow'
col2 = 'grey'
gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, wspace=0.25, hspace=0.25)

ax0 = fig.add_subplot(gs[0,0])
mval=60000.
plt.fill_between([0., mval],[mval,mval],[0.,mval], color=col1, alpha=0.3)
plt.fill_between([0., mval],[0.,mval], [0.,0.], color=col2, alpha=0.3)

plt.xlabel('Observed $\it{CO_2}$, $kg ha^{-1} yr^{-1}$', fontsize =fs )
plt.ylabel('Predicted $\it{CO_2}$, $kg ha^{-1} yr^{-1}$', fontsize =fs)

p = 1.0 #multiply std

obs = np.mean(fin_co2_meas_v1)
obserr=np.std(fin_co2_meas_v1)
pre = np.mean(fin_co2_calc_v1)
preerr=np.std(fin_co2_calc_v1)
plt.plot(obs,pre,'o', markersize=10, marker='o', color='blue', label='Finland, vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='blue', capsize=4)

obs = np.mean(fin_co2_meas_v0)
obserr=np.std(fin_co2_meas_v0)
pre = np.mean(fin_co2_calc_v0)
preerr=np.std(fin_co2_calc_v0)
plt.plot(obs,pre,'o', markersize=10, marker='o', color='blue', mfc='none', label = 'Finland, no vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='blue', capsize=4)

#-----------------------------------
obs = np.mean(es_co2_meas_v1)
obserr=np.std(es_co2_meas_v1)
pre = np.mean(es_co2_calc_v1)
preerr=np.std(es_co2_calc_v1)
plt.plot(obs,pre,'o', markersize=10, marker="D", color='grey',label='Estonia, vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='grey', capsize=4)

obs = np.mean(es_co2_meas_v0)
obserr=np.std(es_co2_meas_v0)
pre = np.mean(es_co2_calc_v0)
preerr=np.std(es_co2_calc_v0)
plt.plot(obs,pre,'o', markersize=10, marker="D", color='grey', mfc='none', label='Estonia, no vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='grey', capsize=4)
#------------------------------------------
obs = np.mean(ir_co2_meas_v1)
obserr=np.std(ir_co2_meas_v1)
pre = np.mean(ir_co2_calc_v1)
preerr=np.std(ir_co2_calc_v1)
plt.plot(obs,pre,'o', markersize=10, marker='v', color='green', label='Ireland, vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='green', capsize=4)

obs = np.mean(ir_co2_meas_v0)
obserr=np.std(ir_co2_meas_v0)
pre = np.mean(ir_co2_calc_v0)
preerr=np.std(ir_co2_calc_v0)
plt.plot(obs,pre,'o', markersize=10, marker='v', color='green', mfc='none',label='Ireland, no vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='green', capsize=4)

#------------------------------
obs = np.mean(swe_co2_meas_v1)
obserr=np.std(swe_co2_meas_v1)
pre = np.mean(swe_co2_calc_v1)
preerr=np.std(swe_co2_calc_v1)
plt.plot(obs,pre,'o', markersize=10, marker='s', color='orange', label='Sweden, vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='orange', capsize=4)

obs = np.mean(swe_co2_meas_v0)
obserr=np.std(swe_co2_meas_v0)
pre = np.mean(swe_co2_calc_v0)
preerr=np.std(swe_co2_calc_v0)
plt.plot(obs,pre,'o', markersize=10, marker='s', color='orange', mfc='none', label='Sweden, vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='orange', capsize=4)


plt.legend(loc='upper right')
for tick in ax0.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax0.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

observed_co2 = np.array([np.mean(fin_co2_meas_v1), np.mean(fin_co2_meas_v0),
                         np.mean(es_co2_meas_v1), np.mean(es_co2_meas_v0),
                         np.mean(ir_co2_meas_v1), np.mean(ir_co2_meas_v0),
                         np.mean(swe_co2_meas_v1), np.mean(swe_co2_meas_v0)])

calculated_co2 = np.array([np.mean(fin_co2_calc_v1), np.mean(fin_co2_calc_v0),
                         np.mean(es_co2_calc_v1), np.mean(es_co2_calc_v0),
                         np.mean(ir_co2_calc_v1), np.mean(ir_co2_calc_v0),
                         np.mean(swe_co2_calc_v1), np.mean(swe_co2_calc_v0)])

diff = observed_co2 - calculated_co2
diff = diff[~np.isnan(diff)]
rmse = np.round(np.sqrt(np.square(diff).mean()),0)
print (rmse)

a, _, _, _ = np.linalg.lstsq(observed_co2[:,np.newaxis],calculated_co2)
#eq = 'y = ' +  str(np.round(a[0],3)) + 'x'
#x =np.arange(0.0 , 50000, 1000)
#plt.plot(x, a*x, 'k--', linewidth=2)
#plt.text(50000, 30000, eq, fontsize = fs-1)
plt.text(40000, 35000, 'RMSE ' + str(rmse), fontsize=fs-1)
ax0.text(0.04, 0.95, 'a', horizontalalignment='left',
               verticalalignment='top', fontsize=18, transform = ax0.transAxes,  fontweight='bold')

plt.xlim([0,mval])
plt.ylim([0,mval])


#---------DOC-------------------------

ax1 = fig.add_subplot(gs[0,1])
mval=300.
mval2 = -100

plt.fill_between([mval2, mval],[mval,mval],[mval2,mval], color=col1, alpha=0.3)
plt.fill_between([mval2, mval],[mval2,mval], [mval2,mval2], color=col2, alpha=0.3)


plt.xlabel('Observed $\Delta DOC$, $mg L^{-1}$', fontsize =fs )
plt.ylabel('Predicted $\Delta DOC$, $mg L^{-1} $', fontsize =fs)


obs = np.mean(fin_doc_meas_v1)
obserr=np.std(fin_doc_meas_v1)
pre = np.mean(fin_doc_calc_v1)
preerr=np.std(fin_doc_calc_v1)
plt.plot(obs,pre,'o', markersize=10, marker='o', color='blue', label='Finland, vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='blue', capsize=4)

obs = np.mean(fin_doc_meas_v0)
obserr=np.std(fin_doc_meas_v0)
pre = np.mean(fin_doc_calc_v0)
preerr=np.std(fin_doc_calc_v0)
plt.plot(obs,pre,'o', markersize=10, marker='o', color='blue', mfc='none', label = 'Finland, no vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='blue', capsize=4)

#-----------------------------------
obs = np.mean(es_doc_meas_v1)
obserr=np.std(es_doc_meas_v1)
pre = np.mean(es_doc_calc_v1)
preerr=np.std(es_doc_calc_v1)
plt.plot(obs,pre,'o', markersize=10, marker="D", color='grey',label='Estonia, vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='grey', capsize=4)

obs = np.mean(es_doc_meas_v0)
obserr=np.std(es_doc_meas_v0)
pre = np.mean(es_doc_calc_v0)
preerr=np.std(es_doc_calc_v0)
plt.plot(obs,pre,'o', markersize=10, marker="D", color='grey', mfc='none', label='Estonia, no vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='grey', capsize=4)
#------------------------------------------
obs = np.mean(ir_doc_meas_v1)
obserr=np.std(ir_doc_meas_v1)
pre = np.mean(ir_doc_calc_v1)
preerr=np.std(ir_doc_calc_v1)
plt.plot(obs,pre,'o', markersize=10, marker='v', color='green', label='Ireland, vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='green', capsize=4)

obs = np.mean(ir_doc_meas_v0)
obserr=np.std(ir_doc_meas_v0)
pre = np.mean(ir_doc_calc_v0)
preerr=np.std(ir_doc_calc_v0)
plt.plot(obs,pre,'o', markersize=10, marker='v', color='green', mfc='none',label='Ireland, no vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='green', capsize=4)

#------------------------------
obs = np.mean(swe_doc_meas_v1)
obserr=np.std(swe_doc_meas_v1)
pre = np.mean(swe_doc_calc_v1)
preerr=np.std(swe_doc_calc_v1)
plt.plot(obs,pre,'o', markersize=10, marker='s', color='orange', label='Sweden, vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='orange', capsize=4)

obs = np.mean(swe_doc_meas_v0)
obserr=np.std(swe_doc_meas_v0)
pre = np.mean(swe_doc_calc_v0)
preerr=np.std(swe_doc_calc_v0)
plt.plot(obs,pre,'o', markersize=10, marker='s', color='orange', mfc='none', label='Sweden, no vegetation')
plt.errorbar(obs,pre,preerr*p,obserr*p, 'none', color='orange', capsize=4)


plt.legend(loc='upper right')
for tick in ax0.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax0.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

observed_doc = np.array([np.mean(fin_doc_meas_v1), np.mean(fin_doc_meas_v0),
                         np.mean(es_doc_meas_v1), np.mean(es_doc_meas_v0),
                         np.mean(ir_doc_meas_v1), np.mean(ir_doc_meas_v0),
                         np.mean(swe_doc_meas_v1), np.mean(swe_doc_meas_v0)])

calculated_doc = np.array([np.mean(fin_doc_calc_v1), np.mean(fin_doc_calc_v0),
                         np.mean(es_doc_calc_v1), np.mean(es_doc_calc_v0),
                         np.mean(ir_doc_calc_v1), np.mean(ir_doc_calc_v0),
                         np.mean(swe_doc_calc_v1), np.mean(swe_doc_calc_v0)])

diff = observed_doc - calculated_doc
diff = diff[~np.isnan(diff)]
rmse = np.round(np.sqrt(np.square(diff).mean()),0)
print (rmse)

#print (np.mean(observed_co2, np.mean(calculated_co2))
       
print (np.mean(observed_doc), np.mean(calculated_doc))

a, _, _, _ = np.linalg.lstsq(observed_doc[:,np.newaxis],calculated_doc)
eq = 'y = ' +  str(np.round(a[0],3)) + 'x'
x =np.arange(0.0 , 300, 10)
#plt.plot(x, a*x, 'k--', linewidth=2)
#plt.text(200, 200, eq, fontsize = fs-1)
plt.text(150, 10, 'RMSE ' + str(rmse), fontsize=fs-1)
plt.xlim([mval2,mval])
plt.ylim([mval2,mval])
ax1.text(0.04, 0.95, 'b', horizontalalignment='left',
               verticalalignment='top', fontsize=18, transform = ax1.transAxes,  fontweight='bold')


#%%

fi_observed_co2 = np.array([np.mean(fin_co2_meas_v1), np.mean(fin_co2_meas_v1),
                         ])

fi_calculated_co2 = np.array([np.mean(fin_co2_calc_v1), np.mean(fin_co2_calc_v0),
                         ])

es_observed_co2 = np.array([
                         np.mean(es_co2_meas_v1), np.mean(es_co2_meas_v0),
                         ])

es_calculated_co2 = np.array([
                         np.mean(es_co2_calc_v1), np.mean(es_co2_calc_v0),
                         ])

ir_observed_co2 = np.array([
                         np.mean(ir_co2_meas_v1), np.mean(ir_co2_meas_v0),
                         ])

ir_calculated_co2 = np.array([
                         np.mean(ir_co2_calc_v1), np.mean(ir_co2_calc_v0),
                        ])


swe_observed_co2 = np.array([
                         np.mean(swe_co2_meas_v1), np.mean(swe_co2_meas_v0)])

swe_calculated_co2 = np.array([
                         np.mean(swe_co2_calc_v1), np.mean(swe_co2_calc_v0)])


fi_diff = np.mean(np.abs(fi_observed_co2-fi_calculated_co2))
print ('Fin co2', fi_diff)

es_diff = np.mean(np.abs(es_observed_co2-es_calculated_co2))
print ('Es co2', es_diff)

ir_diff = np.mean(np.abs(ir_observed_co2-ir_calculated_co2))
print ('Ir co2', ir_diff)

se_diff = np.mean(np.abs(swe_observed_co2-swe_calculated_co2))
print ('se co2', se_diff)



#----------------------------------
fi_observed_doc = np.array([np.mean(fin_doc_meas_v1), np.mean(fin_doc_meas_v0),
                         ])

fi_calculated_doc = np.array([np.mean(fin_doc_calc_v1), np.mean(fin_doc_calc_v0),
                         ])

es_observed_doc = np.array([
                         np.mean(es_doc_meas_v1), np.mean(es_doc_meas_v0),
                         ])

es_calculated_doc = np.array([
                         np.mean(es_doc_calc_v1), np.mean(es_doc_calc_v0),
                         ])

ir_observed_doc = np.array([
                         np.mean(ir_doc_meas_v1), np.mean(ir_doc_meas_v0),
                         ])

ir_calculated_doc = np.array([
                         np.mean(ir_doc_calc_v1), np.mean(ir_doc_calc_v0),
                        ])


swe_observed_doc = np.array([
                         np.mean(swe_doc_meas_v1), np.mean(swe_doc_meas_v0)])

swe_calculated_doc = np.array([
                         np.mean(swe_doc_calc_v1), np.mean(swe_doc_calc_v0)])


fi_diff = np.mean(np.abs(fi_observed_doc-fi_calculated_doc))
print ('Fin doc', fi_diff)

es_diff = np.mean(np.abs(es_observed_doc-es_calculated_doc))
print ('Es doc', es_diff)

ir_diff = np.mean(np.abs(ir_observed_doc-ir_calculated_doc))
print ('Ir doc', ir_diff)

se_diff = np.mean(np.abs(swe_observed_doc-swe_calculated_doc))
print ('se doc', se_diff)




