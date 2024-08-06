# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset  


def create_profile_line(ax, wt, wtmin, sd, cols, ylabel, label, fs, facecolor, colorin, letter, hidex=False, hidey=False, elevation=None):
    if elevation is not None: ax.plot(elevation, color='brown', label = 'soil surface')
    ax.plot(wt, color=colorin, label = label)
    ax.fill_between(range(cols), wt+sd*2, wt-sd, color=colorin, alpha=0.075)
    if elevation is not None: 
        drainnorm = elevation - 0.35
    else:
        drainnorm = -0.35
    ax.hlines(y= drainnorm, xmin=0, xmax = cols, color='red',linestyles='--')
    ax.get_xaxis().set_visible(False) 
    ax.tick_params(axis='y', labelsize=fs)
    if elevation is None: ax.set_ylim([wtmin,0])
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.legend()
    ax.grid(visible=False)
    ax.set_facecolor(facecolor)
    if hidex: 
        ax.get_xaxis().set_visible(False) 
    else:
        ax.tick_params(axis='x', labelsize=fs)
    if hidey: 
        ax.get_yaxis().set_visible(False)
    else:
        ax.get_yaxis().set_visible(True)
    ax.text(0.03, 1, letter, horizontalalignment='left', 
                verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)

    return ax

def create_profile_boxplot(ax, datain, cols, colorin, title, label, fs, facecolor, letter='d', zero=False, hidex=True, hidey=False):
    
    df = pd.DataFrame(data=datain, columns=list(range(cols)))
    df.boxplot(ax = ax,
                  color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
                  boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
                  showfliers=False, grid=False, rot=0)

    if zero: ax.hlines(y= 0, xmin=0, xmax = cols, color='red',linestyles='--')
    meanval = df.mean(axis=0)
    meanval = np.round(meanval.mean(),2)
    title = title + ': mean '+ str(meanval)
    ax.set_title(title, fontsize = fs)
    ax.set_ylabel(label, fontsize=fs)
    ax.set_facecolor(facecolor)
    #ax.get_xaxis().set_visible(False) 
    if hidex: 
        ax.get_xaxis().set_visible(False) 
    else:
        ax.tick_params(axis='x', labelsize=fs)
    if hidey: 
        ax.get_yaxis().set_visible(False)
    else:
        ax.get_yaxis().set_visible(True)
        ax.set_ylabel(label, fontsize=fs)
        ax.tick_params(axis='y', labelsize=fs)
    ax.text(0.03, 1, letter, horizontalalignment='left', 
                verticalalignment='top', fontsize=18,  fontweight='bold', transform=ax.transAxes)

    return ax

def hydrology(ff, scens):
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    
    #-------------- Water tables cross section-------------------------------------
    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num='hydro', figsize=(15,18))   #width, height
    gs = gridspec.GridSpec(ncols=12, nrows=10, figure=fig, wspace=1.2, hspace=0.6)
    
    scen = scens[0]
    wt = np.mean(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    cols = np.shape(wt)[0]
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scens[1],:, :], axis = 0)
    wtmin = min(wtls) -0.7
    
 
    #----------Water tables time series
    
    wt = np.mean(ncf['strip']['dwt'][scen,:, :], axis = 1)
    days = np.shape(wt)[0]
    sd = np.std(ncf['strip']['dwt'][scen,:, :], axis = 1)# sdls = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    
    
    axwtts = fig.add_subplot(gs[0:2, :])
    axwtts.plot(range(days), wt, color='green', label = 'WT')
    #axwtts.fill_between(range(days), wt+sd*2, wt-sd*2, color='green', alpha=0.25)
    axwtts.hlines(y= -0.35, xmin=0, xmax = days, color='red',linestyles='--')
    for c in range(1,cols-1):
        axwtts.plot(range(days), ncf['strip']['dwt'][scen,:, c], alpha=0.2)
        
    #axwtts.get_xaxis().set_visible(False) 
    axwtts.tick_params(axis='y', labelsize=fs)
    axwtts.tick_params(axis='x', labelsize=fs)
    
    axwtts.set_ylim([wtmin,0])
    axwtts.set_ylabel('WT m', fontsize=fs)
    axwtts.legend(loc= 'lower left')
    axwtts.grid(visible=False)
    axwtts.set_facecolor(facecolor)
    axwtts.text(0.01, 1, 'a', horizontalalignment='left', 
                verticalalignment='top', fontsize=18,  fontweight='bold',transform=axwtts.transAxes)
    
    #-------- WT scens comparison-----------------
    wt0 = np.mean(ncf['strip']['dwtyr_growingseason'][scens[0],:, :], axis = 0)
    sd0 = np.std(ncf['strip']['dwtyr_growingseason'][scens[0],:, :], axis = 0)
    wtmin = min(wt0) - 0.7
    cols = np.shape(wt0)[0]
    #------------------------------
    ax = fig.add_subplot(gs[2:4, :4])
    ax = create_profile_line(ax, wt0, wtmin, sd0, cols, 'WT m', 'growing season', fs, facecolor, 'blue', 'b')
    #ax.text(0.0, 0.1, 'b', horizontalalignment='left', 
    #            verticalalignment='top', fontsize=18,  fontweight='bold')
    
    wt1 = np.mean(ncf['strip']['dwtyr_growingseason'][scens[1],:, :], axis = 0)
    sd1 = np.std(ncf['strip']['dwtyr_growingseason'][scens[1],:, :], axis = 0)
 
    ax = fig.add_subplot(gs[2:4, 4:8])
    ax = create_profile_line(ax, wt1, wtmin, sd1, cols, '','growing season', fs, facecolor, 'orange', 'c')

    deltawt = ncf['strip']['dwtyr_growingseason'][scens[1],:, :] - ncf['strip']['dwtyr_growingseason'][scens[0],:, :]
    ax = fig.add_subplot(gs[2:4, 8:])
    ax = create_profile_boxplot(ax, deltawt,cols,'green', 'WT difference', '', fs, facecolor, zero=True, letter = 'd')
 
    #---------------
    #-----------------HMW to ditch------------------------
    hmwtoditch0 = ncf['export']['hmwtoditch'][scens[0],:, :]
    ax = fig.add_subplot(gs[4:6, :4])
    ax = create_profile_boxplot(ax, hmwtoditch0, cols,'blue', 'HMW to ditch', '$kg \ ha^{-1} \ yr^{-1}$', fs, facecolor, letter = 'e')

    hmwtoditch1 = ncf['export']['hmwtoditch'][scens[1],:, :]
    ax = fig.add_subplot(gs[4:6, 4:8])
    ax = create_profile_boxplot(ax, hmwtoditch1, cols,'orange', 'HMW to ditch', '', fs, facecolor, letter = 'f')

    deltahmw = hmwtoditch1 - hmwtoditch0
    ax = fig.add_subplot(gs[4:6, 8:])
    ax = create_profile_boxplot(ax, deltahmw,cols,'green', 'HMW difference', '', fs, facecolor, zero=True, letter = 'g')

    #----------------LMW to ditch-------------------------
    lmwtoditch0 = ncf['export']['lmwtoditch'][scens[0],:, :]
    ax = fig.add_subplot(gs[6:8, :4])
    ax = create_profile_boxplot(ax, lmwtoditch0, cols,'blue', 'LMW to ditch', '$kg \ ha^{-1} \ yr^{-1}$', fs, facecolor, letter = 'h')

    lmwtoditch1 = ncf['export']['lmwtoditch'][scens[1],:, :]
    ax = fig.add_subplot(gs[6:8, 4:8])
    ax = create_profile_boxplot(ax, lmwtoditch1, cols,'orange', 'LMW to ditch', '', fs, facecolor, letter = 'i')

    deltalmw = lmwtoditch1 - lmwtoditch0
    ax = fig.add_subplot(gs[6:8, 8:])
    ax = create_profile_boxplot(ax, deltalmw,cols,'green', 'LMW difference', '', fs, facecolor, zero=True, letter = 'j')

    #------------Stand C balance---------------------------------
    site0 = ncf['balance']['C']['stand_c_balance_c'][scens[0], :, :]    
    ax = fig.add_subplot(gs[8:10, :4])
    ax = create_profile_boxplot(ax, site0, cols,'blue', 'Site C balance', '$kg \ ha^{-1} \ yr^{-1}$', fs, facecolor, letter = 'k')
    
    site1 = ncf['balance']['C']['stand_c_balance_c'][scens[1], :, :]
    ax = fig.add_subplot(gs[8:10, 4:8])
    ax = create_profile_boxplot(ax, site1, cols,'orange', 'Site C balance', '', fs, facecolor, letter = 'l')

    deltasoil = site1 - site0
    ax = fig.add_subplot(gs[8:10, 8:])
    ax = create_profile_boxplot(ax, deltasoil,cols,'green', 'Site C difference', '', fs, facecolor, zero=True, letter = 'm')


    ncf.close()
    

#ff=r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/paroninkorpi_doc/susi_paroninkorpi_60m.nc'

#ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/susi_ullikka_100m.nc'

ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/krycklan/susi_krycklan_100m.nc'

scens = [1,0]
hydrology(ff, scens)







