# -*- coding: utf-8 -*-
"""
Created on Mon May 14 08:47:59 2018

@author: pearseb
"""

#%%

import os 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
import seaborn as sb
sb.set(style='ticks')

# move to working directory
os.chdir("C://Users/pearseb/Dropbox/PhD/My articles/nitrogen-carbon cycles/data_for_publication")


#%% get data files

# grid
grid = nc.Dataset('grid_spec_mk3l_128_112_21_v2.nc','r')
dvts = grid.variables['dvts'][...]
dats = grid.variables['dats'][...]
lats = grid.variables['latts'][...]
lons = grid.variables['lonts'][...]


#%% get abiotic carbon distribution
        
data = nc.Dataset('BGC_Mk3Lpi_Fe-mod_default_N15active.nc','r')
abio = data.variables['dissic_abio'][...]
abio = np.ma.sum(abio*dvts)*12*1e-15


#%% get static simulations

expT_stat = []; dicT_stat = []; sfe_stat = []; spo4_stat = [];

ds = {}
files = ['Fe-25per', 'Fe-50per', 'Fe-75per', 'Fe-mod', 'Fe-150per', 'Fe-300per', 'Fe-500per', 'Fe-2500per']

for fi in files:
    fn = 'BGC_Mk3Lpi_'+fi+'_static22.nc'
    print(fn)
    data = nc.Dataset(fn, 'r')
    ds["%s_no3"%(fi)] = data.variables['no3'][...]*1e3
    mask = np.ma.getmask(ds['%s_no3'%(fi)])
    ds["%s_exp"%(fi)] = (data.variables['epc100'][...] + data.variables['epcalc100'][...])*12*86400*365
    ds["%s_dic"%(fi)] = data.variables['dissic'][...]*1e3
    ds["%s_fe"%(fi)] = data.variables['dfe'][...]*1e6
    ds["%s_po4"%(fi)] = data.variables['po4'][...]*1e3
    
    expT_stat.append(np.ma.sum(ds["%s_exp"%(fi)]*dats)*1e-15)
    dicT_stat.append(np.ma.sum(ds["%s_dic"%(fi)]*dvts)*12*1e-18)
    sfe_stat.append(np.ma.average(ds["%s_fe"%(fi)][0,:,:], weights=dats))
    spo4_stat.append(np.ma.average(ds["%s_po4"%(fi)][0,:,:], weights=dats))

bioT_stat = np.array(dicT_stat)-abio


#%% get phosphate limitation simulations

expT_plim = []; dicT_plim = []; sfe_plim = []; spo4_plim = []; fixT_plim = []

ds = {}
files = ['Fe-50per', 'Fe-75per', 'Fe-mod', 'Fe-150per', 'Fe-500per', 'Fe-2500per']

for fi in files:
    fn = 'BGC_Mk3Lpi_'+fi+'_Kpo40.1forN2fixers.nc'
    print(fn)
    data = nc.Dataset(fn, 'r')
    ds["%s_no3"%(fi)] = data.variables['no3'][...]*1e3
    mask = np.ma.getmask(ds['%s_no3'%(fi)])
    ds["%s_exp"%(fi)] = (data.variables['epc100'][...] + data.variables['epcalc100'][...])*12*86400*365
    ds["%s_dic"%(fi)] = data.variables['dissic'][...]*1e3
    ds["%s_fe"%(fi)] = data.variables['dfe'][...]*1e6
    ds["%s_po4"%(fi)] = data.variables['po4'][...]*1e3
    ds["%s_nfix"%(fi)] = np.ma.masked_where(mask[0,:,:],data.variables['n2fix'][...])

    expT_plim.append(np.ma.sum(ds["%s_exp"%(fi)]*dats)*1e-15)
    dicT_plim.append(np.ma.sum(ds["%s_dic"%(fi)]*dvts)*12*1e-18)
    sfe_plim.append(np.ma.average(ds["%s_fe"%(fi)][0,:,:], weights=dats))
    spo4_plim.append(np.ma.average(ds["%s_po4"%(fi)][0,:,:], weights=dats))
    fixT_plim.append(np.ma.sum(ds["%s_nfix"%(fi)]*dats)*86400*365*14*1e-12)


bioT_plim = np.array(dicT_plim)-abio


#%% get C:P altered simulations

expT_cp165 = []; dicT_cp165 = []; sfe_cp165 = []; spo4_cp165 = []; fixT_cp165 = []

ds = {}
files = ['Fe-50per', 'Fe-75per', 'Fe-mod', 'Fe-150per', 'Fe-500per', 'Fe-2500per']

for fi in files:
    fn = 'BGC_Mk3Lpi_'+fi+'_C2P165forN2fixers.nc'
    print(fn)
    data = nc.Dataset(fn, 'r')
    ds["%s_no3"%(fi)] = data.variables['no3'][...]*1e3
    mask = np.ma.getmask(ds['%s_no3'%(fi)])
    ds["%s_exp"%(fi)] = (data.variables['epc100'][...] + data.variables['epcalc100'][...])*12*86400*365
    ds["%s_dic"%(fi)] = data.variables['dissic'][...]*1e3
    ds["%s_fe"%(fi)] = data.variables['dfe'][...]*1e6
    ds["%s_po4"%(fi)] = data.variables['po4'][...]*1e3
    ds["%s_nfix"%(fi)] = np.ma.masked_where(mask[0,:,:],data.variables['n2fix'][...])
    
    expT_cp165.append(np.ma.sum(ds["%s_exp"%(fi)]*dats)*1e-15)
    dicT_cp165.append(np.ma.sum(ds["%s_dic"%(fi)]*dvts)*12*1e-18)
    sfe_cp165.append(np.ma.average(ds["%s_fe"%(fi)][0,:,:], weights=dats))
    spo4_cp165.append(np.ma.average(ds["%s_po4"%(fi)][0,:,:], weights=dats))
    fixT_cp165.append(np.ma.sum(ds["%s_nfix"%(fi)]*dats)*86400*365*14*1e-12)

bioT_cp165 = np.array(dicT_cp165)-abio


#%% get C:P altered simulations

expT_cp0 = []; dicT_cp0 = []; sfe_cp0 = []; spo4_cp0 = []; fixT_cp0 = []

ds = {}
files = ['Fe-25per', 'Fe-50per', 'Fe-75per', 'Fe-mod', 'Fe-150per', 'Fe-300per', 'Fe-500per', 'Fe-2500per']

for fi in files:
    fn = 'BGC_Mk3Lpi_'+fi+'_noC2PforN2fixers.nc'
    print(fn)
    data = nc.Dataset(fn, 'r')
    ds["%s_no3"%(fi)] = data.variables['no3'][...]*1e3
    mask = np.ma.getmask(ds['%s_no3'%(fi)])
    ds["%s_exp"%(fi)] = (data.variables['epc100'][...] + data.variables['epcalc100'][...])*12*86400*365
    ds["%s_dic"%(fi)] = data.variables['dissic'][...]*1e3
    ds["%s_fe"%(fi)] = data.variables['dfe'][...]*1e6
    ds["%s_po4"%(fi)] = data.variables['po4'][...]*1e3
    ds["%s_nfix"%(fi)] = np.ma.masked_where(mask[0,:,:],data.variables['n2fix'][...])

    expT_cp0.append(np.ma.sum(ds["%s_exp"%(fi)]*dats)*1e-15)
    dicT_cp0.append(np.ma.sum(ds["%s_dic"%(fi)]*dvts)*12*1e-18)
    sfe_cp0.append(np.ma.average(ds["%s_fe"%(fi)][0,:,:], weights=dats))
    spo4_cp0.append(np.ma.average(ds["%s_po4"%(fi)][0,:,:], weights=dats))
    fixT_cp0.append(np.ma.sum(ds["%s_nfix"%(fi)]*dats)*86400*365*14*1e-12)

bioT_cp0 = np.array(dicT_cp0)-abio


#%% get default simulations

expT_def = []; dicT_def = []; sfe_def = []; spo4_def = []; fixT_def = []

ds = {}
files = ['Fe-25per', 'Fe-50per', 'Fe-75per', 'Fe-mod', 'Fe-150per', 'Fe-300per', 'Fe-500per', 'Fe-2500per']

for fi in files:
    if fi == 'Fe-2500per':
        fn = 'BGC_Mk3Lpi_'+fi+'_default_N15active.nc'
    elif fi == 'Fe-mod':
        fn = 'BGC_Mk3Lpi_'+fi+'_default_N15active.nc'
    else:
        fn = 'BGC_Mk3Lpi_'+fi+'_default.nc'
    print(fn)
    data = nc.Dataset(fn, 'r')
    ds["%s_no3"%(fi)] = data.variables['no3'][...]*1e3
    mask = np.ma.getmask(ds['%s_no3'%(fi)])
    ds["%s_exp"%(fi)] = (data.variables['epc100'][...] + data.variables['epcalc100'][...])*12*86400*365
    ds["%s_dic"%(fi)] = data.variables['dissic'][...]*1e3
    ds["%s_fe"%(fi)] = data.variables['dfe'][...]*1e6
    ds["%s_po4"%(fi)] = data.variables['po4'][...]*1e3
    ds["%s_nfix"%(fi)] = np.ma.masked_where(mask[0,:,:],data.variables['n2fix'][...])
    
    expT_def.append(np.ma.sum(ds["%s_exp"%(fi)]*dats)*1e-15)
    dicT_def.append(np.ma.sum(ds["%s_dic"%(fi)]*dvts)*12*1e-18)
    sfe_def.append(np.ma.average(ds["%s_fe"%(fi)][0,:,:], weights=dats))
    spo4_def.append(np.ma.average(ds["%s_po4"%(fi)][0,:,:], weights=dats))
    fixT_def.append(np.ma.sum(ds["%s_nfix"%(fi)]*dats)*86400*365*14*1e-12)

bioT_def = np.array(dicT_def)-abio



#%% source data (fig 3a)

print(sfe_stat,bioT_stat,spo4_stat)
print(sfe_cp0,bioT_cp0, spo4_cp0)
print(sfe_cp165,bioT_cp165, spo4_cp165)
print(sfe_plim,bioT_plim, spo4_plim)
print(sfe_def,bioT_def, spo4_def)

print(fixT_cp0)
print(fixT_cp165)
print(fixT_plim)
print(fixT_def)


#%% figure for paper

fig = plt.figure(facecolor='w', figsize=(14,6))
gs = GridSpec(10,12)

col2 = 'navajowhite'
alf1 = 0.8; alf2 = 0.5
size1 = 200; size2 = 200; size3 = 400
dot1 = 50; dot2 = 50; dot3 = 75


ax00 = plt.subplot(gs[0:8,0:7])
ax00.spines["top"].set_visible(False)
ax00.spines["bottom"].set_visible(True)
ax00.spines["right"].set_visible(False)
ax00.spines["left"].set_visible(True)
ax00.get_xaxis().tick_bottom()
ax00.get_yaxis().tick_left()
ax00.tick_params(axis='both', which='both', bottom='on', labelbottom='on', left='on', labelleft='on')

plt.plot(sfe_stat,bioT_stat,color='k', linewidth=1, zorder=0)
plt.plot(sfe_cp0,bioT_cp0,color='k', linewidth=1, zorder=0)
plt.plot(sfe_cp165,bioT_cp165,color='k', linewidth=1, zorder=0)
plt.plot(sfe_plim,bioT_plim,color='k', linewidth=1, zorder=0)
plt.plot(sfe_def,bioT_def,color='k', linewidth=1, zorder=0)

plt.scatter(sfe_stat,bioT_stat, cmap=cmo.amp, c=spo4_stat, s=size2, marker='o', zorder=1, alpha=alf1, vmin=0.25, vmax=0.75)
plt.scatter(sfe_stat,bioT_stat, c='k', s=dot1, marker='$1$', zorder=2, alpha=alf1)
plt.scatter(sfe_cp0,bioT_cp0, cmap=cmo.amp, c=spo4_cp0, s=size2, marker='o', zorder=1, alpha=alf1, vmin=0.25, vmax=0.75)
plt.scatter(sfe_cp0,bioT_cp0, c='k', s=dot1, marker='$2$', zorder=2, alpha=alf1)
plt.scatter(sfe_cp165,bioT_cp165, cmap=cmo.amp, c=spo4_cp165, s=size2, marker='o', zorder=1, alpha=alf1, vmin=0.25, vmax=0.75)
plt.scatter(sfe_cp165,bioT_cp165, c='k', s=dot1, marker='$3$', zorder=2, alpha=alf1)
plt.scatter(sfe_plim,bioT_plim, cmap=cmo.amp, c=spo4_plim, s=size2, marker='o', zorder=1, alpha=alf1, vmin=0.25, vmax=0.75)
plt.scatter(sfe_plim,bioT_plim, c='k', s=dot1, marker='$4$', zorder=2, alpha=alf1)
sc = plt.scatter(sfe_def,bioT_def, cmap=cmo.amp, c=spo4_def, s=size2, marker='o', zorder=1, alpha=alf1, vmin=0.25, vmax=0.75)
plt.scatter(sfe_def,bioT_def, c='k', s=dot1, marker='$5$', zorder=2, alpha=alf1)

plt.xlabel('global mean surface Fe ($\mu$mol m$^{-3}$)', fontsize=13, family='sans-serif')
plt.xticks(np.arange(0,1.41,0.2),np.array([0.0,0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4]),family='sans-serif', fontsize=13)
plt.xlim(0,1.4)

plt.ylabel('Respired C inventory (Pg C)', fontsize=13, family='sans-serif')
plt.yticks(np.arange(2000,3600,200),np.arange(2000,3600,200),family='sans-serif', fontsize=13)
plt.ylim(2000,3400)


plt.plot((sfe_def[3],sfe_def[3]), (2000,3400), linestyle='--', alpha=0.5, color='k',zorder=0)
plt.plot((sfe_def[-2],sfe_def[-2]), (2000,3400), linestyle='--', alpha=0.5, color='k',zorder=0)
plt.text(sfe_def[3],3450, 'Modern Fe flux', ha='center',va='center', family='sans-serif', fontsize=11, zorder=0) 
plt.text(sfe_def[-2],3450, 'Glacial Fe flux', ha='center',va='center', family='sans-serif', fontsize=11, zorder=0) 


plt.text(sfe_stat[-1]+0.05,bioT_stat[-1], 'No N$_2$ fixers', family='sans-serif', fontsize=11, ha='left', va='center')
plt.text(sfe_cp165[-1]+0.05,bioT_cp165[-1], 'N$_2$ fixers (C:P=165,  $K^D_{PO_4}$=10$^{-10}$)', family='sans-serif', fontsize=11, ha='left', va='center')
plt.text(sfe_plim[-1]+0.05,bioT_plim[-1], 'N$_2$ fixers (C:P=331,  $K^D_{PO_4}$=10$^{-1}$)', family='sans-serif', fontsize=11, ha='left', va='center')
plt.text(sfe_cp0[-1]+0.05,bioT_cp0[-1], 'N$_2$ fixers (C:P=0,  $K^D_{PO_4}$=10$^{-10}$)', family='sans-serif', fontsize=11, ha='left', va='center')
plt.text(sfe_def[-1]+0.05,bioT_def[-1], 'N$_2$ fixers (C:P=331,  $K^D_{PO_4}$=10$^{-10}$)', family='sans-serif', fontsize=11, ha='left', va='center')


#%%


ax10 = plt.subplot(gs[0:8,9::])
ax10.spines["top"].set_visible(False)
ax10.spines["bottom"].set_visible(True)
ax10.spines["right"].set_visible(True)
ax10.spines["left"].set_visible(False)
ax10.tick_params(axis='both', which='both', bottom=True, labelbottom=True, right=True, labelright=True, left=False, labelleft=False)


plt.scatter(fixT_cp0,bioT_cp0, cmap=cmo.amp, c=spo4_cp0, s=size2, marker='o', zorder=1, alpha=alf1, vmin=0.25, vmax=0.75)
plt.scatter(fixT_cp0,bioT_cp0, c='k', s=dot1, marker='$2$', zorder=2, alpha=alf1)
plt.scatter(fixT_cp165,bioT_cp165, cmap=cmo.amp, c=spo4_cp165, s=size2, marker='o', zorder=1, alpha=alf1, vmin=0.25, vmax=0.75)
plt.scatter(fixT_cp165,bioT_cp165, c='k', s=dot1, marker='$3$', zorder=2, alpha=alf1)
plt.scatter(fixT_plim,bioT_plim, cmap=cmo.amp, c=spo4_plim, s=size2, marker='o', zorder=1, alpha=alf1, vmin=0.25, vmax=0.75)
plt.scatter(fixT_plim,bioT_plim, c='k', s=dot1, marker='$4$', zorder=2, alpha=alf1)
colour = plt.scatter(fixT_def,bioT_def, cmap=cmo.amp, c=spo4_def, s=size2, marker='o', zorder=1, alpha=alf1, vmin=0.25, vmax=0.75)
plt.scatter(fixT_def,bioT_def, c='k', s=dot1, marker='$5$', zorder=2, alpha=alf1)

plt.plot(fixT_cp0,bioT_cp0, c='k', zorder=0, alpha=alf1)
plt.plot(fixT_cp165,bioT_cp165, c='k', zorder=0, alpha=alf1)
plt.plot(fixT_plim,bioT_plim, c='k', zorder=0, alpha=alf1)
plt.plot(fixT_def,bioT_def, c='k', zorder=0, alpha=alf1)


plt.plot((fixT_def[3],fixT_def[3]), (1800.,3400), linestyle='--', alpha=alf2, color='black' ,zorder=0)
plt.text(fixT_def[3]+5, 2050, 'PO$_4$ limited', fontsize=11, family='sans-serif', rotation=90, ha='left', va='bottom')
plt.text(fixT_def[3]-2, 2050, 'NO$_3$ limited', fontsize=11, family='sans-serif', rotation=90, ha='right', va='bottom')

xx, yy = np.meshgrid(np.arange(50,201,1), np.arange(2000,3401,10))
zz = xx-fixT_def[3]
plt.contourf(xx,yy,zz, cmap=cmo.delta_r, alpha=0.3, zorder=0, levels=np.arange(-150,151,1))

plt.ylim(2000,3400)
plt.xlim(50,200)

plt.ylabel('Respired C inventory (Pg C)', fontsize=13, family='sans-serif')
ax10.yaxis.set_label_position('right')
plt.yticks(np.arange(2200,3401,200),np.arange(2200,3401,200),family='sans-serif', fontsize=13)
plt.xlabel('N$_2$ fixation (Tg N yr $^{-1}$)', fontsize=13, family='sans-serif')
plt.xticks(np.arange(60,200,40),np.arange(60,201,40),family='sans-serif', fontsize=13)


plt.text(0.05,0.925, 'a)', transform=ax00.transAxes, fontsize=15, family='sans-serif')
plt.text(0.05,0.925, 'b)', transform=ax10.transAxes, fontsize=15, family='sans-serif')


plt.subplots_adjust(top=0.90, left=0.08, right=0.92, bottom=0.2)


cbax = fig.add_axes([0.1, 0.12, 0.8, 0.04])
cbar = plt.colorbar(sc, cax=cbax, orientation='horizontal', ticks=np.arange(0.3,0.8,0.1))
cbar.ax.set_xlabel('global mean surface PO$_4$ (mmol m$^{-3}$)', family='sans-serif', fontsize=13)
plt.xticks(np.arange(0.3,0.8,0.1), np.array([0.3,0.4,0.5,0.6,0.7]), family='sans-serif', fontsize=12)


#%%

fig.savefig('figures_for_publication/fig3.pdf', dpi=300, bbox_inches='tight')
fig.savefig('figures_for_publication/fig3.png', dpi=300, bbox_inches='tight')

