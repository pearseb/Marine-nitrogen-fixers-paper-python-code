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
import matplotlib as mpl
import seaborn as sb
sb.set(style='ticks')

os.chdir('C:\\Users\\pearseb\\Dropbox\\PhD\\My articles\\nitrogen-carbon cycles\\data_for_publication')


#%% get grid information

grid = nc.Dataset('grid_spec_mk3l_128_112_21_v2.nc','r')
dvts = grid.variables['dvts'][...]
dats = grid.variables['dats'][...]
lats = grid.variables['latts'][...]
lons = grid.variables['lonts'][...]


#%% collect final nfixation rates and pco2 concentrations

fixT = []
spo4 = []
sfe = []
dicT = []
adicT = []

d = {}
states = ['GFDLpi_', 'Mk3Lpi_', 'HadGEMpi_', 'Mk3Llgm_']
numbers = [ 'Fe-50per_default','Fe-mod_KfeHigh','Fe-80per_default','Fe-mod_default','Fe-mod_KfeLow','Fe-500per_default','Fe-2500per_default']

for aaa in states:
    for num in numbers:
        fn = 'BGC_'+aaa+num+'_openCO2.nc'
        print(fn)
        data = nc.Dataset(fn, 'r')
        d["%s%s_no3"%(aaa,num)] = data.variables['no3'][...]*1e3
        mask = np.ma.getmask(d['%s%s_no3'%(aaa,num)])
        d["%s%s_nfix"%(aaa,num)] = data.variables['n2fix'][...]
        d["%s%s_po4"%(aaa,num)] = data.variables['po4'][...]*1e3
        d["%s%s_fe"%(aaa,num)] = data.variables['dfe'][...]*1e6
        d["%s%s_dic"%(aaa,num)] = data.variables['dissic'][...]*1e3
        d["%s%s_adic"%(aaa,num)] = data.variables['dissic_abio'][...]*1e3
        
        fix = np.ma.masked_where(mask[0,:,:], d["%s%s_nfix"%(aaa,num)])
        fixT.append(np.ma.sum(fix*dats)*86400*365*14*1e-12)
        dicT.append(np.ma.sum(d["%s%s_dic"%(aaa,num)]*dvts)*12*1e-18)
        adicT.append(np.ma.sum(d["%s%s_adic"%(aaa,num)]*dvts)*12*1e-18)
        spo4.append(np.ma.average(d["%s%s_po4"%(aaa,num)][0,:,:], weights=dats))
        sfe.append(np.ma.average(d["%s%s_fe"%(aaa,num)][0,:,:], weights=dats))


pco2 = []

for aaa in states:
    for num in numbers:
        fn = 'transientCO2_'+aaa+num+'.txt'
        print(fn)
        data = np.genfromtxt(fn, usecols=(3,4))
        pco2.append(data[-1,1])
        
bioC = np.array(dicT)-np.array(adicT)


#%% fit lines

from scipy.optimize import curve_fit

def func(x,a,b):
    return a*x + b


#%% convert lists to arrays


fixT = np.array(fixT)
pco2 = np.array(pco2)
sfe = np.array(sfe)


#%% make figure

fig = plt.figure(facecolor='w', figsize=(12,8))

col2 = 'navajowhite'
alf1 = 0.8; alf2 = 0.5
size1 = 200; size2 = 300; size3 = 400; size4=1250
dot1 = 50; dot2 = 50; dot3 = 75

cmap = cmo.amp
# extract all colors from the map
cmaplist = [cmap(i) for i in range(cmap.N)]
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.array([0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.7, 1.3])
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

ax10 = plt.gca()
ax10.spines["top"].set_visible(False)
ax10.spines["bottom"].set_visible(True)
ax10.spines["right"].set_visible(False)
ax10.spines["left"].set_visible(True)
ax10.get_xaxis().tick_bottom()
ax10.get_yaxis().tick_left()
ax10.tick_params(axis='both', which='both', bottom='on', labelbottom='on', left='on', labelleft='on')

femin=0.2; femax=1.4

plt.scatter(fixT[0],pco2[0], c=cmap(norm(sfe[0])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='chocolate', linewidths=4.0)
plt.scatter(fixT[1],pco2[1], c=cmap(norm(sfe[1])), s=size3, marker='^', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='chocolate', linewidths=4.0)
plt.scatter(fixT[2],pco2[2], c=cmap(norm(sfe[2])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='chocolate', linewidths=4.0)
plt.scatter(fixT[3],pco2[3], c=cmap(norm(sfe[3])), s=size4, marker='*', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='chocolate', linewidths=4.0)
plt.scatter(fixT[4],pco2[4], c=cmap(norm(sfe[4])), s=size3, marker='v', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='chocolate', linewidths=4.0)
plt.scatter(fixT[5],pco2[5], c=cmap(norm(sfe[5])), s=size4, marker='P', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='chocolate', linewidths=4.0)
plt.scatter(fixT[6],pco2[6], c=cmap(norm(sfe[6])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='chocolate', linewidths=4.0)

plt.scatter(fixT[7], pco2[7], c=cmap(norm(sfe[7])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='k', linewidths=4.0)
plt.scatter(fixT[8], pco2[8], c=cmap(norm(sfe[8])), s=size3, marker='^', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='k', linewidths=4.0)
plt.scatter(fixT[9], pco2[9], c=cmap(norm(sfe[9])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='k', linewidths=4.0)
plt.scatter(fixT[10],pco2[10], c=cmap(norm(sfe[10])), s=size4, marker='*', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='k', linewidths=4.0)
plt.scatter(fixT[11],pco2[11], c=cmap(norm(sfe[11])), s=size3, marker='v', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='k', linewidths=4.0)
plt.scatter(fixT[12],pco2[12], c=cmap(norm(sfe[12])), s=size4, marker='P', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='k', linewidths=4.0)
plt.scatter(fixT[13],pco2[13], c=cmap(norm(sfe[13])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='k', linewidths=4.0)

plt.scatter(fixT[14],pco2[14], c=cmap(norm(sfe[14])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='dimgrey', linewidths=4.0)
plt.scatter(fixT[15],pco2[15], c=cmap(norm(sfe[15])), s=size3, marker='^', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='dimgrey', linewidths=4.0)
plt.scatter(fixT[16],pco2[16], c=cmap(norm(sfe[16])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='dimgrey', linewidths=4.0)
plt.scatter(fixT[17],pco2[17], c=cmap(norm(sfe[17])), s=size4, marker='*', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='dimgrey', linewidths=4.0)
plt.scatter(fixT[18],pco2[18], c=cmap(norm(sfe[18])), s=size3, marker='v', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='dimgrey', linewidths=4.0)
plt.scatter(fixT[19],pco2[19], c=cmap(norm(sfe[19])), s=size4, marker='P', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='dimgrey', linewidths=4.0)
plt.scatter(fixT[20],pco2[20], c=cmap(norm(sfe[20])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='dimgrey', linewidths=4.0)

plt.scatter(fixT[21],pco2[21], c=cmap(norm(sfe[21])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='royalblue', linewidths=4.0)
plt.scatter(fixT[22],pco2[22], c=cmap(norm(sfe[22])), s=size3, marker='^', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='royalblue', linewidths=4.0, label='Fe limitation stronger ($K_{Fe}^{D}$=0.5)')
plt.scatter(fixT[23],pco2[23], c=cmap(norm(sfe[22])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='royalblue', linewidths=4.0)
plt.scatter(fixT[24],pco2[24], c=cmap(norm(sfe[23])), s=size4, marker='*', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='royalblue', linewidths=4.0, label='Modern Fe flux')
plt.scatter(fixT[25],pco2[25], c=cmap(norm(sfe[24])), s=size3, marker='v', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='royalblue', linewidths=4.0, label='Fe limitation weaker ($K_{Fe}^{D}$=0.1)')
plt.scatter(fixT[26],pco2[26], c=cmap(norm(sfe[25])), s=size4, marker='P', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='royalblue', linewidths=4.0, label='Glacial Fe flux')
sc = plt.scatter(fixT[27],pco2[27], c=cmap(norm(sfe[26])), s=size3, marker='o', zorder=1, alpha=alf1, vmin=femin, vmax=femax, edgecolor='royalblue', linewidths=4.0)


# create lines
slopes=[];error=[]

popt, pcov = curve_fit(func, np.array(fixT[0:7]), np.array(pco2[0:7]))
print("a = ",popt[0], " b = ",popt[1])
print("one standard deviation = ", np.sqrt(np.diag(pcov)))
slopes.append(popt[0]); error.append(np.sqrt(np.diag(pcov))[0])
x = np.arange(40, 240.0, 0.5)
plt.plot(x, func(x,popt[0],popt[1]),color='chocolate',alpha=1.0,zorder=0, linewidth=2.0)

popt, pcov = curve_fit(func, np.array(fixT[7:14]), np.array(pco2[7:14]))
print( "a = ",popt[0], " b = ",popt[1])
print( "one standard deviation = ", np.sqrt(np.diag(pcov)))
slopes.append(popt[0]); error.append(np.sqrt(np.diag(pcov))[0])
x = np.arange(40, 240.0, 0.5)
plt.plot(x, func(x,popt[0],popt[1]),color='black',alpha=1.0,zorder=0, linewidth=2.0)

popt, pcov = curve_fit(func, np.array(fixT[14:21]), np.array(pco2[14:21]))
print("a = ",popt[0], " b = ",popt[1])
print("one standard deviation = ", np.sqrt(np.diag(pcov)))
slopes.append(popt[0]); error.append(np.sqrt(np.diag(pcov))[0])
x = np.arange(40, 240.0, 0.5)
plt.plot(x, func(x,popt[0],popt[1]),color='dimgrey',alpha=1.0,zorder=0, linewidth=2.0)

popt, pcov = curve_fit(func, np.array(fixT[21:28]), np.array(pco2[21:28]))
print("a = ",popt[0], " b = ",popt[1])
print("one standard deviation = ", np.sqrt(np.diag(pcov)))
slopes.append(popt[0]); error.append(np.sqrt(np.diag(pcov))[0])
x = np.arange(40, 240.0, 0.5)
plt.plot(x, func(x,popt[0],popt[1]),color='royalblue',alpha=1.0,zorder=0, linewidth=2.0)

plt.text(63,279,r'Glacial Fe-induced $\Delta CO_2$', fontweight='bold', family='sans-serif', fontsize=15, color='k', ha='left', va='center', rotation=0)
plt.plot((63,100),(277.75,277.75),'k')
plt.text(63,275,r'GFDL$^{warm}$ = %.1f ppm'%(pco2[5]-pco2[3]), fontweight='bold', family='sans-serif', fontsize=15, color='chocolate', ha='left', va='center', rotation=0)
plt.text(63,271,r'Mk3L$^{mild}$ = %.1f ppm'%(pco2[12]-pco2[10]), fontweight='bold', family='sans-serif', fontsize=15, color='k', ha='left', va='center', rotation=0)
plt.text(63,267,r'HadGEM$^{cool}$ = %.1f ppm'%(pco2[19]-pco2[17]), fontweight='bold', family='sans-serif', fontsize=15, color='dimgrey', ha='left', va='center', rotation=0)
plt.text(63,263,r'Mk3L$^{cold}$ = %.1f ppm'%(pco2[26]-pco2[24]), fontweight='bold', family='sans-serif', fontsize=15, color='royalblue', ha='left', va='center', rotation=0)


plt.xlabel('N$_2$fixation (Tg N yr$^{-1}$)', fontsize=18, family='sans-serif')
plt.xticks(np.arange(60,201,20),np.arange(60,201,20),family='sans-serif', fontsize=15)
plt.xlim(60,200)

plt.ylabel('Atmospheric pCO$_2$ (ppm)', fontsize=18, family='sans-serif')
plt.yticks(np.arange(260,321,10),np.arange(260,321,10),family='sans-serif', fontsize=15)
plt.ylim(260,320)


# show mean slope of N2fix-CO2 relationship
print(np.mean(slopes))
print(np.mean(error))


# make legend
from matplotlib.lines import Line2D
legend_elements = [Line2D([0],[0], marker='^', markersize=15, color='w', markerfacecolor='k', linewidth=4.0, label='Fe limitation stronger ($K_{Fe}^{D}$ = 0.5)'),
                   Line2D([0],[0], marker='v', markersize=15, color='w', markerfacecolor='k', linewidth=4.0, label='Fe limitation weaker ($K_{Fe}^{D}$ = 0.1)'),
                   Line2D([0],[0], marker='*', markersize=25, color='w', markerfacecolor='k', linewidth=4.0, label='Modern Fe supply'),
                   Line2D([0],[0], marker='P', markersize=20, color='w', markerfacecolor='k', linewidth=4.0, label='Glacial Fe supply (500%)'),
                   Line2D([0],[0], marker='o', markersize=15, color='w', markerfacecolor='k', linewidth=4.0, label='Other Fe supply (50%, 80% or 2500%)')]
plt.legend(handles=legend_elements, loc='upper right', frameon=False, fontsize=13, scatterpoints=1, handletextpad=0.2)


plt.subplots_adjust(bottom=0.25, top=0.92, right=0.9)

# colorbar
cbax = fig.add_axes([0.125, 0.1, 0.775, 0.03])
cbar = mpl.colorbar.ColorbarBase(cbax, cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds, spacing='uniform', orientation='horizontal')
#cbar = plt.colorbar(sc, cax=cbax, orientation='horizontal')
cbar.ax.set_xlabel('global mean surface Fe ($\mu$mol m$^{-3}$)', family='sans-serif', fontsize=18)
plt.xticks(bounds, bounds, family='sans-serif', fontsize=15)


#%%

fig.savefig('figures_for_publication/fig2.pdf', dpi=300, bbox_inches='tight')
fig.savefig('figures_for_publication/fig2.png', dpi=300, bbox_inches='tight')



