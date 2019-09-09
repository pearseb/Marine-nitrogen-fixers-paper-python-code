# -*- coding: utf-8 -*-
"""
Created on Wed May 16 10:17:32 2018

@author: pearseb
"""

#%% imporst
 
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmocean.cm as cmo
import seaborn as sb
sb.set(style='ticks')
import mpl_toolkits.basemap as bm
import pandas as pd 


# move to working directory
os.chdir("C://Users/pearseb/Dropbox/PhD/My articles/nitrogen-carbon cycles/data_for_publication")


#%% get data

data = nc.Dataset('BGC_Mk3Lpi_Fe-mod_default_N15active.nc', 'r')
mk3lpi_no3 = data.variables['no3'][...]
mk3lpi_n15 = data.variables['no3_15'][...]
mk3lpi_no3 = np.ma.masked_where(mk3lpi_no3<0.1, mk3lpi_no3)
mk3lpi_n15 = np.ma.masked_where(mk3lpi_no3<0.1, mk3lpi_n15)
mk3lpi_d15n = (mk3lpi_n15/(mk3lpi_no3-mk3lpi_n15)-1)*1000
mk3lpi_d15org = data.variables['sed_d15n'][...]

data = nc.Dataset('BGC_Mk3Lpi_Fe-2500per_default_N15active.nc', 'r')
mk3lpidust_no3 = data.variables['no3'][...]
mk3lpidust_n15 = data.variables['no3_15'][...]
mk3lpidust_no3 = np.ma.masked_where(mk3lpidust_no3<0.1, mk3lpidust_no3)
mk3lpidust_n15 = np.ma.masked_where(mk3lpidust_no3<0.1, mk3lpidust_n15)
mk3lpidust_d15n = (mk3lpidust_n15/(mk3lpidust_no3-mk3lpidust_n15)-1)*1000
mk3lpidust_d15org = data.variables['sed_d15n'][...]

data = nc.Dataset('BGC_Mk3Llgm_Fe-mod_default_N15active.nc', 'r')
mk3llgm_no3 = data.variables['no3'][...]
mk3llgm_n15 = data.variables['no3_15'][...]
mk3llgm_no3 = np.ma.masked_where(mk3llgm_no3<0.1, mk3llgm_no3)
mk3llgm_n15 = np.ma.masked_where(mk3llgm_no3<0.1, mk3llgm_n15)
mk3llgm_d15n = (mk3llgm_n15/(mk3llgm_no3-mk3llgm_n15)-1)*1000
mk3llgm_d15org = data.variables['sed_d15n'][...]

data = nc.Dataset('BGC_Mk3Llgm_Fe-2500per_default_N15active.nc', 'r')
mk3llgmdust_no3 = data.variables['no3'][...]
mk3llgmdust_n15 = data.variables['no3_15'][...]
mk3llgmdust_no3 = np.ma.masked_where(mk3llgmdust_no3<0.1, mk3llgmdust_no3)
mk3llgmdust_n15 = np.ma.masked_where(mk3llgmdust_no3<0.1, mk3llgmdust_n15)
mk3llgmdust_d15n = (mk3llgmdust_n15/(mk3llgmdust_no3-mk3llgmdust_n15)-1)*1000
mk3llgmdust_d15org = data.variables['sed_d15n'][...]

grid = nc.Dataset('grid_spec_mk3l_128_112_21_v2.nc', 'r')
dvts = grid.variables['dvts'][...]
dats = grid.variables['dats'][...]
lats = grid.variables['latts'][...]
lat_bnds = grid.variables['latts_bnds'][...]
lons = grid.variables['lonts'][...]
lon_bnds = grid.variables['lonts_bnds'][...]
lon_bnds[:,0] += 360.0
deps = grid.variables['zts'][...]
dep_bnds = grid.variables['zts_bnds'][...]
zts = dvts/dats



#%% apply depth correction to d15N of organic matter (see Robinson et al., 2012, Paleoceanography)

deps3d = np.cumsum(zts,axis=0)

# correction
mk3lpi_d15org_cor = mk3lpi_d15org + 0.9*(deps3d*1e-3)
mk3lpidust_d15org_cor = mk3lpidust_d15org + 0.9*(deps3d*1e-3)

mk3llgm_d15org_cor = mk3llgm_d15org + 0.9*(deps3d*1e-3)
mk3llgmdust_d15org_cor = mk3llgmdust_d15org + 0.9*(deps3d*1e-3)


# average over all depths
mk3lpi_d15org_corz = np.ma.average(mk3lpi_d15org_cor, axis=0, weights=zts)
mk3lpidust_d15org_corz = np.ma.average(mk3lpidust_d15org_cor, axis=0, weights=zts)

mk3llgm_d15org_corz = np.ma.average(mk3llgm_d15org_cor, axis=0, weights=zts)
mk3llgmdust_d15org_corz = np.ma.average(mk3llgmdust_d15org_cor, axis=0, weights=zts)



#%% collect prepared compilation of sedimentary d15N records

df = pd.read_csv('metadata_for_cores.csv')
print(df)

records = df[~np.isnan(df['d15n_LateH'])]
bulk_records = records[records['type']=='bulk']
bound_records = records[records['type']!='bulk']


#%%

lat_labs = ['80$^{\circ}$S', '60$^{\circ}$S', '40$^{\circ}$S', '20$^{\circ}$S', '0$^{\circ}$', \
        '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N', '80$^{\circ}$N']
lon_labs = ['0$^{\circ}$E', '50$^{\circ}$E', '100$^{\circ}$E', '150$^{\circ}$E', '200$^{\circ}$E', \
        '250$^{\circ}$E', '300$^{\circ}$E', '350$^{\circ}$E']

domain = [-45,0,45,355]                
domain_draw = [-40,0,40,360]
dlat=20
dlon=60

xx,yy = np.meshgrid(lons, lats)


#%%

levs = np.arange(-5,5.1,0.5)
conts = [-1,1]

fig = plt.figure(facecolor='w', figsize=(10,4))

plt.title('Glacial minus Late Holocene change in $\delta^{15}$N$_{org}$', family='sans-serif', fontsize=12)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='c')
lonproj, latproj = proj(xx, yy)

bulk_x, bulk_y = proj(np.array(bulk_records['lon']),np.array(bulk_records['lat']))
bound_x, bound_y = proj(np.array(bound_records['lon']),np.array(bound_records['lat']))

proj.drawcoastlines(linewidth=0.5, color='k')
proj.fillcontinents(color='grey')
p3 = plt.contourf(lonproj, latproj, mk3lpidust_d15org_corz-mk3lpi_d15org_corz, cmap=cmo.balance, corner_mask=False, \
             levels=levs, vmin=np.ma.min(levs), vmax=np.ma.max(levs), extend='both')
c3 = plt.contour(lonproj, latproj, mk3lpidust_d15org_corz-mk3lpi_d15org_corz, colors='black', levels=conts, alpha=0.8, linewidths=0.5, linestyle='-')
s31 = plt.scatter(bound_x, bound_y, s=150, c=bound_records['d15n_LGM']-bound_records['d15n_LateH'], \
                 marker='*', vmin=np.ma.min(levs), vmax=np.ma.max(levs), cmap=cmo.balance, \
                 alpha=0.75, edgecolor='k', linewidths=1.0, zorder=3)
s32 = plt.scatter(bulk_x, bulk_y, s=40, c=bulk_records['d15n_LGM']-bulk_records['d15n_LateH'], \
                 marker='o', vmin=np.ma.min(levs), vmax=np.ma.max(levs), cmap=cmo.balance, \
                 alpha=0.75, edgecolor='k', linewidths=1.0, zorder=2)

proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], color=(.3,.3,.3), linewidth=0, fontsize=12, family='sans-serif')
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[True,False,False,True], color=(.3,.3,.3), linewidth=0, fontsize=12)


from matplotlib.lines import Line2D
elements = [Line2D([0], [0], marker='o', markerfacecolor='w', markeredgecolor='k', color='w', linewidth=1.0, markersize=15, label='Bulk organic N'),\
            Line2D([0], [0], marker='*', markerfacecolor='w', markeredgecolor='k', color='w', linewidth=1.0, markersize=20, label='Bound organic N')]

plt.legend(handles=elements, loc='center', bbox_to_anchor=(0.5,-0.25), ncol=2, frameon=False, scatterpoints=1)

plt.subplots_adjust(bottom=0.1, top=0.95, left=0.075, right=0.85)
cbax = fig.add_axes([0.88, 0.25, 0.03, 0.55])
cbar = plt.colorbar(p3, cax=cbax, orientation='vertical')
cbar.ax.set_ylabel(u'$\delta^{15}$N ‰ vs air', fontsize=12, family='sans-serif')

plt.clabel(c3, manual=True, fmt='%i', fontsize=10, colors='k', inline=True)


fig.savefig('figures_for_publication/fig4.pdf',dpi=300,bbox_inches='tight')
fig.savefig('figures_for_publication/fig4.png',dpi=300,bbox_inches='tight')



