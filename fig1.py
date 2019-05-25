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
import matplotlib
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
import mpl_toolkits.basemap as bm
import seaborn as sb
sb.set(style='ticks')

# move to working directory
os.chdir("C://Users/pearseb/Dropbox/PhD/My articles/nitrogen-carbon cycles/data_for_publication")


#%% get data files

# get grid
grid = nc.Dataset('grid_spec_mk3l_128_112_21_v2.nc','r')
dvts = grid.variables['dvts'][...]
dats = grid.variables['dats'][...]
lats = grid.variables['latts'][...]
lons = grid.variables['lonts'][...]


#%% get the non-static N cycle datafiles

data = nc.Dataset('BGC_Mk3Lpi_Fe-mod_default_N15active.nc')

mod_no3 = data.variables['no3'][...]*1e3
mask = np.ma.getmask(mod_no3)

# collect variables
mod_nfix = data.variables['n2fix'][...]*1e3*86400*365
mod_wcden = data.variables['wcden'][...]*1e3*86400*365
mod_exp = (data.variables['epc100'][...] + data.variables['epcalc100'][...])*12*86400*365
mod_epp = data.variables['epp100'][...]*1e3*86400*365
mod_dic = data.variables['dissic'][...]*1e3
mod_po4 = data.variables['po4'][...]*1e3
mod_oxy = data.variables['o2'][...]*1e3
mod_cp = data.variables['CtoP'][...]

# calculate suboxic zone depth
mod_subox = np.ma.masked_where(mod_oxy>=10.0, mod_oxy)
mod_subox = mod_subox*0+1
mod_subox = np.ma.sum(mod_subox*(dvts/dats),axis=0)

# mask variables
mod_nfix = np.ma.masked_where(mask[0,:,:], mod_nfix)
mod_wcden = np.ma.masked_where(mask[0,:,:], mod_wcden)
mod_cp = np.ma.masked_where(mask[0,:,:], mod_cp)

# calculate totals and averages
mod_nfixT = np.ma.sum(mod_nfix*dats)*14*1e-3*1e-12              # Tg N yr-1
mod_nfixC = np.ma.sum(mod_nfix*dats)*(331/50.0)*12*1e-3*1e-15   # Pg C yr-1
mod_wcdenT = np.ma.sum(mod_wcden*dats)*14*1e-3*1e-12            # Tg N yr-1
mod_expT = np.ma.sum(mod_exp*dats)*1e-15 + mod_nfixC            # Pg C yr-1
mod_dicZ = np.ma.sum(mod_dic*(dvts/dats), axis=0)*12*1e-6       # tonnes C m-2
mod_dicT = np.ma.sum(mod_dic*dvts)*12*1e-18                     # Pg C
mod_spo4 = np.ma.average(mod_po4[0,:,:], weights=dats)          # mmol PO4 / m3
mod_cpT = np.ma.average(mod_cp, weights=dats)



data = nc.Dataset('BGC_Mk3Lpi_Fe-500per_default.nc')

fe500_no3 = data.variables['no3'][...]*1e3
mask = np.ma.getmask(fe500_no3)

# collect variables
fe500_nfix = data.variables['n2fix'][...]*1e3*86400*365
fe500_wcden = data.variables['wcden'][...]*1e3*86400*365
fe500_exp = (data.variables['epc100'][...] + data.variables['epcalc100'][...])*12*86400*365
fe500_epp = data.variables['epp100'][...]*1e3*86400*365
fe500_dic = data.variables['dissic'][...]*1e3
fe500_po4 = data.variables['po4'][...]*1e3
fe500_oxy = data.variables['o2'][...]*1e3
fe500_cp = data.variables['CtoP'][...]

# calculate suboxic zone depth
fe500_subox = np.ma.masked_where(fe500_oxy>=10.0, fe500_oxy)
fe500_subox = fe500_subox*0+1
fe500_subox = np.ma.sum(fe500_subox*(dvts/dats),axis=0)

# mask variables
fe500_nfix = np.ma.masked_where(mask[0,:,:], fe500_nfix)
fe500_wcden = np.ma.masked_where(mask[0,:,:], fe500_wcden)
fe500_cp = np.ma.masked_where(mask[0,:,:], fe500_cp)

# calculate totals and averages
fe500_nfixT = np.ma.sum(fe500_nfix*dats)*14*1e-3*1e-12              # Tg N yr-1
fe500_nfixC = np.ma.sum(fe500_nfix*dats)*(331/50.0)*12*1e-3*1e-15   # Pg C yr-1
fe500_wcdenT = np.ma.sum(fe500_wcden*dats)*14*1e-3*1e-12            # Tg N yr-1
fe500_expT = np.ma.sum(fe500_exp*dats)*1e-15 + fe500_nfixC            # Pg C yr-1
fe500_dicZ = np.ma.sum(fe500_dic*(dvts/dats), axis=0)*12*1e-6       # tonnes C m-2
fe500_dicT = np.ma.sum(fe500_dic*dvts)*12*1e-18                     # Pg C
fe500_spo4 = np.ma.average(fe500_po4[0,:,:], weights=dats)          # mmol PO4 / m3
fe500_cpT = np.ma.average(fe500_cp, weights=dats)

        
data = nc.Dataset('AGE_Mk3Lpi.nc','r')
age = data.variables['age'][...]


#%% create non-linear colormap for nfixation

from matplotlib.colors import LinearSegmentedColormap

class nlcmap(LinearSegmentedColormap):
    """A nonlinear colormap"""

    name = 'nlcmap'

    def __init__(self, cmap, levels):
        self.cmap = cmap
        self.monochrome = self.cmap.monochrome
        self.levels = np.asarray(levels, dtype='float64')
        self._x = self.levels/ self.levels.max()
        self.levmax = self.levels.max()
        self.levmin = self.levels.min()
        self._y = np.linspace(self.levmin, self.levmax, len(self.levels))

    def __call__(self, xi, alpha=1.0, **kw):
        yi = np.interp(xi, self._x, self._y)
        return self.cmap(yi/self.levmax, alpha)
        
        
levs = [0, 0.5, 1, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 350, 500]

cmap_lin = cm.viridis
cmap_nonlin = nlcmap(cmap_lin, levs)


#%% set the bounds of the figures

lat_labs = ['80$^{\circ}$S', '60$^{\circ}$S', '40$^{\circ}$S', '20$^{\circ}$S', '0$^{\circ}$', \
        '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N', '80$^{\circ}$N']
lon_labs = ['0$^{\circ}$E', '50$^{\circ}$E', '100$^{\circ}$E', '150$^{\circ}$E', '200$^{\circ}$E', \
        '250$^{\circ}$E', '300$^{\circ}$E', '350$^{\circ}$E']

domain = [-65,0,65,355]                
domain_draw = [-50,0,50,360]
dlat=20
dlon=60

xx,yy = np.meshgrid(lons, lats)


#%% make figure

fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, facecolor='w', figsize=(16,8))
plt.subplots_adjust(wspace=0.1)

levs0 = np.arange(-300,301,10)

ax00 = plt.axes(ax[0,0])
ax00.tick_params(top='off', right='off')
plt.title('$\Delta_{total}$ N$_2$ fixation = %i Tg N yr$^{-1}$'%(fe500_nfixT-mod_nfixT), family='sans-serif', fontsize=18)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='c')
lonproj, latproj = proj(xx, yy)

proj.drawcoastlines(linewidth=0.5, color='k')
proj.fillcontinents(color='grey')
plot0 = plt.contourf(lonproj, latproj, fe500_nfix-mod_nfix, cmap=cmo.balance, levels=levs0, vmin=np.min(levs0), vmax=np.max(levs0), extend='both')
stip = plt.contourf(lonproj, latproj, np.log10(fe500_wcden), colors='none', hatches=['....'])
lonproj, latproj = proj(xx[24::,:], yy[24::,:])
conts = plt.contour(lonproj, latproj, age[2,24::,:], levels=np.array([20]), colors='k', linewidths=1.0)

levs1 = np.arange(-0.3,0.31,0.02)

ax01 = plt.axes(ax[0,1])
ax01.tick_params(top='off', right='off')
plt.title('$\Delta_{ave}$ surface PO$_4$ = %.2f mmol m$^{-3}$'%(fe500_spo4-mod_spo4), family='sans-serif', fontsize=18)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='c')
lonproj, latproj = proj(xx, yy)

proj.drawcoastlines(linewidth=0.5, color='k')
proj.fillcontinents(color='grey')
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], color=(.3,.3,.3), linewidth=0, fontsize=15, family='sans-serif')
plot1 = plt.contourf(lonproj, latproj, fe500_po4[0,:,:]-mod_po4[0,:,:], cmap=cmo.balance, levels=levs1, vmin=np.min(levs1), vmax=np.max(levs1), extend='both')
cont1 = plt.contour(lonproj, latproj, fe500_cp-mod_cp, linewidths=1.0, levels=np.arange(-50,51,10), colors='k')
plt.clabel(cont1, inline=True, fmt='%i', manual=True, fontsize=14)


levs3 = np.arange(-50,50,2)

ax10 = plt.axes(ax[1,0])
ax10.tick_params(top='off', right='off')
plt.title('$\Delta_{total}$ C export = %.1f Pg C yr$^{-1}$'%(fe500_expT-mod_expT), family='sans-serif', fontsize=18)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='c')
lonproj, latproj = proj(xx, yy)
proj.drawcoastlines(linewidth=0.5, color='k')
proj.fillcontinents(color='grey')
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[True,False,False,True], color=(.3,.3,.3), linewidth=0, fontsize=15)
plot3 = plt.contourf(lonproj, latproj, fe500_exp-mod_exp, cmap=cmo.balance, levels=levs3, vmin=np.min(levs3), vmax=np.max(levs3), extend='both')
#cont3a = plt.contour(lonproj, latproj, fe500_epp-mod_epp, linewidths=1.0, levels=[-1], linestyles='--', colors='k')
#cont3b = plt.contour(lonproj, latproj, fe500_epp-mod_epp, linewidths=1.0, levels=[1], linestyles='-', colors='k')
stip = plt.contourf(lonproj, latproj, fe500_subox-mod_subox, levels=np.array([-10000,-500,500,10000]), colors='none', hatches=['/////', ' ', '....'])


levs5 = np.arange(-10,10,0.2)

ax11 = plt.axes(ax[1,1])
ax11.tick_params(top='off', right='off')
plt.title('$\Delta_{total}$ respired C = %i Pg C'%(fe500_dicT-mod_dicT), family='sans-serif', fontsize=18)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='c')
lonproj, latproj = proj(xx, yy)
proj.drawcoastlines(linewidth=0.5, color='k')
proj.fillcontinents(color='grey')
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], color=(.3,.3,.3), linewidth=0, fontsize=15, family='sans-serif')
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[True,False,False,True], color=(.3,.3,.3), linewidth=0, fontsize=15)
plot5 = plt.contourf(lonproj, latproj, fe500_dicZ-mod_dicZ, cmap=cmo.balance, levels=levs5, vmin=np.min(levs5), vmax=np.max(levs5), extend='both')

plt.subplots_adjust(wspace=0.05)


plt.text(0.025,0.035, 'coupling between N$_2$ fixation, upwelling (contours) and denitrification (dots)', \
         fontsize=12, style='italic', family='sans-serif', transform=ax00.transAxes)
plt.text(0.025,0.035, 'global phosphorus drawdown and increase in C:P ratios (contours)', \
         fontsize=12, style='italic', family='sans-serif', transform=ax01.transAxes)
plt.text(0.025,0.035, 'C-enriched organic matter exported through deeper suboxic zones (dots)', \
         fontsize=12, style='italic', family='sans-serif', transform=ax10.transAxes)
plt.text(0.025,0.035, 'increased storage of respired C in Pacific basin', \
         fontsize=12, style='italic', family='sans-serif', transform=ax11.transAxes)

plt.text(0.025,1.035, 'a)', \
         fontsize=18, family='sans-serif', transform=ax00.transAxes)
plt.text(0.025,1.035, 'b)', \
         fontsize=18, family='sans-serif', transform=ax01.transAxes)
plt.text(0.025,1.035, 'c)', \
         fontsize=18, family='sans-serif', transform=ax10.transAxes)
plt.text(0.025,1.035, 'd)', \
         fontsize=18, family='sans-serif', transform=ax11.transAxes)


#%%

cbax0 = fig.add_axes([0.085, 0.535, 0.02, 0.365])
cbar0 = plt.colorbar(plot0, cax=cbax0, orientation='vertical')
cbar0.ax.set_ylabel('N$_2$ fixation (mmol N m$^{-2}$ yr$^{-1}$)', family='sans-serif', fontsize=18)
cbar0.ax.yaxis.set_ticks_position('left')
cbar0.ax.yaxis.set_label_position('left')
plt.yticks(np.arange(-280,281,70),np.arange(-280,281,70),family='sans-serif', fontsize=15)

cbax1 = fig.add_axes([0.92, 0.535, 0.02, 0.365])
cbar1 = plt.colorbar(plot1, cax=cbax1, orientation='vertical')
cbar1.ax.set_ylabel('surface PO$_4$ (mmol m$^{-3}$)', family='sans-serif', fontsize=18)
lev = np.arange(-0.24,0.241,0.08); lev[3] = 0.0
plt.yticks(lev,np.array([-0.24, -0.16, -0.08, 0.0, 0.08, 0.16, 0.24]),family='sans-serif', fontsize=15)

cbax3 = fig.add_axes([0.085, 0.1, 0.02, 0.365])
cbar3 = plt.colorbar(plot3, cax=cbax3, orientation='vertical', ticks=np.arange(-50,51,10))
cbar3.ax.set_ylabel('C export (g C m$^{-2}$ yr$^{-1}$)', family='sans-serif', fontsize=18)
cbar3.ax.yaxis.set_ticks_position('left')
cbar3.ax.yaxis.set_label_position('left')
plt.yticks(np.arange(-50,51,10),np.arange(-50,51,10),family='sans-serif', fontsize=15)

cbax5 = fig.add_axes([0.92, 0.1, 0.02, 0.365])
cbar5 = plt.colorbar(plot5, cax=cbax5, orientation='vertical',ticks=np.arange(-10,11,2))
cbar5.ax.set_ylabel('respired C (Kg m$^{-2}$)', family='sans-serif', fontsize=18)
plt.yticks(np.arange(-10,11,2),np.arange(-10,11,2),family='sans-serif', fontsize=15)


#%%

fig.savefig('figures_for_publication/fig1.pdf', dpi=300, bbox_inches='tight')
fig.savefig('figures_for_publication/fig1.png', dpi=300, bbox_inches='tight')
