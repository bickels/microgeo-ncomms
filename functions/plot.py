# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 18:47:57 2018

@author: bickels
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.interpolate import griddata
#FIXME: UGLY hack
import os
os.environ['PROJ_LIB'] = r'C:\Users\bickels\Anaconda\envs\mgeo\pkgs\proj4-4.9.3-vc9_5\Library\share'

from mpl_toolkits.basemap import Basemap#, maskoceans, shiftgrid
from mpl_toolkits.axes_grid import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

#%%
def show_map(lat, lon, a, projection='robin', points=None, cbar=True, bg=False, **kwargs):
#    ax = plt.gca()
#    ax.set_rasterization_zorder(0)
    if bg:
        global_map = Basemap(projection=projection, lat_0=0, lon_0=0, resolution=None,llcrnrlat=-90,urcrnrlat=90)
#        global_map.bluemarble()
        global_map.shadedrelief()
    else:
        global_map = Basemap(projection=projection, lat_0=0, lon_0=0, resolution='l',llcrnrlat=-90,urcrnrlat=90)
        global_map.drawcoastlines()
        global_map.drawmapboundary()
        x,y = global_map(*np.meshgrid(lon,lat))
        global_map.pcolormesh(x, y, a, **kwargs)

    global_map.drawmeridians(np.arange(-180, 180, 30),labels=[1,0,0,0],labelstyle='+/-',linewidth=0.1)
    global_map.drawparallels(np.arange(-90, 90, 30),labels=[1,0,0,0],labelstyle='+/-',linewidth=0.1)
    
    if cbar:
        if 'vmax' in kwargs and 'vmin' in kwargs:
            extend = 'both'
        elif 'vmin' in kwargs:
            extend = 'min'
        elif 'vmax' in kwargs:
            extend = 'max'
        else:
            extend = 'neither'
        global_map.colorbar(extend=extend, size='2%')
    
    if points:
        for (lon,lat),label in points:
            px,py = global_map(lon,lat)
            global_map.scatter(px,py,label=label,s=3,edgecolors='w',linewidths=0.2, zorder=10)
            plt.legend(fontsize='small')
    
    plt.tight_layout()

def jointplot(x, y, c='k', m='o', labels=('x','y','c'), limits=False, s=20, bins=(50,50), alpha=1., logx=False, logy=False,cmap=None, **kwargs):
    fig = plt.figure()#figsize=(8,8))
    gs = mpl.gridspec.GridSpec(2, 3, height_ratios=[0.5,5], width_ratios=[5,0.5,0.5])
    
    ax2 = fig.add_subplot(gs[1,0])
    if logx:
       ax2.set_xscale('log')
    if logy:
       ax2.set_yscale('log')
    
    xy = ax2.scatter(x, y, c=c, s=s, marker=m, alpha=alpha,lw=0,cmap=cmap, **kwargs)
    
    ax0 = fig.add_subplot(gs[0,0], sharex=ax2)
    ax0.axis('off')
    ax0.hist(x[~np.isnan(x)], color='k', normed=True, histtype='stepfilled',stacked=True, alpha=0.25, lw=0, bins=bins[0])
    
    ax3 = fig.add_subplot(gs[1,1], sharey=ax2)
    ax3.axis('off')
    ax3.hist(y[~np.isnan(y)], color='k', normed=True, histtype='stepfilled',stacked=True, alpha=0.25, lw=0, bins=bins[1],orientation='horizontal')
    
    plt.subplots_adjust(hspace=.0, wspace=.0)
    
    xlab,ylab,clab = labels
    ax2.set_xlabel(xlab)
    ax2.set_ylabel(ylab)
    
    if limits:
        xl,xh,yl,yh = limits
        ax2.set_xlim(xl,xh)
        ax2.set_ylim(yl,yh)
        
    cax = fig.add_subplot(gs[1,2])
    cax.axis('off')
    divider = make_axes_locatable(cax)
    cax = divider.append_axes('right', size='50%', pad=0.05)
    if cmap is not None:
        cb = fig.colorbar(xy, cax=cax, orientation='vertical')
        cb.set_label(clab)
        cb.set_alpha(1)
        cb.draw_all()

def imdist(a, labels=('x','y','z'), cmap='viridis', norm=None, vmax=None, vmin=None, extent=None):
    nx,ny=a.shape
    x=np.arange(0,nx)
    y=np.arange(0,ny)
    asp = nx/float(ny)
    
    fig = plt.figure(figsize=(8*asp,8))
    gs = mpl.gridspec.GridSpec(2, 2, height_ratios=[0.5,2], width_ratios=[2,asp*2])
        
    ax2 = fig.add_subplot(gs[2])

    xy = ax2.imshow(a.T, cmap=cmap, norm=norm, vmax=vmax, vmin=vmin, extent=None, interpolation=None, origin='lower')
    plt.xlim(0,nx)
    plt.ylim(0,ny)
    
    ax0 = fig.add_subplot(gs[0], sharex=ax2)
    ax0.axis('off')
    ax0.plot(x,np.nansum(a,axis=1), color='k')
    
    ax3 = fig.add_subplot(gs[3], sharey=ax2)
    ax3.axis('off')
    ax3.plot(np.nansum(a,axis=0), y, color='k')
    
    plt.subplots_adjust(hspace=.0, wspace=.0)
    
    xlab,ylab,clab = labels
    ax2.set_xlabel(xlab)
    ax2.set_ylabel(ylab)
    
    cax = fig.add_subplot(gs[1])
    cax.axis('off')
    divider = make_axes_locatable(cax)
    cax = divider.append_axes('right', size='50%', pad=0.05)
    
    if cmap is not None:
        cb = fig.colorbar(xy, cax=cax, orientation='vertical')
        cb.set_label(clab)
        cb.set_alpha(1)
        cb.draw_all()

def pshow(land0, v, cmap='viridis', norm=None, interpolation='None', extent=None, **kwargs):
    a = np.full_like(land0, np.nan, dtype='float')
    a[land0] = v
    plt.imshow(np.roll(a,a.shape[1]/2),cmap=cmap,norm=norm, extent=extent, **kwargs)

def plot_clust(labeled,cmap='nipy_spectral',**kwargs):    
    b = np.arange(labeled.max() + 1)
    np.random.shuffle(b)
    shuffled_lab = b[labeled]
    fg = labeled != 0
    plt.imshow(shuffled_lab*fg, cmap=cmap, interpolation='None',**kwargs)

def add_scalebar(ax, length, text, location=4, **kwargs):
    scalebar = AnchoredSizeBar(ax.transData, length, text, location, **kwargs)
    ax.add_artist(scalebar)

def kdecolor(x,y):
    xy = np.vstack([x,y])
    xy[np.isnan(xy)] = 0.
    xy[~np.isfinite(xy)] = 0.
    return gaussian_kde(xy)(xy)

def gradient_fill_between(x0, y0, y1, color='k', min_alpha=0.05, max_alpha=0.95, npts=1024):
    x = np.r_[x0,x0[::-1]]
    y = np.r_[y0,y1[::-1]]
    t = np.r_[np.full(len(y0),min_alpha),np.full(len(y1),max_alpha)]
    xx=np.linspace(x.min(),x.max(),npts)
    yy=np.linspace(y.min(),y.max(),npts)
    X,Y = np.meshgrid(xx/x.max(),yy/y.max())
    T = griddata((x/x.max(),y/y.max()),t,(X,Y),method='linear')
    path = mpl.path.Path(np.array([x,y]).transpose())
    patch = mpl.patches.PathPatch(path, facecolor='None', edgecolor='None')
    ax = plt.gca()
    ax.add_patch(patch)
    ax.set_rasterization_zorder(0)
    z=np.empty(T.shape+(4,))
    rgb = mpl.colors.colorConverter.to_rgb(color)
    z[:,:,:3] = rgb
    z[:,:,-1] = T
    plt.imshow(z,interpolation='bilinear', 
               origin='lower',
               extent=[x.min(),x.max(),y.min(),y.max()],
               aspect='auto', 
               clip_path=patch, 
               clip_on=True,
               zorder=-1)
