# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 18:47:57 2018

@author: bickels
"""
import numpy as np
import scipy.ndimage as ndi

#%%
def rle(a):
    """ run length encoding. Partial credit to R rle function.
        Multi datatype arrays catered for including non Numpy
        returns: tuple (runlengths, startpositions, values) """
    ia = np.array(a)                  # force numpy
    n = len(ia)
    if n == 0:
        return (None, None, None)
    else:
        y = np.array(ia[1:] != ia[:-1])     # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)   # must include last element posi
        z = np.diff(np.append(-1, i))       # run lengths
        p = np.cumsum(np.append(0, z))[:-1] # positions
        return (z, p, ia[i])

def lin_gapfill(x):
    try:
        ind = np.arange(len(x))
        not_nan = ~np.isnan(x)
        return np.interp(ind,ind[not_nan], x[not_nan])
    except:
        return x

def sum_chunk(x, chunk_size, axis=-1):
    shape = x.shape
    if axis < 0:
        axis += x.ndim
    shape = shape[:axis] + (-1, chunk_size) + shape[axis+1:]
    x = x.reshape(shape)
    return x.sum(axis=axis+1)

def moving_mean(x,y,n, dropna=False):
    xy = np.array([x,y])
    if dropna:
        xy = xy[:,xy[0].argsort()]
        bna = np.isnan(xy).any(axis=0)
        xy = xy[:,~bna]
    else:
        xy = xy[:,xy[0].argsort()]
    return xy[0], np.convolve(xy[1],np.ones(n)/n,mode='same')

def rank(v):
    ranks = np.empty(len(v),dtype='int')
    ranks[np.argsort(-v)] = np.arange(1,len(v)+1)
    return ranks

def radial_mask(a, b, nx, ny, r=1):
    xm, ym = np.ogrid[-a:nx-a, -b:ny-b]
    return xm*xm + ym*ym <= r*r

def radial_mean(a, pos=(0,0),bins=20):
    sx, sy = a.shape
    ox,oy = pos
    X,Y = np.ogrid[0:sx,0:sy]
    r = np.hypot(X - ox, Y - oy)
    rbin = (bins* r/r.max()).astype(np.int)
    return ndi.mean(a, labels=rbin, index=np.arange(1, rbin.max() +1))

def rms(sim,obs):
    return np.mean((sim-obs)**2)**0.5

#%% FILE HANDLING
def query_spatial(data, desc, loc, out='colrow', dataout=None):
    loc = np.asarray(loc)
    height = desc['height']
    width = desc['width']
    af = desc['transform']
    locout = np.zeros_like(loc)
    if out == 'colrow':
        for i, l in enumerate(loc):
            locout[i] = ~af*l
    elif out == 'xy':
        for i, l in enumerate(loc):
            locout[i] = af*l
    col,row = locout.T
    bNA = np.isnan(col) | np.isnan(row) | (row>height) | (col>width)
    if dataout == None:
        dataout = np.zeros_like(bNA,dtype=float)
    dataout[~bNA] = data[row[~bNA].astype(int),col[~bNA].astype(int)]
    return dataout, bNA

#def read_soilgrids(f,loc):
#    with rio.open(f) as X:
#        rows,cols = np.array([X.index(lo-360,lt) for lo,lt in loc]).T
#        x = X.read(1)[rows,cols]
#    x = x.astype('float')
#    x[x==X.meta['nodata']] = np.nan
#    return x

#%% CONSTANTS
def constants(key,val=None,mode='r'):
    with open('CONSTANTS.txt', mode) as f:
        if val:
            f.write('\ndict('+key+'='+str(val)+')')
            print r'saved: {'+key+':'+str(val)+'}'
        else:
            dicts = [eval(line) for line in f]
            d0 = {}
            for d1 in dicts:
                d0 = dict(d0, **d1)
            constant = d0
            return constant[key]

