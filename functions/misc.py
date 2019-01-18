# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 18:47:57 2018

@author: bickels
"""
import numpy as np

#%%
def rank(v):
    """ returns ranks of values

    Parameters
    ----------
    v : list or ndarray
        values to be ranked

    Returns
    -------
    ndarray
        ranks of values (descending)
    """
    ranks = np.empty(len(v),dtype='int')
    ranks[np.argsort(-v)] = np.arange(1,len(v)+1)
    return ranks

#%% CONSTANTS
def constants(key,val=None,mode='r'):
	""" helper function to import/export constants to textfile as python dictionaries

    Parameters
    ----------
	key : str
        key for value

	val : value
        value to retrieve. None (default) or value to save if mode 'a'.
	
    mode: str
        either 'r'(default) for read or 'a' for append

    Returns
    -------
    value
        value in dictionary
    """
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

