# -*- coding: utf-8 -*-
"""
Created on Thu Nov  15 16:17:44 2018

@author: bickels
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import lognorm

#%% CUSTOM FUNCTIONS
from functions.model import cells_from_npp, theta_ptf, climatic_watercontent, sample_length, textural_length, qdiversity, aqueous_clusters, mSAD
from functions.misc import rank, constants

#%% CONSTANTS
fC = constants('fC') #carbon content per cell
mt = constants('mt') #maintenance rate
muTpars = constants('muTpars') #parametrization of temperature dependency (Schoolfield)
lsh,lsc = constants('zpars') #parametrization of vertical distribution in soil profile

#%% DATA
covariates = pd.read_csv(r'./data/EMP_snippet.csv')
relative_abundance = pd.read_csv(r'./data/EMP_SAD_snippet.csv')

# perform calculations of bacterial diversity for first sample (.iloc[0])
ID, lat, lon, depth, BLD, SND, CLY, SLT, MAT, PET, DRY, NPP, SOC = covariates.iloc[0]
sample_relative_abundance = relative_abundance.select_dtypes('number').iloc[0]

#%% CALCULATIONS
# estimate cell density/carrying capacity based on 25% of NPP
fz = lognorm.pdf(depth, lsh, loc=0, scale=lsc)
cell_density = cells_from_npp(0.25*NPP, MAT, mt, fC, muTpars)

# estimate climatic water content
por = theta_ptf(CLY*1e2, SLT*1e2, BLD*1e-6, SOC*1e-3, depth)
cwc = climatic_watercontent(DRY, 1e-3*PET, por)

# estimate sample length
L = sample_length(BLD, msmp=0.175)

# estimate characteristic textural length
dx = textural_length(SND, SLT, CLY)

# calculate diversity of sample
D0 = qdiversity(sample_relative_abundance, q=0)
D1 = qdiversity(sample_relative_abundance, q=1)
evenness = D1/D0

#%% PREDICTION
# calculate size and number of aqueous habitats
s = np.arange(1,len(sample_relative_abundance)+1,dtype='int64')
sc, _, ncl = aqueous_clusters(cwc, por, L, 3, dx, s, pc=0.52)

# calculate species abundance distribution
# one species per habitat
model_relative_abundance  = mSAD(sc, ncl, np.array([cell_density]), dx, 3, lod=1e4, fnsp=1)
model_relative_abundance /= np.sum(model_relative_abundance)

# multiple species per habitat
model_relative_abundance_multi  = mSAD(sc, ncl, np.array([cell_density]), dx, 3, lod=1e4, fnsp=s**(1/3.))
model_relative_abundance_multi /= np.sum(model_relative_abundance_multi)

#%% PLOT
print 'Sample richness={:.2f} and evenness={:.2f}'.format(D0,evenness)

plt.plot(rank(sample_relative_abundance), sample_relative_abundance,'^', label='observation')
plt.plot(rank(model_relative_abundance), model_relative_abundance,'.', label='one species')
plt.plot(rank(model_relative_abundance_multi), model_relative_abundance_multi,'.', label='multiple species')
plt.legend()
plt.title('Observed and modeled species abundance distribution')
plt.xlabel('species rank')
plt.ylabel('relative abundance')
plt.loglog()
plt.show()
