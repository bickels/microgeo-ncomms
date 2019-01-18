# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 18:47:57 2018

@author: bickels
"""
import numpy as np
import scipy.constants as const

#%% SOIL
def climatic_watercontent(tau, pet, awc, d=1):
    awcFC = 0.5*awc
    alpha = pet/(awcFC*d)
    wc = awcFC*np.exp(-alpha*tau)
    return wc

def sample_length(bld, msmp=1.):
    Vsoil = msmp/bld 
    L = Vsoil**(1./3)
    return L

def textural_length(snd,slt,cly):
    #https://www.tandfonline.com/doi/full/10.1080/03650340903099676?scroll=top&needAccess=true
    dx = 1e-3*np.exp((snd*np.log(1.025)+slt*np.log(0.026)+cly*np.log(0.001)))
    return dx

#%% PTF
def theta_ptf(OC,Cl,Si,BD,depth):
    '''Toth et al 2015 https://onlinelibrary.wiley.com/doi/pdf/10.1111/ejss.12192'''
    T = depth <= 0.3
    theta_S = 0.6819-0.06480*(1./(OC+1))-0.11900*BD**2-0.02668*T+0.001489*Cl+0.0008031*Si+0.02321*(1./(OC+1))*BD**2+0.01908*BD**2*T-0.0011090*Cl*T-0.00002315*Si*Cl-0.0001197*Si*BD**2-0.0001068*Cl*BD**2
    return theta_S

#%% NPP
def miamiNPP(MAT,MAP,a,b,c,d):
    '''https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2006GB002705'''
#    a = 0.5
#    b = 3.75
#    c = 6.34
#    d = 7.5
#    e = -0.13
    nppT = a*(1.+np.exp(b-c*MAT))**-1
    nppP = a*(1.-np.exp(d*MAP))
    return np.min(np.stack([nppT,nppP]),axis=0)

#%% MICROBIAL
def muTs(T, T0, TL, TH, dH0, dHL, dHH):
    '''schoolfield'''
    T,T0,TH,TL = const.convert_temperature([T,T0,TH,TL],'C','K')
    return ((T/T0)*np.exp((dH0/const.R)*(1./T0-1./T)))/(1+np.exp((dHL/const.R)*(1./TL-1./T))+np.exp((dHH/const.R)*(1./TH-1./T)))

def satf(x,v,K):
    return v*x/(K+x)

def cells_from_npp(npp, mat, m, fC, muTpars=(25,9.85,41.5,-5.43e3,-141.1e3,687.9e3), d=1.):
    mt = m*muTs(mat,*muTpars)
    Ncells = npp/(mt*fC*d)
    return Ncells

#%% percolation
def ns(s,se,tau, log=False):
    if log:
        return -tau*np.log(s)+s/se
    else:
#        return mp_pow(s,-tau)*mp_exp(-s/se)
        return s**-tau*np.exp(-s/se)

def Pp_logistic(p, pc, k, e=1., sm=1.):
#    return p/mp_pow(1.+mp_exp(-sm*k*(p-pc)),e)
    return p/(1.+np.exp(-sm*k*(p-pc)))**e

def scalepar(L,a,b):
    return L/(a+b*L)
#    return a*L**b

def Pplo(p, pc, L, ks=[0.39325394,  0.01192326],sm=1, e=1):#, es=[1.48357842, 0.21059964],ps=[0.86668178, 0.01764669]):
#    pcL = scalepar(L,*ps)
    k = scalepar(L,*ks)
#    e = scalepar(L,*es)
    return Pp_logistic(p, pc, k, sm=sm,e=e)

def aqueous_clusters(p, por, L, d, dx, s, sm=1., pc=False):
    L0 = int(L/dx)
    N0 = int(L**d/dx**d)
    
    if d==2:
        if pc:
            pc *= por
        else:
            # site percolation square grid
            pc = 0.592746*por
        v = 4/3.
        beta = 5/36.
        sigma = 36/91.
        gamma = 43/18.
        df = 91/48.
        tau = (d/df)+1.
    
    if d==3:
        if pc:
            pc *= por
        else:
            # site percolation simple cubic
            pc = 0.31160*por
        v = 0.87
        beta = 0.4181
        sigma = 0.445
        gamma = 1.8
        df = 2.52
        tau = (d/df)+1.
    
    se = abs(p-pc)**(-1/sigma)

    Pp = Pplo(p, pc, L0, sm=sm)

    nswc = ns(s[:,None], se, tau, log=False)
    
    ns0 = (p-Pp)/np.sum(s[:,None]*nswc,axis=0)
    
    Ncs = ns0*s[:,None]*nswc*N0
    NcsP = Pp*N0
    bNcsP = (NcsP > Ncs[-1]) #& bapc
    bsc = NcsP > s[-1]
#    bNcsP = NcsP < np.finfo('float').resolution*N0
    try:
        sc = np.tile(s[:,None],len(p))
        Ncs[-1,bNcsP] = NcsP[bNcsP]
        sc[-1,bsc] = NcsP[bsc]
        ncl = Ncs/sc
    except:
        sc = s.copy()
        if bNcsP:
            Ncs[-1] = NcsP
        if bsc:
            sc[-1] = NcsP
        ncl = Ncs[:,0]/sc
#        print 'single value'
    return sc, Ncs, ncl

def qdiversity(p,q,n=1,exp=True):
    p /= np.sum(n*p)
    if q == 1:
        if exp:
#            return mp_exp(-np.sum(n*p*mp_log(p)))
            return np.exp(-np.sum(n*p*np.log(p)))
        else:
#            return -np.sum(n*p*mp_log(p))
            return -np.sum(n*p*np.log(p))
    else:
#        return np.sum(n*mp_pow(p,q))**(1./(1.-q))
        return np.sum(n*p**q)**(1./(1.-q))

def habitat_diversity(sc, Ncs, Cd, dx, d, q, ncells=1):
    bCd = dx**d*Ncs*Cd < ncells
    Cns = dx**d*Ncs*Cd
    Cns[bCd] = np.nan
    Dq = np.array([qdiversity(Cn[~np.isnan(Cn)]/np.sum(Cn[~np.isnan(Cn)]), q) for Cn in Cns.T])
    return Dq

def species_diversity(sc, ncl, Cd, dx, d, q, lod=1e-6, fnsp=1., legacy=False):
#    divCd = (dx**d*Cd)[:,None]*sc[::-1]*np.cumsum((1./fnsp)*ncl[::-1],axis=0)/2.
#    bCdiv = (dx**d*Cd)[:,None]*sc[::-1] < lod
#    divCd = (dx**d*Cd)[:,None]*sc*np.cumsum((1./fnsp)*ncl[::-1],axis=0)/2.
    Ncells = (dx**d*Cd)[:,None]*sc
    bCdiv =  Ncells < lod
    Ncells[bCdiv] = np.nan
    if legacy:
        divCd = Ncells*fnsp*ncl#np.cumsum(ncl,axis=0)/2.
        divCd[bCdiv] = np.nan
        Dq = np.array([qdiversity(Cn[~np.isnan(Cn)]/np.sum(Cn[~np.isnan(Cn)]), q) for Cn in divCd])
        return divCd, Dq
    else:
        Dq = np.array([qdiversity(Cn[~np.isnan(Cn)], q, n=(fnsp*ncl)[~np.isnan(Cn)]) for Cn in Ncells])
#        if q==0:
#            Dq = np.array([fnsp*np.sum(ncl[~bC]) for bC in bCdiv])
#        else:
#            Dq = np.array([qdiversity(Cn[~np.isnan(Cn)]/np.sum(Cn[~np.isnan(Cn)]), q, n=fnsp*ncl[~np.isnan(Cn)], exp=False) for Cn in Ncells])
        return Dq

def mSAD(sc, ncl, Cd, dx, d, lod=1e-6, fnsp=1.):
#    Ncells = (dx**d*Cd)[:,None]*sc[::-1]
#    sp = fnsp*np.cumsum(ncl[::-1],axis=0)/2.
#    if len(sp.shape) == 1:
#        sad = np.array([np.interp(np.arange(1,nsp), sp, Nc) for Nc in Ncells])
#    else:
#        sad = np.array([np.interp(np.arange(1,nsp), spp, Nc) for spp,Nc in zip(sp,Ncells)])
    Ncells = (dx**d*Cd)[:,None]*sc
    bCdiv =  Ncells < lod
    Ncells[bCdiv] = np.nan
    spabu = Ncells
    sprnk = np.cumsum(fnsp*ncl,axis=0)/2.
    spncl = fnsp*ncl
    return sprnk, spabu, spncl


def Dq(p, por, L, dx, smax, sm, Cd, lod, fnsp, q):
    d = 3
    s = np.arange(1.,smax)
    sc, Ncs, ncl = aqueous_clusters(p, por, L, d, dx, s, sm)
    Dq = np.array([species_diversity(scc, ncll, np.array(Cd), dx, d, q, lod, fnsp) for scc, ncll in zip(sc.T,ncl.T)])
    return Dq[:,0]