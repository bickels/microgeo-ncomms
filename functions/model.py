# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 18:47:57 2018

@author: bickels
"""
import numpy as np
import scipy.constants as const

#%% SOIL
def climatic_watercontent(tau, pet, awc, d=1):
    """ estimation of climatic soil water content assuming soil reaches field capacity after rainfall event and dries during consecutive dry days (tau) at a constant rate.

    Parameters
    ----------
    tau : ndarray
        time between rainfall events

    pet : ndarray
        potential evapotranspiration

    awc : ndarray
        available water holding capacity

    d : ndarray
        soil depth (default 1m)

    Returns
    -------
    ndarray
        climatic soil water content
    """
    awcFC = 0.5*awc
    alpha = pet/(awcFC*d)
    wc = awcFC*np.exp(-alpha*tau)
    return wc

def sample_length(bld, msmp=1.):
    """ linear length of soil sample.

    Parameters
    ----------
    bld : ndarray
        bulk density

    msmp : ndarray
        mass of the soil sample (default 1g)
    
    Returns
    -------
    ndarray
        linear length of soil sample
    """
    Vsoil = msmp/bld 
    L = Vsoil**(1./3)
    return L

def textural_length(snd,slt,cly):
    """ Soil characteristic textural length.

    Parameters
    ----------
    snd : ndarray
        sand content

    slt : ndarray
        silt content

    cly : ndarray
        clay content

    Returns
    -------
    ndarray
        textural length scale

    Source
    ------     
    https://www.tandfonline.com/doi/full/10.1080/03650340903099676?scroll=top&needAccess=true
    https://dl.sciencesocieties.org/publications/sssaj/pdfs/48/1/SS0480010142
    """
    dx = 1e-3*np.exp((snd*np.log(1.025)+slt*np.log(0.026)+cly*np.log(0.001)))
    return dx

#%% PTF
def theta_ptf(cly,slt,bld,soc,depth):
    """ Pedotransferfunction for estimating porosity (available water holding capacity).

    Parameters
    ----------
    cly : ndarray
        clay content

    slt : ndarray
        silt content

    bld : ndarray
        bulk density

    soc : ndarray
        soil organic carbon content

    Returns
    -------
    ndarray
        porosity

    Source
    ------
    https://onlinelibrary.wiley.com/doi/pdf/10.1111/ejss.12192
    """
    T = depth <= 0.3
    theta_S = 0.6819-0.06480*(1./(soc+1))-0.11900*bld**2-0.02668*T+0.001489*cly+0.0008031*slt+0.02321*(1./(soc+1))*bld**2+0.01908*bld**2*T-0.0011090*cly*T-0.00002315*slt*cly-0.0001197*slt*bld**2-0.0001068*cly*bld**2
    return theta_S

#%% MICROBIAL
def muTs(T, T0, TL, TH, dH0, dHL, dHH):
    """Schoolfield model of temperature dependency of reaction rates.

    Parameters
    ----------
    T : ndarray
        temperature

    T0 : float or int
        reference temperature

    TL : float or int
         low temperature associated with inactivation 

    TH : float or int
         high temperature associated with inactivation 

    dH0 : float or int
        activation enthalpy of reaction

    dHL : float or int
        change in enthalpy associated with low temperature inactivation

    dHH : float or int
        change in enthalpy associated with high temperature inactivation

    Returns
    -------
    ndarray
        factor for reduction of reaction rate fT

    Source
    ------
    http://linkinghub.elsevier.com/retrieve/pii/0022519381902460
    """
    T,T0,TH,TL = const.convert_temperature([T,T0,TH,TL],'C','K')
    return ((T/T0)*np.exp((dH0/const.R)*(1./T0-1./T)))/(1+np.exp((dHL/const.R)*(1./TL-1./T))+np.exp((dHH/const.R)*(1./TH-1./T)))

def cells_from_npp(npp, mat, m, fC, muTpars=(25, 9.85, 41.5, -5.43e3, -141.1e3, 687.9e3), d=1.):
    """ Conversion of net primary productivity to cell density/carrying capacity considering temperature dependent maintenance rates.
 
    Parameters
    ----------
    npp : ndarray
        net primary productivity (or fraction thereof)

    mat : ndarray
        mean annual temperature

    m : float or int
         maintenance rate

    fC : float or int
         carbon content per cell

    muTpars : tuple
        Schoolfield model parameters

    d : float or int
        soil depth (default 1m)

    Returns
    -------
    ndarray
        cell density at carrying capacity
    """
    mt = m*muTs(mat,*muTpars)
    Ncells = npp/(mt*fC*d)
    return Ncells

#%% percolation
def ns(s,se,tau):
    """ Size distribution of percolation clusters.

    Parameters
    ----------
    s : ndarray
        size of cluster

    se : ndarray
        size of characteristic cluster

    tau : float
         scaling exponent
    
    Returns
    -------
    ndarray
        cluster numbers
    
    Source
    ------
    http://www.sciencedirect.com/science/article/pii/0370157379900607
    """    
    return s**-tau*np.exp(-s/se)

def Pp_logistic(p, pc, k):
    """ Bounded general logistic function to approximate percolation transition.

    Parameters
    ----------
    p : ndarray
        occupancy probability

    pc : float
        critical point

    k : float
        steepness at inflection point (smoothness of transition)
    
    Returns
    -------
    ndarray
        bounded logistic function
    """
    return p/(1.+np.exp(-k*(p-pc)))

def aqueous_clusters(p, por, L, d, dx, s, k=1., pc=False):
    """ Number of aqueous habitats (percolation clusters) as function of occupancy probability/water contents based on percolation theory.

    Parameters
    ----------
    p : ndarray
        occupancy probability

    por : float or int
        porosity/void fraction

    L : float or int
        linear length of soil sample

    d : float or int
        dimensionality, 2 or 3
    
    dx : float or int
        characteristic textural length
    
    s : ndarray
        cluster size

    k : float or int
        smoothness of percolation transition
    
    pc : float or int
        critical point
    
    Returns
    -------
    sc : ndarray
        sizes of clusters

    Ncs : ndarray
        number of sites occupied by cluster of size s

    ncl : ndarray
        number of clusters of size s
    """
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

    Pp = Pp_logistic(p, pc, k)

    nswc = ns(s[:,None], se, tau)
    
    ns0 = (p-Pp)/np.sum(s[:,None]*nswc,axis=0)
    
    Ncs = ns0*s[:,None]*nswc*N0
    NcsP = Pp*N0
    bNcsP = (NcsP > Ncs[-1])
    bsc = NcsP > s[-1]

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
    return sc, Ncs, ncl

def qdiversity(p,q,n=1):
    """ Generalized diversity of order q.

    Parameters
    ----------
    p : ndarray
        probability

    q : float or int
        order of diversity
    
    n : float or ndarray
        weights

    Returns
    -------
    ndarray 
        diversity of order q (q^D)

    Source:
    -------
    http://doi.wiley.com/10.2307/1934352
    """
    p /= np.sum(n*p)
    if q == 1:
        return np.exp(-np.sum(n*p*np.log(p)))
    else:
        return np.sum(n*p**q)**(1./(1.-q))


def mSAD(sc, ncl, Cd, dx, d, lod, fnsp=1.):
    """ Species abundance distribution.

    Parameters
    ----------
    sc : ndarray
        sizes of clusters/habitats
    
    ncl : ndarray
        number of clusters of size s

    Cd : ndarray
        cell density
    
    dx : ndarray
        characteristic textural length
    
    d : int or float
        dimensionality

    lod : int or float
        limit of detection
    
    fnsp : ndarray or int or float
        number of species in a single aqueous habitat
    
    Returns
    -------
    spa : ndarray
        species abundance
    """
    Ncells = (dx**d*Cd)[:,None]*sc
    bCdiv =  Ncells < lod
    Ncells[bCdiv] = np.nan
    spabu = Ncells
    spncl = fnsp*ncl
    spabu[np.isnan(spabu)] = 0
    spa = (spabu*spncl).T
    spa[-1] = 0. #remove largest cluster
    return spa[:,0]
