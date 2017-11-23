#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 15:00:16 2017

@author: dagoret
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysynphot as S
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import astropy.units as u

filtercolor=['blue','green','red','orange','grey','black']
WLMIN=3000.
WLMAX=11000.

NBINS=10000
BinWidth=(WLMAX-WLMIN)/float(NBINS)
WL=np.linspace(WLMIN,WLMAX,NBINS)

LSST_COLL_SURF=35*(u.m)**2/(u.cm)**2  # LSST collectif surface
S.refs.setref(area=LSST_COLL_SURF.decompose(), waveset=None)
S.refs.set_default_waveset(minwave=3000, maxwave=11000, num=8000, delta=1, log=False)
S.refs.showref()

# to enlarge the sizes
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (6, 4),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)

#---------------------------------------------------------------------------------
def CountRate(wl,fl):
    dlambda=BinWidth 
    df=wl*fl*LSST_COLL_SURF/(S.units.C*S.units.H)*dlambda
    # (erg/s/cm2/Angstrom) x (Angstrom)  x  (cm^2) / (erg . s . Angstrom /s) * Angstrom
    # units :  s-1
    count=df.sum()
    return count
#---------------------------------------------------------------------------------
    
#---------------------------------------------------------------------------------
def InstrumMag(wl,fl):
    m=-2.5*np.log10(CountRate(wl,fl))
    return m
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
def ComputeColor(wl,fl1,fl2):
    m1=InstrumMag(wl,fl1)
    m2=InstrumMag(wl,fl2)
    return m1-m2
#---------------------------------------------------------------------------------
    
#-------------------------------------------------------------------------------------
class Atmosphere(object):
    '''
    Class Atmosphere(object):
        Container for several atmospheric simulation.
        It is supposed all simulation have the same wavelengths.
        Wavelength are only known at run-time.
        Atmospheric transmission are stored in a numpy 2D array
    '''
    def __init__(self,name):
        '''
        Atmosphere.__init__:
            Initialization of containers
        '''
        self.atmcodename = name  # name for the atm simulation
        self.nbevents=0           # number of atm transm entries 
        self.NBWL=0              # number of Wavelength bins
        self.array = np.array([])# numpy array container 
        self.pys_pb = []         # pysynphot container of transmissions 
        
    def fill_array(self,newarray):
        '''
        function to enter atmospheric simulation one by one.
        - first entry : the wavelength array
        - next entries : the atmospheric transmissions
        '''
        
        nbwl=len(newarray)       # get the size of the array
        
        if self.NBWL==0:         # enter the wavelength array
            self.NBWL=nbwl       # number of wavelength bins
            self.array=np.array([]).reshape(0,self.NBWL)
            
        if self.NBWL==nbwl:      #check the size of the input array
            self.array = np.r_[self.array, [newarray]]
            self.nbevents+=1
        else:
            print 'error in the number of wavelength bins' 
       
    def print_array(self):  
        print self.array
        print type(self.array)
        
    def get_array(self):
        return self.array
    
    def make_pys_pb(self):
        '''
        function to transform the atmospheric sim data into a pysynphot list
        of transmission
        '''
        self.pys_pb = []         # overwrite the container 
        wl_atm=self.array[0,:]   # extract the wavelength array 
        nbevents=self.nbevents
        for eventnum in np.arange(1,nbevents):         
            atmname='ATM_'+str(eventnum)
            tr_atm=self.array[eventnum,:]
            bp_atm= S.ArrayBandpass(wl_atm*10.,tr_atm, name=atmname)
            self.pys_pb.append(bp_atm)
            
        return self.pys_pb
    
    def get_pys_pb(self):
        if len(self.pys_pb)== 0:
             self.make_pys_pb()
        return self.pys_pb
     
    def plot_pys_bp(self):
        if len(self.pys_pb) == 0:
            self.make_pys_pb()
            
        plt.figure()
        for bp in self.pys_pb:
            plt.plot(bp.wave,bp.throughput)
        plt.grid()
        plt.title("Atmosphere sim")
        plt.xlabel("$\lambda$ (Angstrom)")
        plt.ylabel("transmission")
        plt.savefig("atm-transm.png")
        
        
#---------------------------------------------------------------------------        



        
#------------------------------------------------------------------------------------
        
    
class LSSTTransmission(object):
    '''
    class LSSTTransmission(object)
    
    Compute the product of atmospheric transmissions with filter bands
    Fill a 2D array of transmission
    
    '''
    def __init__(self,name):
        self.name = name
        self.NBBANDS = 0
        self.NBEVENTS = 0
        self.array = []          # container for the product of atm x filters
        self.det_pb= []
        self.atm_pb= []
        
        
    def fill_det_allbands(self,all_bands):
        self.det_pb= all_bands
        self.NBBANDS=len(all_bands)
       
    def fill_atm_allevents(self,all_atms):
        self.atm_pb= all_atms
        self.NBEVENTS=len(all_atms)
        
    def make_transmissions(self):
        if len(self.array) != 0:
            return self.array
        for event in np.arange(self.NBEVENTS):
            transmission_event=[]
            for band in np.arange(self.NBBANDS):               
                transmission_event.append(self.atm_pb[event]*self.det_pb[band])
            self.array.append(transmission_event) # add row by row
        return self.array
                
    def get_transmissions(self):
        if len(self.array) == 0:
            self.make_transmissions()
        return self.array
    
    def plot_transmissions(self):
        if(len(self.array))== 0:
            self.make_transmissions()
            

        for event in np.arange(self.NBEVENTS):
            all_bands=self.array[event]
            ib=0
            for bp in all_bands:
                plt.plot(bp.wave,bp.throughput,color=filtercolor[ib])
                ib+=1
        plt.title("all transmissions")
        plt.xlabel( '$\lambda$ (Angstrom)')
        plt.ylabel('transmission')
        plt.grid()
            
#------------------------------------------------------------------------------------        

       
class LSSTObservation(object):
    '''
    class LSSTObservation(object)
    
    Compute the product of the effective LSST transmission (already filter x atmosphere)
    by the SED
    
    '''
    def __init__(self,name):
        self.name = name
        self.NBBANDS = 0
        self.NBEVENTS = 0
        self.NBSED = 0
        self.obsarray = []       # container for the product of atm x filters x SED
        self.obssamplarray = []  # sampled array for magnitude calculation
        self.all_sed = []        # must be a pysynphot source
        self.all_transmission = []   # must be a pysynphot passband
        
        
    def fill_sed(self,all_sed):
        self.all_sed= all_sed
        self.NBSED=len(all_sed)
        
    def fill_transmission(self,all_transm):
        self.all_transmission=all_transm
        self.NBEVENTS=len(all_transm)
        self.NBBANDS=len(all_transm[0])
       
    def make_observations(self):
        if len(self.obsarray)!=0:
            return self.obsarray
        # loop on all SED
        self.obsarray=[]
        for sed in self.all_sed:
            # loop on atmospheric events
            all_obs_persed=[]
            for transmission in self.all_transmission:
                all_obs_perevent= []
                #loop on all bands
                for band in transmission:
                    obs= S.Observation(sed,band)   # do OBS = SED x Transmission
                    all_obs_perevent.append(obs)
                all_obs_persed.append(all_obs_perevent)
            self.obsarray.append(all_obs_persed)
        return self.obsarray
    
    def get_observations(self):
        if len(self.obsarray) ==0:
            return self.make_observations()
        else:
            return self.obsarray
        
    def plot_observations(self,sednum):
        if len(self.obsarray) ==0:
            self.make_observations()
            
        if (sednum>=0 and sednum <self.NBSED):
            theobservation=self.obsarray[sednum]
            for event in np.arange(self.NBEVENTS):
                all_bands=theobservation[event]
                ib=0
                for bp in all_bands:
                    plt.plot(bp.wave,bp.flux,color=filtercolor[ib])
                    ib+=1
            plt.title("all observations")
            plt.xlabel( '$\lambda$ (Angstrom)')
            plt.ylabel('flux')
            plt.grid()
            plt.xlim(WLMIN,WLMAX)
            
    def make_samplobservations(self):
        '''
        Resample the flux with equi-width bins
        '''
        if len(self.obssamplarray)!=0:
            return self.obssamplarray
        
        self.obssamplarray=[]                  # output contaioner
        for sedsource in self.obsarray:        # loop on input SED sources
            # loop on atmospheric events
            all_obssampl_persedsource=[]
            for obs_per_event in sedsource:    # loop on all event for that sed             
                all_obssampl_bands= []
                #loop on all bands
                for obsband in obs_per_event:  # loop on bands
                    # do interpolation
                    func=interp1d(obsband.wave,obsband.flux,kind='cubic')
                    flux=func(WL)
                    all_obssampl_bands.append(flux) # save the band
                all_obssampl_persedsource.append(all_obssampl_bands) # save that event   
            self.obssamplarray.append(all_obssampl_persedsource) # save all the event for that sed
        return self.obssamplarray
    
    def plot_samplobservations(self,sednum):
        if len(self.obssamplarray) ==0:
            self.make_samplobservations()
            
        if (sednum>=0 and sednum <self.NBSED):
            theobservation=self.obssamplarray[sednum]
            for event in np.arange(self.NBEVENTS):
                all_bands=theobservation[event]
                ib=0
                for flux in all_bands:
                    plt.plot(WL,flux,color=filtercolor[ib])
                    ib+=1
            plt.title("all sampled observations")
            plt.xlabel( '$\lambda$ (Angstrom)')
            plt.ylabel('flux')
            plt.grid()
            plt.xlim(WLMIN,WLMAX)
 
#-----------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    print 'hello'
    
    NBEVENTS=20
    NBCOLUMNS=5
    
    MyAtm=Atmosphere('libradtran')
    
    
    for event in np.arange(NBEVENTS):
        newarray=np.random.random_sample(NBCOLUMNS)
        MyAtm.fill_array(newarray)
 
    
    
    MyAtm.print_array()
    
    theatm=MyAtm.get_array()
    
    print theatm
    
    