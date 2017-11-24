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

NBBANDS=6
band_to_number={'u':0,'g':1,'r':2,'i':3,'z':4,'y4':5}
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
          'figure.figsize': (12, 8),
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
        plt.title("Atmosphere sim",weight="bold")
        plt.xlabel("$\lambda$ (Angstrom)",weight="bold")
        plt.ylabel("transmission",weight="bold")
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
        plt.title("all transmissions",weight="bold")
        plt.xlabel( '$\lambda$ (Angstrom)',weight="bold")
        plt.ylabel('transmission',weight="bold")
        plt.grid()
            
#------------------------------------------------------------------------------------        

       
class LSSTObservation(object):
    '''
    class LSSTObservation(object)
    
    Compute the product of the effective LSST transmission (already filter x atmosphere)
    by the SED
    
    '''
    def __init__(self,name):
        self.name = name        # name given to the instance
        self.NBBANDS = 0        # number of bands
        self.NBEVENTS = 0       # number of different atmosphere called events
        self.NBSED = 0          # number of different input SED
        self.obsarray = []       # container for the product of atm x filters x SED
        self.obssamplarray = []  # sampled array for magnitude calculation
        self.all_sed = []        # must be a pysynphot source
        self.all_transmission = []   # must be a pysynphot passband
        self.counts = []         # number of counts
        self.magnitude = []     # instrumental magnitude
        
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
        plt.figure()   
        if (sednum>=0 and sednum <self.NBSED):
            theobservation=self.obsarray[sednum]
            for event in np.arange(self.NBEVENTS):
                all_bands=theobservation[event]
                ib=0
                for bp in all_bands:
                    plt.plot(bp.wave,bp.flux,color=filtercolor[ib])
                    ib+=1
            plt.title("all observations",weight="bold")
            plt.xlabel('$\lambda$ (Angstrom)',weight="bold")
            plt.ylabel('flux',weight="bold")
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
        plt.figure() 
        #selection of the SED    
        if (sednum>=0 and sednum <self.NBSED):
            theobservation=self.obssamplarray[sednum]
            #loop on event
            for event in np.arange(self.NBEVENTS):
                all_bands=theobservation[event]
                ib=0
                #loop on band
                for flux in all_bands:
                    plt.plot(WL,flux,color=filtercolor[ib])
                    ib+=1
            plt.title("all sampled observations",weight="bold")
            plt.xlabel( '$\lambda$ (Angstrom)',weight="bold")
            plt.ylabel('flux',weight="bold")
            plt.grid()
            plt.xlim(WLMIN,WLMAX)
            
    def compute_counts(self):
        if len(self.obssamplarray) ==0:
            self.make_samplobservations()
        # loop on SED  
        self.counts=[]
        for sed in self.obssamplarray:   
            # loop on each event of a sed 
            all_obs_counts = []
            for obs in sed:
                #loop on band
                all_band_counts = []
                for band in obs:
                    counts=CountRate(WL,band) # my personnal implementation of counts
                    all_band_counts.append(counts)
                all_obs_counts.append(all_band_counts) 
            self.counts.append(all_obs_counts) 
            self.counts=np.array(self.counts)  # at the end, the array is converted in numpy
        return self.counts
    
    def compute_magnitudeold(self):
        if len(self.counts) == 0:
            self.compute_counts()
                # loop on SED  
        self.magnitude=[]
        for sed in self.counts:   
            # loop on each event of a sed 
            all_obs_mag = []
            for obs in sed:
                #loop on band
                all_band_mag = []
                for band_counts in obs:
                    mag=-2.5*np.log10(band_counts)
                    all_band_mag.append(mag)
                all_obs_mag.append(all_band_mag) 
            self.magnitude.append(all_obs_mag) 
        return self.magnitude
    
    def compute_magnitude(self):
        if len(self.counts) == 0:
            self.compute_counts()
                # loop on SED  
        self.magnitude=-2.5*np.log10(self.counts)
        return self.magnitude   
    
    
    def plot_countsold(self,sednum):
        if len(self.counts) == 0:
            self.compute_counts()
              
        if (sednum>=0 and sednum <self.NBSED):
            thecounts=self.counts[sednum]
            plt.figure() 
            for event in np.arange(self.NBEVENTS):
                all_bands_counts=thecounts[event]
                ib=0
                #loop on band
                for bandcounts in all_bands_counts:
                    plt.plot([event],[bandcounts],'o',color=filtercolor[ib])
                    ib+=1
            plt.title("all counts",weight="bold")
            plt.xlabel( 'event number',weight="bold")
            plt.ylabel('counts',weight="bold")
            plt.grid()
            
    def plot_counts(self,sednum):
        if len(self.counts) == 0:
            self.compute_counts()
              
        if (sednum>=0 and sednum <self.NBSED):
            plt.figure() 
            for ib in np.arange(NBBANDS):
                plt.plot(self.counts[sednum,:,ib],'-',color=filtercolor[ib])
            plt.title("all counts",weight="bold")
            plt.xlabel( 'event number',weight="bold")
            plt.ylabel('counts',weight="bold")
            plt.grid()
            
    def plot_magnitudesold(self,sednum):
        if len(self.magnitude) == 0:
            self.compute_magnitude()
              
        if (sednum>=0 and sednum <self.NBSED):
            themagnitudes=self.magnitude[sednum]
        
            for event in np.arange(self.NBEVENTS):
                all_bands_mag=themagnitudes[event]
                ib=0
                #loop on band
                for bandmag in all_bands_mag:
                    plt.plot([event],[bandmag],'o',color=filtercolor[ib])
                    ib+=1
            plt.title("all instrumental magnitudes",weight="bold")
            plt.xlabel( 'event number',weight="bold")
            plt.ylabel('magnitude',weight="bold")
            plt.grid()
            
    def plot_magnitudes(self,sednum):
        if len(self.magnitude) == 0:
            self.compute_magnitude()
              
        if (sednum>=0 and sednum <self.NBSED):
            plt.figure()
            for ib in np.arange(NBBANDS):
                plt.plot(self.magntitudes[sednum,:,ib],'-',color=filtercolor[ib])
            plt.title("all instrumental magnitudes",weight="bold")
            plt.xlabel( 'event number',weight="bold")
            plt.ylabel('magnitude',weight="bold")
            plt.grid()
            
    def get_magnitudeforfilternum(self,sednum,filternum):
        if len(self.magnitude) == 0:
            self.compute_magnitude()
            
        if (sednum>=0 and sednum <self.NBSED): 
            return self.magnitude[sednum,:,filternum]
        else:
            return None
        
    def get_magnitudeforfiltername(self,sednum,filtername):
        if len(self.magnitude) == 0:
            self.compute_magnitude()
            
        filternum=band_to_number[filtername]
        if (sednum>=0 and sednum <self.NBSED):        
            return self.get_magnitudeforfilternum(sednum,filternum)
        else:
            return None
            
 
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
    
    