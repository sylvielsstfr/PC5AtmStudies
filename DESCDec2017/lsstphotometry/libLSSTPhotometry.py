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

import matplotlib.cm as cm
import matplotlib as mpl
cmap = cm.jet

import pysynphot as S
from scipy.interpolate import interp1d
import astropy.units as u

NBBANDS=6
band_to_number={'u':0,'g':1,'r':2,'i':3,'z':4,'y4':5}
number_to_band={0:'u',1:'g',2:'r',3:'i',4:'z',5:'y4'}
filtercolor=['blue','green','red','orange','grey','black']
NBCOLORS=NBBANDS-1
number_to_color={0:'U-G',1:'G-R',2:'R-I',3:'I-Z',4:'Z-Y'}
color_to_number={'U-G':0,'G-R':1,'R-I':2,'I-Z':3,'Z-Y':4}
mpl_colors_col=['b','g','r','y','k']

WLMIN=3000. # Minimum wavelength : PySynPhot works with Angstrom
WLMAX=11000. # Minimum wavelength : PySynPhot works with Angstrom

NBINS=10000 # Number of bins between WLMIN and WLMAX
BinWidth=(WLMAX-WLMIN)/float(NBINS) # Bin width in Angstrom
WL=np.linspace(WLMIN,WLMAX,NBINS)   # Array of wavelength in Angstrom

LSST_COLL_SURF=35*(u.m)**2/(u.cm)**2  # LSST collectif surface
S.refs.setref(area=LSST_COLL_SURF.decompose(), waveset=None)
S.refs.set_default_waveset(minwave=WLMIN, maxwave=WLMAX, num=NBINS, delta=BinWidth, log=False)
S.refs.showref()

EXPOSURE=30.0                      # LSST Exposure time

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
    '''
    CountRate(wl,fl):
        This Count rate is calculated in a BandWidth
        Input :
            wl : wavelength array
            fl : flux
        Output :
            count
            
        CountRate is normalized to the number of photoelectrons per sec.
        It takes into account the LSST collection surface and the unit conversion into photoel.
        Notice one does not need ADU units.
    '''
    dlambda=BinWidth 
    df=wl*fl*LSST_COLL_SURF/(S.units.C*S.units.H)*dlambda
    # (erg/s/cm2/Angstrom) x (Angstrom)  x  (cm^2) / (erg . s . Angstrom /s) * Angstrom
    # units :  s-1
    count=df.sum()
    return count
#---------------------------------------------------------------------------------
    
#---------------------------------------------------------------------------------
def InstrumMag(countrate,dt=EXPOSURE):
    '''
    InstrumMag(countrate,dt=EXPOSURE):
        Compute Instrumental Magnitude given the Exposure time
        Notice that exposure time is taken into account in magnitudes, not before !
        countrate is a numpy array.
        Both magnitude and magnitude error are numpy arrays.
    '''

    m=np.where(countrate>0,-2.5*np.log10(countrate*dt),0 )        # the magnitude
    dm=np.where(countrate>0,-2.5/2.3/np.sqrt(countrate*dt),0)     # the error on magnitude 
    
    inst_mag=np.array([m,dm])              # return both magntitude and its error.
    return inst_mag
#---------------------------------------------------------------------------------



#---------------------------------------------------------------------------------
def ComputeColor(m1_e,m2_e):
    '''
    ComputeColor(m1,m2):
        Compute Color and its error
    '''
    C=m1_e[0]-m2_e[0] # Color
    dC=np.sqrt(m1_e[1]*m1_e[1]+m2_e[1]*m2_e[1]) # error on color
   
    return np.array([C,dC])
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
            plt.plot(bp.wave,bp.throughput,lw=2)
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
                plt.plot(bp.wave,bp.throughput,color=filtercolor[ib],lw=2)
                ib+=1
        plt.title("all transmissions",weight="bold")
        plt.xlabel( '$\lambda$ (Angstrom)',weight="bold")
        plt.ylabel('transmission',weight="bold")
        plt.grid()
            
#------------------------------------------------------------------------------------        

#------------------------------------------------------------------------------------------------       
class LSSTObservation(object):
    '''
    class LSSTObservation(object)
    
    Compute the product of the effective LSST transmission (already filter x atmosphere)
    by the SED
    
        self.name = name        # name given to the instance
        self.NBBANDS = 0        # number of bands
        self.NBEVENTS = 0       # number of different atmosphere called events
        self.NBSED = 0          # number of different input SED
        self.obsarray = []       # container for the product of atm x filters x SED
        self.obssamplarray = []  # sampled array for magnitude calculation
        self.all_sed = []        # must be a pysynphot source
        self.all_transmission = []   # must be a pysynphot passband
        self.counts = []         # number of counts
        self.magnitude = []     # instrumental magnitude for each SED, each atmosphere, each band
        self.magnit_zeropt = [] # zero point
        self.magnit_zeropt_err = []  # zero point error
        self.colors = []              # colors U-G, G-R, R-I, I-Z, Z-Y
        self.colors_err = []          # color errors
    
    
    '''
    def __init__(self,name):
        '''
        Initialize variable containers
        '''
        self.name = name        # name given to the instance
        self.NBBANDS = 0        # number of bands
        self.NBEVENTS = 0       # number of different atmosphere called events
        self.NBSED = 0          # number of different input SED
        self.obsarray = []       # container for the product of atm x filters x SED
        self.obssamplarray = []  # sampled array for magnitude calculation
        self.all_sed = []        # must be a pysynphot source
        self.all_transmission = []   # must be a pysynphot passband
        self.counts = []         # number of counts rate, per second
        self.magnitude = []      # instrumental magnitude for each SED, each atmosphere, each band
        self.magnitude_err= []   # error on Magnitude
        self.magnit_zeropt = []  # zero point
        self.magnit_zeropt_err = []  # zero point error
        self.magnit_bias =  []        # magnitude bias wrt zero point
        self.magnit_bias_err =  []        # magnitude bias wrt zero point
        self.colors = []              # colors U-G, G-R, R-I, I-Z, Z-Y
        self.colors_err = []          # color errors
        self.colors_bias = []              # colors bias (-0pt) U-G, G-R, R-I, I-Z, Z-Y
        self.colors_bias_err = []          # color bias errors (-0pt)
    #----------------------------------------------------------------------------        
    def get_NBSED(self):
        '''
        get_NBSED() getter for the number of SED
        
        '''
        return self.NBSED
    #----------------------------------------------------------------------------        
    def fill_sed(self,all_sed):
        '''
        fill_sed(all_sed) :
            Initialization
            fill all SED array of all objects among which the zero point will be evaluated
        '''
        self.all_sed= all_sed
        self.NBSED=len(all_sed)
   #----------------------------------------------------------------------------            
    def fill_transmission(self,all_transm):
        '''
        fill_transmission(all_transm) :
            Initialization
            fill transmission for each filter (including atmosphere and LSST throughput)
            Usually transmission varies with one parameter at a time, say PWV or VAOD, or z_am
        '''
        self.all_transmission=all_transm
        self.NBEVENTS=len(all_transm)
        self.NBBANDS=len(all_transm[0])
    #----------------------------------------------------------------------------       
    def make_observations(self):
        '''
        make_observations():
            1) First stage of calculation, compute the flux in photoelectrons
            ----------------------------------------------------------------
            For each SED in the SED collector-
               For each atmospheric transmission conditions and various airmasses
                   For each of the six LSST filter band
                       Compute of each wl in the SED the corresponding Flux
        '''
        if len(self.obsarray)!=0:
            return self.obsarray
        # loop on all SED
        self.obsarray=[]
        # loop on all kind of SED of the catalog
        for sed in self.all_sed:
            # loop on atmospheric events
            
            sed.convert('flam') # to be sure every spectrum is in flam unit
                                #----------------------------------------------
            
            all_obs_persed=[]
            # for each SED, loop on all stransmissions
            for transmission in self.all_transmission:
                all_obs_perevent= []
                #loop on all bands of LSST
                for band in transmission:
                    # force=[extrap|taper] <---  check if OK
                    # Normalement pysynphot se débrouille avec les unités de la SED
                    obs= S.Observation(sed,band,force='extrap')   # do OBS = SED x Transmission
                    all_obs_perevent.append(obs)
                all_obs_persed.append(all_obs_perevent)
            self.obsarray.append(all_obs_persed)
        return self.obsarray
    #----------------------------------------------------------------------------  
    def get_observations(self):
        '''
        get_observations():
            Getter of all observations
        '''
        if len(self.obsarray) ==0:
            return self.make_observations()
        else:
            return self.obsarray
    #----------------------------------------------------------------------------
    def get_observationsforSED(self,sednum):
        '''
            get_observationsforSED(sednum):
                Get all observations for a single SED by giving the sed number
        '''
        if len(self.obsarray) ==0:
            self.make_observations()
    
        if (sednum>=0 and sednum <self.NBSED):
            return self.obsarray[sednum]
        else:
            print 'get_observationsforSED :: bad SED number',sednum
            return None
    #----------------------------------------------------------------------------        
    def plot_observations(self,sednum):
        '''
            plot_observations(sednum):
                control plot all observations for a single SED selected by sednum
        '''
        if len(self.obsarray) ==0:
            print 'plot_observations :: len(self.obsarray) = ',len(self.obsarray)
            print ' plot_observations :: ==> self.make_observations()'
            self.make_observations()
        plt.figure()   
        if (sednum>=0 and sednum <self.NBSED):
            theobservation=self.obsarray[sednum]
            for event in np.arange(self.NBEVENTS):
                all_bands=theobservation[event]
                ib=0
                for bp in all_bands:
                    plt.plot(bp.wave,bp.flux,color=filtercolor[ib],lw=2)
                    ib+=1
            plt.title("all observations",weight="bold")
            plt.xlabel('$\lambda$ (Angstrom)',weight="bold")
            plt.ylabel('flux',weight="bold")
            plt.grid()
            plt.xlim(WLMIN,WLMAX)
    #----------------------------------------------------------------------------            
    def make_samplobservations(self):
        '''
        make_samplobservations():
              2nd Stage : Resample the observed flux with equi-width bins
              --------------------------------------------------------------
              This is nesserary to compute later magnitude by integrating the experimental flux
              over the wavelength range.
              The original flux which was filled acording the SED binning is replaced
              by an interpolated  flux in a regular wl binning.
              This is mandatory to calculate magnitudes later.
              Be sure make_observations() has been called previously
        '''
        if len(self.obssamplarray)!=0:        # if already computed return it
            return self.obssamplarray
        
        if len(self.obsarray)==0:
            self.make_observations()
        
        
        self.obssamplarray=[]                  # init output contaioner
        for sedsource in self.obsarray:        # loop on input SED sources
            # loop on atmospheric events
            all_obssampl_persedsource=[]
            for obs_per_event in sedsource:    # loop on all observation-event for that sed             
                all_obssampl_bands= []
                #loop on all bands U G R I Z Y
                for obsband in obs_per_event:  # loop on bands
                    # do interpolation inside each band
                    #------------------------------------
                    func=interp1d(obsband.wave,obsband.flux,kind='cubic')
                    flux=func(WL)              # flux in each bin in array WL
                    all_obssampl_bands.append(flux) # save the band
                all_obssampl_persedsource.append(all_obssampl_bands) # save that event   
            self.obssamplarray.append(all_obssampl_persedsource) # save all the event for that sed
        return self.obssamplarray
    #----------------------------------------------------------------------------
    def get_samplobservationsforSED(self,sednum):
        '''
        get_samplobservationsforSED(sednum):
            Getter for sampled observations of a single SED having the number sednum.
            Check obssamplarray has been calculated.
        '''
        if len(self.obssamplarray) ==0:
            self.make_samplobservations()
    
        if (sednum>=0 and sednum <self.NBSED):
            return self.obssamplarray[sednum]
        else:
            print 'get_samplobservationsforSED :: bad SED number',sednum
            return None
    #---------------------------------------------------------------------------------    
    def plot_samplobservations(self,sednum):
        '''
        plot_samplobservations(sednum):
            Plot sampled observation for one SED of number sednum
            Plot Flux x  vs wl
            Check obssamplarray has been calculated
        '''
        if len(self.obssamplarray) ==0:
            print 'plot_samplobservations :: len(self.obssamplarray) = ',len(self.obssamplarray)
            print ' plot_samplobservations :: ==> self.make_samplobservations()'
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
                    plt.plot(WL,flux,color=filtercolor[ib],lw=2)
                    ib+=1
            plt.title("all sampled observations",weight="bold")
            plt.xlabel( '$\lambda$ (Angstrom)',weight="bold")
            plt.ylabel('flux',weight="bold")
            plt.grid()
            plt.xlim(WLMIN,WLMAX)
    #---------------------------------------------------------------------------------        
    def plot_samplobservationsflux(self,sednum):
        '''
        plot_samplobservationsflux(sednum):
            Plot Flux x wl vs wl
             
        '''
        if len(self.obssamplarray) ==0:
            print 'plot_samplobservations :: len(self.obssamplarray) = ',len(self.obssamplarray)
            print ' plot_samplobservations :: ==> self.make_samplobservations()'
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
                    plt.plot(WL,flux*WL,color=filtercolor[ib],lw=2)
                    ib+=1
            plt.title("all sampled observations flux ",weight="bold")
            plt.xlabel( '$\lambda$ (Angstrom)',weight="bold")
            plt.ylabel('flux * wl',weight="bold")
            plt.grid()
            plt.xlim(WLMIN,WLMAX)        
    #-------------------------------------------------------------------------------        
    def compute_counts(self):
        '''
            compute_counts(): 
            3rd Stage : Compute the photoelectron rate in each LSST band
            ------------------------------------------------------------
            External CountRate() function is used
            Check obssamplarray has been calculated
            Notice the very important fact that self.counts is converted into a 3D numpy array:
                - first index : sednum
                - second index : atm/ z_am condition
                - third index : obs band 0,... 5    
        '''
        if len(self.obssamplarray) ==0:
            print 'compute_counts :: len(self.obssamplarray) = ',len(self.obssamplarray)
            print 'compute_counts :: ==> self.make_samplobservations()'
            self.make_samplobservations()
            
            
        # loop on SED  
        self.counts=[]                            # init count
        for sed in self.obssamplarray:            # loop on SED 
            # loop on each event of a sed 
            all_obs_counts = []
            for obs in sed:                       # loop on all sampled obs of each SED
                #loop on bands U G R I Z Y
                all_band_counts = []              # init band count container
                for band in obs:                  # loop on each band of sampled obs of that SED
                    counts=CountRate(WL,band)     # my personnal implementation of counts
                    all_band_counts.append(counts)
                all_obs_counts.append(all_band_counts) 
            self.counts.append(all_obs_counts) 
        self.counts=np.array(self.counts)  # at the end, the array is converted in numpy array
        return self.counts
    #-------------------------------------------------------------------------------------   
    def compute_magnitude(self,dt=EXPOSURE):
        '''
        compute_magnitude() : 
            4th Stage :  Compute Magnitude and Magnitude error
            ----------------------------------------------------
            External InstrumMag() function is used.
            Notice Exposure must be provided at this stage.
            Check counts has been calculated.
            Notice the very important fact that self.magnitude and elf.magnitude_err are
            converted into a 3D numpy array:
                - first index : sednum
                - second index : atm/ z_am condition
                - third index : obs band 0,... 5
                
            Temporary: the error on magnitude is not returned.
        '''
        if len(self.counts) == 0:
            print 'compute_magnitude :: len(self.counts) = ',self.counts
            print 'compute_magnitude :: ==> self.counts()'
            self.compute_counts()
        self.magnitude, self.magnitude_err= InstrumMag(self.counts,dt) # do not consider here the error
        return self.magnitude  
    #-----------------------------------------------------------------------------------------
    # Notice calibration analysis start here (could be moved to another class)
    #------------------------------------------------------------------------------------------
    def compute_magnit_zeropt(self,dt=EXPOSURE):
        '''
        compute_magnit_zeropt(): 
            5th stage : Compute magnititude zero point and its error magnitude error over all SED
            ------------------------------------------------------------------------
            BECAREFULL : What if the observation are at different airmasses ?
            MAYBE we should separate different airmasses in different objects ?
            THIS POINT OF DIFFERENT AIRMASSES MUST BE STUDIED !!!!!!!!!!!
            Check if magntiudes has been computed.
        '''
        if len(self.magnitude) == 0:
            print 'magnit_zeropt :: len(self.magnitude) = ',self.magnitude
            print 'magnit_zeropt :: ==> self.counts()'
            self.compute_magnitude(dt)
            
            
        # CHOOSE AVERAGE OR MEDIAN OVER ALL SED for each atm conditions and   
        self.magnit_zeropt=np.average(self.magnitude,axis=0) # zero point average over all SED
                                                             # axis=0 means for all SED
        
        self.magnit_zeropt_err=np.sqrt(np.average(self.magnitude_err*self.magnitude_err,axis=0)) #zero point average error
        return self.magnit_zeropt
    #----------------------------------------------------------------------------------------------                 
    def plot_counts(self,sednum):
        '''
        plot_counts(sednum):
            Plot the counts rate for a given SED
        '''
        if len(self.counts) == 0:
            self.compute_counts()
              
        if (sednum>=0 and sednum <self.NBSED):
            plt.figure() 
            for ib in np.arange(NBBANDS):
                plt.plot(self.counts[sednum,:,ib],'-',color=filtercolor[ib],lw=2)
            plt.title("all counts",weight="bold")
            plt.xlabel( 'event number',weight="bold")
            plt.ylabel('counts',weight="bold")
            plt.grid()
    #---------------------------------------------------------------------------------------     
    def plot_magnitudes(self,sednum,dt=EXPOSURE):
        '''
        plot_magnitudes(sednum,dt=EXPOSURE):
            Plot the magnitudes
            
        '''
        if len(self.magnitude) == 0:                # check magnitudes has been computed
            self.compute_magnitude(dt)
              
        if (sednum>=0 and sednum <self.NBSED):      # check SED number
            plt.figure()
            for ib in np.arange(NBBANDS):
                plt.plot(self.magnitude[sednum,:,ib],'-',color=filtercolor[ib],lw=2)
            plt.title("all instrumental magnitudes",weight="bold")
            plt.xlabel( 'event number',weight="bold")
            plt.ylabel('magnitude',weight="bold")
            plt.grid()
        else:
            print " plot_magnitudes() : bad SED number",sednum
   #-----------------------------------------------------------------------------------------         
    def plot_magnit_zeropt(self,sednum,dt=EXPOSURE):
        '''
        plot_magnit_zeropt(sednum,dt=EXPOSURE):
            Plot the zero-point magnitudes.
            
        '''

        if len(self.magnit_zeropt) == 0:              # compute zero-point
            self.compute_magnit_zeropt(dt)

        if (sednum>=0 and sednum <self.NBSED):        # check zero-point
            plt.figure()
            for ib in np.arange(NBBANDS):
                delta_mag=self.magnitude[sednum,:,ib]-self.magnit_zeropt[:,ib]
                plt.plot(delta_mag,'-',color=filtercolor[ib],lw=2)
            plt.title("all instrumental Dleta mag (zeropt - substraction) ",weight="bold")
            plt.xlabel( 'event number',weight="bold")
            plt.ylabel('magnitude',weight="bold")
            plt.grid()
        else:
            print " plot__magnit_zeropt() : bad SED number",sednum
   #---------------------------------------------------------------------------------------------     
            
    def get_magnitudeforfilternum(self,sednum,filternum,dt=EXPOSURE):
        '''
        get_magnitudeforfilternum(sednum,filternum,dt=EXPOSURE):
            getter for magnitude
        '''
        if len(self.magnitude) == 0:
            self.compute_magnitude(dt)
            
        if (sednum>=0 and sednum <self.NBSED and filternum>0 and filternum< self.NBBANDS): 
            return self.magnitude[sednum,:,filternum]
        else:
            print "get_magnitudeforfilternum(sednum,filternum) bad number argument :: "
            print "sednum = ",sednum
            print "filternum =", filternum
            return None
    #-----------------------------------------------------------------------------------------
    def get_magnitudeforfiltername(self,sednum,filtername,dt=EXPOSURE):
        '''
        get_magnitudeforfiltername(self,sednum,filtername,dt=EXPOSURE):
            getter for magnitude
        '''
        if len(self.magnitude) == 0:
            self.compute_magnitude(dt)
            
        filternum=band_to_number[filtername]
        if (sednum>=0 and sednum <self.NBSED):        
            return self.get_magnitudeforfilternum(sednum,filternum)
        else:
            print "get_magnitudeforfiltername(sednum,filtername) bad number argument :: "
            print "sednum = ",sednum
            print "filtername =", filtername  
            return None
    #--------------------------------------------------------------------------------------------
        
    def get_magnitzeroptforfilternum(self,sednum,filternum,dt=EXPOSURE):
        '''
        get_magnitzeroptforfilternum(sednum,filternum,dt=EXPOSURE):
            getter for magnitude zero-point
        '''
 
        if len(self.magnit_zeropt) == 0:
            self.compute_magnit_zeropt(dt)
            
        if (sednum>=0 and sednum <self.NBSED and filternum>0 and filternum< self.NBBANDS): 
            return self.magnitude[sednum,:,filternum]-self.magnit_zeropt[:,filternum]
        else:
            print "get_magnitzeroptforfilternum(sednum,filternum) bad number argument :: "
            print "sednum = ",sednum
            print "filternum =", filternum
            
            return None
    #-------------------------------------------------------------------------------------------       
    def get_magnitzeroptforfiltername(self,sednum,filtername,dt=EXPOSURE):
        '''
        get_magnitzeroptforfiltername(self,sednum,filtername,dt=EXPOSURE):
            getter for magnitude zero-point
        '''
       
        if len(self.magnit_zeropt) == 0:
            self.compute_magnit_zeropt(dt)
            
        filternum=band_to_number[filtername]
        
        if (sednum>=0 and sednum <self.NBSED):        
            return self.get_magnitzeroptforfilternum(sednum,filternum)
        else:
            print "get_magnitzeroptforfiltername(sednum,filtername) bad number argument :: "
            print "sednum = ",sednum
            print "filtername =", filtername  
            return None
    #---------------------------------------------------------------------------------------      
    def compute_magnitude_bias(self,dt=EXPOSURE):
        '''
        compute_magnitude_bias(self,dt=EXPOSURE):
            6th Stage : Compute the magnitude bias
            ----------------------------------------      
            Compute magnitude bias with respect to zero point for each atmosphere event
            Notice no error computed by now.
            
            In the end, the magnitude_bias array is transformed in numpy array of 3 indexes
        
        '''
  
        if len(self.magnit_zeropt) == 0:                # check zero-point has been computed
            self.compute_magnit_zeropt(dt)
            
        self.magnit_bias =  []                          # init bias
        self.magnit_bias_err =  [] 
        
        # LOOP ON all SED
        for sednum in  np.arange(self.NBSED): 
            # loop on each event of a sed 
            all_magnitudes_all_events_thatsed=self.magnitude[sednum]
            all_event_mag_biased = [] 
            
            # LOOP  on atmospheric events magnitudes
            NbEVT=-1
            for mag_event in all_magnitudes_all_events_thatsed :
                NbEVT+=1
                all_band_mag_biased = []
                #LOOP on bands U G R I Z Y
                for iband in np.arange(NBBANDS):
                    # subtract the magnitude to that of the zero point
                    #----------------------------------------------------------
                    thebias=mag_event[iband]-self.magnit_zeropt[NbEVT,iband]
                    
                    all_band_mag_biased.append(thebias)
                all_event_mag_biased.append(all_band_mag_biased) 
            self.magnit_bias.append(all_event_mag_biased) 
        self.magnit_bias=np.array(self.magnit_bias)  # at the end, the array is converted in numpy array
        return self.magnit_bias
    #---------------------------------------------------------------------------------------------------------- 
        
    #---------------------------------------------------------------------------------------------------------- 
    def compute_colors(self,dt=EXPOSURE):
        '''
        compute_colors(self,dt=EXPOSURE):
            7th stage : Compute the colors :
            --------------------------------
            Compute the colors for all SED-Atm-Events
            
        '''
        if len(self.magnitude) == 0:                  # check if magnitudes exists
            print 'compute_colors :: len(self.magnitude) = ',self.magnitude
            print 'compute_color :: ==> self.magnitude()'
            self.compute_magnitude(dt)
            
        self.colors=[]            # initialisation of colors
        
        
        # LOOP on SED
        for sednum in  np.arange(self.NBSED): 
            all_magnitudes_all_events_thatsed=self.magnitude[sednum]
            all_event_colors = [] 
            # LOOP  on atmospheric events magntitues
            for mag_event in all_magnitudes_all_events_thatsed :
                all_band_colors = []
                # LOOP on BANDS U, G, R, I, Z, and not Y   
                for iband in np.arange(NBCOLORS):
                    # Subtract the magnitudes
                    thecolor=mag_event[iband]-mag_event[iband+1]
                    all_band_colors.append(thecolor)
                all_event_colors.append(all_band_colors) 
            self.colors.append(all_event_colors) 
        self.colors=np.array(self.colors)  # at the end, the array is converted in numpy array
        return self.colors
    #----------------------------------------------------------------------------------------------------------
    
    def get_colors(self):
        return self.colors
    
    #----------------------------------------------------------------------------------------------------------
    def show_colors_0pt(self,index0,xarray,title,xtitle,figname,dt=EXPOSURE):
            '''
            show_color_0pt(index0,xarray,title,xtitle,figname,dt=EXPOSURE):
            
            input:
                index0 : index of colors with the smallest wavelength
                xarray: array in X ex  VAOD, PWV, or even airmass
            output:
                Show the zero-pont color difference wrt the reference atm-am condition
                
            
            '''
            
            if len(self.magnit_zeropt) == 0:       # check if magnit_zeropt has been calculated
                self.compute_magnit_zeropt(dt)
            
            # LOOP ON colors
            for icol in np.arange(NBCOLORS):
                thelabel=number_to_color[icol]
                # zero-point color 
                col=self.magnit_zeropt[:,icol]- self.magnit_zeropt[:,icol+1] # all colors for all atm event and all colors
                # zero-point color for reference atmospheric-am conditions
                col0=self.magnit_zeropt[index0,icol]-self.magnit_zeropt[index0,icol+1]
                
                deltacol=col-col0 # color difference relative to reference condition
            
                plt.plot(xarray,deltacol,color=mpl_colors_col[icol],lw=2,label=thelabel)
                  
            plt.plot([xarray[0],xarray[-1]],[0.005,0.005],'k:',lw=3)
            plt.plot([xarray[0],xarray[-1]],[-0.005,-0.005],'k:',lw=3)
   
            plt.title(title,weight='bold')
            plt.xlabel(xtitle,weight='bold')
            plt.ylabel("$\Delta$ col (mag)",weight='bold')
            plt.legend(loc=2)
            plt.grid()
            plt.savefig(figname)   
    #----------------------------------------------------------------------------------------------------------
            
    
    
    
    
    #----------------------------------------------------------------------------------------------------------
    def show_colors(self,index0,xarray,title,xtitle,figname,dt=EXPOSURE):
            '''
            show_color(self,index0,xarray,title,xtitle,figname,dt=EXPOSURE)
            
            input:
                index0 : index of colors with the smallest wavelength
                xarray: array in X ex  VAOD, PWV or Airmass
            
            '''
            
            # loop on all sed
            if len(self.colors) == 0:
                print 'show_color :: len(self.colors) = ',self.colors
                print 'show_color :: ==> self.colors()'
                self.compute_colors(dt)
            
            # LOOP on SED
            for ised in np.arange(self.NBSED):
                col=self.colors[ised,:,:] # all colors for all atm event and all colors
                col0=self.colors[ised,index0,:]  # all colors for index0 atm event and all colors
                
                # LOOP on colors for plotting
                for icol in np.arange(NBCOLORS):
                    thelabel=number_to_color[icol]
                    deltacol=col[:,icol]-col0[icol]   # color difference relative to reference condition
                    if(ised==0):
                        plt.plot(xarray,deltacol,color=mpl_colors_col[icol],lw=1,label=thelabel)
                    else:
                        plt.plot(xarray,deltacol,color=mpl_colors_col[icol],lw=1)
            plt.plot([xarray[0],xarray[-1]],[0.005,0.005],'k:',lw=3)
            plt.plot([xarray[0],xarray[-1]],[-0.005,-0.005],'k:',lw=3)
   
            plt.title(title,weight='bold')
            plt.xlabel(xtitle,weight='bold')
            plt.ylabel("$\Delta$ col (mag)",weight='bold')
            plt.legend(loc=2)
            plt.grid()
            plt.savefig(figname)   
    #----------------------------------------------------------------------------------------------------------
            
    
    
    #----------------------------------------------------------------------------------------------------------
    def compute_colors_bias(self,dt=EXPOSURE):
        '''
        compute_colors_bias(self,dt=EXPOSURE):
            8th stage : Compute the colors bias :
            ---------------------------------------
            Compute biased colors
            The bias is always taken relative to the zero-point.
        
        
        '''
        # check if mangitude bias has been calculated
        if len(self.magnit_bias) == 0:
            print 'compute_colors_bias :: len(self.magnit_bias) = ',self.magnit_bias
            print 'compute_color_bias :: ==> self.compute_magnit_bias()'
            self.compute_magnit_bias(dt)
  

          
        self.colors_bias=[]                  # init color bias
        # LOOP on SED
        for sednum in  np.arange(self.NBSED): 
            all_magnit_bias_all_events_thatsed=self.magnit_bias[sednum] # select biased magnitudes
            all_event_colors_bias = [] 
            # LOOP on atmospheric-am events biased magnitudes
            for mag_event in all_magnit_bias_all_events_thatsed :
                all_band_colors_bias = []
                # loop on COLORS U,G,R?i?Z and not Y
                for iband in np.arange(NBCOLORS):
                     # Compute biased colors by magnitude subtraction one band wrt next one
                     #-----------------------------------------------------------------------------
                    thecolor_bias=mag_event[iband]-mag_event[iband+1]
                    all_band_colors_bias.append(thecolor_bias)
                all_event_colors_bias.append(all_band_colors_bias) 
            self.colors_bias.append(all_event_colors_bias) 
        self.colors_bias=np.array(self.colors_bias)  # at the end, the array is converted in numpy array
        return self.colors_bias
    #----------------------------------------------------------------------------------------------------------
    
    #----------------------------------------------------------------------------------------------------------
    def show_colors_bias(self,index0,xarray,title,xtitle,figname,dt=EXPOSURE):
            '''
            show_color_bias(self,index0,xarray,title,xtitle,figname,dt=EXPOSURE)
            
            input:
                index0 : index of colors with the smallest wavelength
                xarray: array in X ex  VAOD, PWV
            output:
                Show the colors wrt atm-am condition reference
                in X : the conditions
                in Y : Delta-Cor
            '''
            
            # Check of colors_bias has been computed
            if len(self.colors_bias) == 0:
                print 'show_color_bias :: len(self.colors) = ',self.colors_bias
                print 'show_color_bias :: ==> self.colors_bias()'
                self.compute_colors_bias(dt)
            
            # LOOP on all SED:
            for ised in np.arange(self.NBSED):
                col=self.colors_bias[ised,:,:] # all colors for all atm event and all colors
                col0=self.colors_bias[ised,index0,:]  # all colors for index0 of reference in  atm event 
                
                # LOOP on colors for plotting
                for icol in np.arange(NBCOLORS):
                    thelabel=number_to_color[icol]
                    deltacol=col[:,icol]-col0[icol] # Delta Col in Y
                    if(ised==0):
                        plt.plot(xarray,deltacol,color=mpl_colors_col[icol],lw=1,label=thelabel)
                    else:
                        plt.plot(xarray,deltacol,color=mpl_colors_col[icol],lw=1)
            plt.plot([xarray[0],xarray[-1]],[0.005,0.005],'k:',lw=3)
            plt.plot([xarray[0],xarray[-1]],[-0.005,-0.005],'k:',lw=3)
   
            plt.title(title,weight='bold')
            plt.xlabel(xtitle,weight='bold')
            plt.ylabel("$\Delta$ col (mag)",weight='bold')
            plt.legend(loc=2)
            plt.grid()
            plt.savefig(figname)   
    #----------------------------------------------------------------------------------------------------------
           

    #----------------------------------------------------------------------------------------------------------          
    def ShowColorBiasTrajectory(self,index0,Index_Start,xarray,zlabel,title,figname,dt=EXPOSURE):
        '''
        ShowColorBiasTrajectory(index0,Index_Start,xarray,xlabel,ylabel,title,figname,dt=EXPOSURE):
            input:
                index0: reference index on atm-am conditions
                Index_Start : index of colors with the smallest wavelength
                xarray: array in X ex  VAOD, PWV
            output:
                2D-plot of colors_bias , with each point corresponding to 1 SED-1-conditions.
                The colors of plots is linked to the condition.
                example 
                Plot of  G-R vs U-G points, one per SED-condition, the colors of points is related to cond
            
        
        
        '''
    
        # define the max color
        Index_Stop=Index_Start+3
    
        # reserve a table of ColorDiff in 3 indexes
        # - Index 1 : first or second color
        # - Index 2 : SED number
        # - Index 3 : each atmospheric-am condition
        ColorDiff=np.zeros([2,self.NBSED,len(xarray)])
        
        # the x-array indeed will be colors in this 2D plot
        x= xarray
        for i in range(len(x)):
            c = cmap(int(np.rint(x[i] / x.max() * 255)))

        norm = mpl.colors.Normalize(vmin=x.min(), vmax=x.max())

        fig = plt.figure(figsize=(10, 10))    
    
        # LOOP on SED to fill array ColorDiff[icol,ised,icond]
        for ised in np.arange(self.NBSED):
            # LOOP on the two selected colors for plotting along X and Y axis
            for icol in np.arange(Index_Start,Index_Stop-1):
                ColorIndex=icol-Index_Start
                              
                col=self.colors_bias[ised,:,:]        # get biased colors for all atm event and all colors
                col0=self.colors_bias[ised,index0,:]  # get biased colors for the index0 reference atm event 
                deltacol=col[:,icol]-col0[icol]       # compute the deltacol
                
                ColorDiff[ColorIndex,ised,:]=deltacol # fill the array

        # LOOP on SED to plot colors in the 2D plot of colors
        for ised in np.arange(self.NBSED):
            i=0  # index of atm-sed condition
            # LOOP on atm-am-conditions 
            for xel in xarray:
                c = cmap(int(np.rint(x[i] / x.max() * 255)))   # define the color of the point in 2D plot
                plt.plot(ColorDiff[0,ised,i],ColorDiff[1,ised,i],'o', markersize=20, mfc=c, mec=c)
                i+=1   # increase index on atm-sed condition 
            
            plt.plot([-0.005,0.005],[0.005,0.005],'k:',lw=2)
            plt.plot([-0.005,0.005],[-0.005,-0.005],'k:',lw=2)
            plt.plot([0.005,0.005],[-0.005,0.005],'k:',lw=2)
            plt.plot([-0.005,-0.005],[-0.005,0.005],'k:',lw=2)
    
        plt.grid()
        xtitle=number_to_color[Index_Start]+' color bias'
        ytitle=number_to_color[Index_Start+1]+' color bias'
        plt.xlabel(xtitle,fontsize=20,weight='bold')
        plt.ylabel(ytitle,fontsize=20,weight='bold')
        plt.title(title,fontsize=30,weight='bold')

        rect = 0.91,0.1,0.02,0.8   #  4-length sequence of [left, bottom, width, height] quantities.     
        ax1 = fig.add_axes(rect,label=zlabel)
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
        cb1.set_label(zlabel, rotation=90,fontsize=20,weight='bold')
        plt.savefig(figname)      

    #----------------------------------------------------------------------------------------------------------
 
    #----------------------------------------------------------------------------------------------------------          
    def ShowColorBiasTrajectoryVect(self,index0,Index_Start,xarray,zlabel,title,figname,dt=EXPOSURE):
        '''
        ShowColorBiasTrajectoryVect(index0,Index_Start,xarray,xlabel,ylabel,title,figname,dt=EXPOSURE):
            input:
                index0: reference index on atm-am conditions
                Index_Start : index of colors with the smallest wavelength
                xarray: array in X ex  VAOD, PWV
            output:
                2D-plot of colors_bias , with each point corresponding to 1 SED-1-conditions.
                The colors of plots is linked to the condition.
                example 
                Plot of  G-R vs U-G points, one per SED-condition, the colors of points is related to cond
            
        
        '''
    
        # define the max color
        Index_Stop=Index_Start+3
    
        # reserve a table of ColorDiff in 3 indexes
        # - Index 1 : first or second color
        # - Index 2 : SED number
        # - Index 3 : each atmospheric-am condition
        ColorDiff=np.zeros([2,self.NBSED,len(xarray)])
        
        # the x-array indeed will be colors in this 2D plot
        x= xarray
        for i in range(len(x)):
            c = cmap(int(np.rint(x[i] / x.max() * 255)))

        norm = mpl.colors.Normalize(vmin=x.min(), vmax=x.max())

        fig = plt.figure(figsize=(10, 10))    
    
        # LOOP on SED to fill array ColorDiff[icol,ised,icond]
        for ised in np.arange(self.NBSED):
            # LOOP on the two selected colors for plotting along X and Y axis
            for icol in np.arange(Index_Start,Index_Stop-1):
                ColorIndex=icol-Index_Start
                              
                col=self.colors_bias[ised,:,:]        # get biased colors for all atm event and all colors
                col0=self.colors_bias[ised,index0,:]  # get biased colors for the index0 reference atm event 
                deltacol=col[:,icol]-col0[icol]       # compute the deltacol
                
                ColorDiff[ColorIndex,ised,:]=deltacol # fill the array

        # LOOP on SED to plot colors in the 2D plot of colors
        for ised in np.arange(self.NBSED):
            
            # LOOP on atm-am-conditions 
            for i in np.arange(len(xarray)-1):
                c = cmap(int(np.rint(x[i] / x.max() * 255)))   # define the color of the point in 2D plot
                
                X=ColorDiff[0,ised,i]
                nextX=ColorDiff[0,ised,i+1]
                Y=ColorDiff[1,ised,i]
                nextY=ColorDiff[1,ised,i+1]
                U=nextX-X
                V=nextY-Y
                plt.quiver(X,Y, U, V,scale_units='xy', angles='xy', scale=1, color=c)
                
                #plt.plot(ColorDiff[0,ised,i],ColorDiff[1,ised,i],'o', markersize=20, mfc=c, mec=c)
            
            
            plt.plot([-0.005,0.005],[0.005,0.005],'k:',lw=2)
            plt.plot([-0.005,0.005],[-0.005,-0.005],'k:',lw=2)
            plt.plot([0.005,0.005],[-0.005,0.005],'k:',lw=2)
            plt.plot([-0.005,-0.005],[-0.005,0.005],'k:',lw=2)
    
        plt.grid()
        xtitle=number_to_color[Index_Start]+' color bias'
        ytitle=number_to_color[Index_Start+1]+' color bias'
        plt.xlabel(xtitle,fontsize=20,weight='bold')
        plt.ylabel(ytitle,fontsize=20,weight='bold')
        plt.title(title,fontsize=30,weight='bold')

        rect = 0.91,0.1,0.02,0.8   #  4-length sequence of [left, bottom, width, height] quantities.     
        ax1 = fig.add_axes(rect,label=zlabel)
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
        cb1.set_label(zlabel, rotation=90,fontsize=20,weight='bold')
        plt.savefig(figname)      

    #----------------------------------------------------------------------------------------------------------
    
    
    #----------------------------------------------------------------------------------------------------------
    def ShowColor2DPlot(self,index0,Index_Start,xarray,zlabel,title,figname,dt=EXPOSURE):
        '''
        ShowColor2DPlot(self,index0,Index_Start,xarray,zlabel,title,figname,dt=EXPOSURE):
            input:
                index0: reference index on atm-am conditions
                Index_Start : index of colors with the smallest wavelength
                xarray: array in X ex  VAOD, PWV
                zoom :
            output:
                2D-plot of colors, with each point corresponding to 1 SED-1-conditions.
                The colors of plots is linked to the condition.
                example 
                Plot of  G-R vs U-G points, one per SED-condition, the colors of points is related to cond
        
        '''
         # define the max color
        Index_Stop=Index_Start+3
    
        # reserve TWO tables of ColorDiff in 3 indexes
        # - Index 1 : first or second color
        # - Index 2 : SED number
        # - Index 3 : each atmospheric-am condition
        
        ColorDiff=np.zeros([2,self.NBSED,len(xarray)])
        ColorPosition=np.zeros([2,self.NBSED,len(xarray)])
    
         # the x-array indeed will be colors in this 2D plot
        x= xarray
        for i in range(len(x)):
            c = cmap(int(np.rint(x[i] / x.max() * 255)))

        norm = mpl.colors.Normalize(vmin=x.min(), vmax=x.max())

        fig = plt.figure(figsize=(10, 10))    
    
        # LOOP on SED
        for ised in np.arange(self.NBSED):
            # LOOP on the tow colors to fill the two Tables
            for icol in np.arange(Index_Start,Index_Stop-1):
            
                ColorIndex=icol-Index_Start
                
                
                col=self.colors_bias[ised,:,:]         # get biased colors for all atm event and all colors
                col0=self.colors_bias[ised,index0,:]   # get biased colors for the index0 reference atm event 
                deltacol=col[:,icol]-col0[icol]        # compute the deltacol
                ColorDiff[ColorIndex,ised,:]=deltacol           # exagerate the deltacol              
                ColorPosition[ColorIndex,ised,:]=self.colors[ised,:,icol] # original color of the source
         
        #LOOP on SED to plot colors in the 2D plot of colors
        for ised in np.arange(self.NBSED):
            i=0 # index of atm condition
            # LOOP on atm-am-conditions 
            for xel in xarray:
                c = cmap(int(np.rint(x[i] / x.max() * 255)))
                plt.plot(ColorPosition[0,ised,i],ColorPosition[1,ised,i],'o', markersize=5, mfc=c, mec=c)
                i+=1    # increase index on sed condition
            
    
        plt.grid()

        xtitle=number_to_color[Index_Start]+' color (mag)'
        ytitle=number_to_color[Index_Start+1]+' color (mag)'

        plt.xlabel(xtitle,fontsize=20,weight='bold')
        plt.ylabel(ytitle,fontsize=20,weight='bold')        
        plt.title(title,fontsize=30,weight='bold')
        #plt.xlim(-5,5)
        #plt.ylim(-5,5)


        rect = 0.91,0.1,0.02,0.8   #  4-length sequence of [left, bottom, width, height] quantities.     
        ax1 = fig.add_axes(rect,label=zlabel)
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='vertical')
        cb1.set_label(zlabel, rotation=90,fontsize=20,weight='bold')
        plt.savefig(figname)      
    #----------------------------------------------------------------------------------------------------------   
    
    
     #----------------------------------------------------------------------------------------------------------
    def ShowColor2DPlotOrigin(self,index0,Index_Start,xarray,zlabel,title,figname,dt=EXPOSURE):
        '''
        ShowColor2DPlotOrigin(self,index0,Index_Start,xarray,zlabel,title,figname,dt=EXPOSURE):
            input:
                index0: reference index on atm-am conditions
                Index_Start : index of colors with the smallest wavelength
                xarray: array in X ex  VAOD, PWV
                zoom :
            output:
                2D-plot of colors recentered, with each point corresponding to 1 SED-1-conditions.
                The colors of plots is linked to the condition.
                example 
                Plot of  G-R vs U-G points, one per SED-condition, the colors of points is related to cond
        
        '''
         # define the max color
        Index_Stop=Index_Start+3
    
        # reserve TWO tables of ColorDiff in 3 indexes
        # - Index 1 : first or second color
        # - Index 2 : SED number
        # - Index 3 : each atmospheric-am condition
        
        ColorDiff=np.zeros([2,self.NBSED,len(xarray)])
        ColorPosition=np.zeros([2,self.NBSED,len(xarray)])
        ColorShift=np.zeros([2,self.NBSED,len(xarray)])
    
         # the x-array indeed will be colors in this 2D plot
        x= xarray
        for i in range(len(x)):
            c = cmap(int(np.rint(x[i] / x.max() * 255)))

        norm = mpl.colors.Normalize(vmin=x.min(), vmax=x.max())

        fig = plt.figure(figsize=(10, 10))    
    
        # LOOP on SED
        for ised in np.arange(self.NBSED):
            # LOOP on the tow colors to fill the two Tables
            for icol in np.arange(Index_Start,Index_Stop-1):
            
                ColorIndex=icol-Index_Start
                
                
                col=self.colors_bias[ised,:,:]         # get biased colors for all atm event and all colors
                col0=self.colors_bias[ised,index0,:]   # get biased colors for the index0 reference atm event 
                deltacol=col[:,icol]-col0[icol]        # compute the deltacol
                ColorDiff[ColorIndex,ised,:]=deltacol           # exagerate the deltacol              
                ColorPosition[ColorIndex,ised,:]=self.colors[ised,:,icol] # original color of the source      
                # compute the relative shift
                ColorShift[ColorIndex,ised,:]= ColorPosition[ColorIndex,ised,:] - ColorPosition[ColorIndex,ised,index0]
        
        
        #LOOP on SED to plot colors in the 2D plot of colors
        for ised in np.arange(self.NBSED):
            i=0 # index of atm condition
            # LOOP on atm-am-conditions 
            for xel in xarray:
                c = cmap(int(np.rint(x[i] / x.max() * 255)))
                plt.plot(ColorShift[0,ised,i],ColorShift[1,ised,i],'o', markersize=5, mfc=c, mec=c)
                i+=1    # increase index on sed condition
            
    
        plt.grid()

        xtitle=number_to_color[Index_Start]+' relative color shift (mag)'
        ytitle=number_to_color[Index_Start+1]+' relative color shift (mag)'

        plt.xlabel(xtitle,fontsize=20,weight='bold')
        plt.ylabel(ytitle,fontsize=20,weight='bold')        
        plt.title(title,fontsize=25,weight='bold')
        #plt.xlim(-5,5)
        #plt.ylim(-5,5)


        rect = 0.91,0.1,0.02,0.8   #  4-length sequence of [left, bottom, width, height] quantities.     
        ax1 = fig.add_axes(rect,label=zlabel)
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='vertical')
        cb1.set_label(zlabel, rotation=90,fontsize=20,weight='bold')
        plt.savefig(figname)      
        #--------------------------------------
    
    
    #----------------------------------------------------------------------------------------------------------
    def ShowColorTrajectory(self,index0,Index_Start,xarray,zoom,zlabel,title,figname,dt=EXPOSURE):
        '''
        ShowColorTrajectory(self,index0,Index_Start,xarray,zlabel,title,figname,dt=EXPOSURE):
            input:
                index0: reference index on atm-am conditions
                Index_Start : index of colors with the smallest wavelength
                xarray: array in X ex  VAOD, PWV
                zoom :
            output:
                2D-plot of colors shifted by color-bias, with each point corresponding to 1 SED-1-conditions.
                The colors of plots is linked to the condition.
                example 
                Plot of  G-R vs U-G points, one per SED-condition, the colors of points is related to cond
        
        '''
         # define the max color
        Index_Stop=Index_Start+3
    
        # reserve TWO tables of ColorDiff in 3 indexes
        # - Index 1 : first or second color
        # - Index 2 : SED number
        # - Index 3 : each atmospheric-am condition
        
        ColorDiff=np.zeros([2,self.NBSED,len(xarray)])
        ColorDiff2=np.zeros([2,self.NBSED,len(xarray)])
    
         # the x-array indeed will be colors in this 2D plot
        x= xarray
        for i in range(len(x)):
            c = cmap(int(np.rint(x[i] / x.max() * 255)))

        norm = mpl.colors.Normalize(vmin=x.min(), vmax=x.max())

        fig = plt.figure(figsize=(10, 10))    
    
        # LOOP on SED
        for ised in np.arange(self.NBSED):
            # LOOP on the tow colors to fill the two Tables
            for icol in np.arange(Index_Start,Index_Stop-1):
            
                ColorIndex=icol-Index_Start
                
                
                col=self.colors_bias[ised,:,:]         # get biased colors for all atm event and all colors
                col0=self.colors_bias[ised,index0,:]   # get biased colors for the index0 reference atm event 
                deltacol=col[:,icol]-col0[icol]        # compute the deltacol
                ColorDiff[ColorIndex,ised,:]=deltacol*zoom             # exagerate the deltacol              
                ColorDiff2[ColorIndex,ised,:]=self.colors[ised,:,icol] # original color of the source
         
        #LOOP on SED to plot colors in the 2D plot of colors
        for ised in np.arange(self.NBSED):
            i=0 # index of atm condition
            # LOOP on atm-am-conditions 
            for xel in xarray:
                c = cmap(int(np.rint(x[i] / x.max() * 255)))
                plt.plot(ColorDiff[0,ised,i]+ColorDiff2[0,ised,0],ColorDiff[1,ised,i]+ColorDiff2[1,ised,0],'o', markersize=5, mfc=c, mec=c)
                i+=1    # increase index on sed condition
            
    
        plt.grid()

        xtitle=number_to_color[Index_Start]+' color (mag)'
        ytitle=number_to_color[Index_Start+1]+' color (mag)'

        plt.xlabel(xtitle,fontsize=20,weight='bold')
        plt.ylabel(ytitle,fontsize=20,weight='bold')        
        plt.title(title,fontsize=25,weight='bold')
        #plt.xlim(-5,5)
        #plt.ylim(-5,5)


        rect = 0.91,0.1,0.02,0.8   #  4-length sequence of [left, bottom, width, height] quantities.     
        ax1 = fig.add_axes(rect,label=zlabel)
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='vertical')
        cb1.set_label(zlabel, rotation=90,fontsize=20,weight='bold')
        plt.savefig(figname)      
        #----------------------------------------------------------------------------------------------------------   
    
  

    def ShowColorTrajectoryVect(self,index0,Index_Start,xarray,zoom,zlabel,title,figname,dt=EXPOSURE):
        '''
        ShowColorTrajectoryVect(self,index0,Index_Start,xarray,zlabel,title,figname,dt=EXPOSURE):
            input:
                index0: reference index on atm-am conditions
                Index_Start : index of colors with the smallest wavelength
                xarray: array in X ex  VAOD, PWV
                zoom :
            output:
                2D-plot of colors shifted by color-bias, with each point corresponding to 1 SED-1-conditions.
                The colors of plots is linked to the condition.
                example 
                Plot of  G-R vs U-G points, one per SED-condition, the colors of points is related to cond
        
        '''
         # define the max color
        Index_Stop=Index_Start+3
    
        # reserve TWO tables of ColorDiff in 3 indexes
        # - Index 1 : first or second color
        # - Index 2 : SED number
        # - Index 3 : each atmospheric-am condition
        
        ColorDiff=np.zeros([2,self.NBSED,len(xarray)])
        ColorDiff2=np.zeros([2,self.NBSED,len(xarray)])
    
         # the x-array indeed will be colors in this 2D plot
        x= xarray
        for i in range(len(x)):
            c = cmap(int(np.rint(x[i] / x.max() * 255)))

        norm = mpl.colors.Normalize(vmin=x.min(), vmax=x.max())

        fig = plt.figure(figsize=(10, 10))    
    
        # LOOP on SED
        for ised in np.arange(self.NBSED):
            # LOOP on the tow colors to fill the two Tables
            for icol in np.arange(Index_Start,Index_Stop-1):
            
                ColorIndex=icol-Index_Start
                
                
                col=self.colors_bias[ised,:,:]         # get biased colors for all atm event and all colors
                col0=self.colors_bias[ised,index0,:]   # get biased colors for the index0 reference atm event 
                deltacol=col[:,icol]-col0[icol]        # compute the deltacol
                ColorDiff[ColorIndex,ised,:]=deltacol*zoom             # exagerate the deltacol              
                ColorDiff2[ColorIndex,ised,:]=self.colors[ised,:,icol] # original color of the source
         
        #LOOP on SED to plot colors in the 2D plot of colors
        for ised in np.arange(self.NBSED):
            # LOOP on atm-am-conditions 
            for i in np.arange(len(xarray)-1):
                c = cmap(int(np.rint(x[i] / x.max() * 255)))
                X=ColorDiff[0,ised,i]+ColorDiff2[0,ised,0]
                nextX=ColorDiff[0,ised,i+1]+ColorDiff2[0,ised,0]
                Y=ColorDiff[1,ised,i]+ColorDiff2[1,ised,0]
                nextY=ColorDiff[1,ised,i+1]+ColorDiff2[1,ised,0]
                U=nextX-X
                V=nextY-Y
                plt.quiver(X,Y, U, V,scale_units='xy', angles='xy', scale=1, color=c)
                
    
        plt.grid()

        xtitle=number_to_color[Index_Start]+' color (mag)'
        ytitle=number_to_color[Index_Start+1]+' color (mag)'

        plt.xlabel(xtitle,fontsize=20,weight='bold')
        plt.ylabel(ytitle,fontsize=20,weight='bold')        
        plt.title(title,fontsize=25,weight='bold')
        #plt.xlim(-5,5)
        #plt.ylim(-5,5)


        rect = 0.91,0.1,0.02,0.8   #  4-length sequence of [left, bottom, width, height] quantities.     
        ax1 = fig.add_axes(rect,label=zlabel)
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='vertical')
        cb1.set_label(zlabel, rotation=90,fontsize=20,weight='bold')
        plt.savefig(figname)      
        #----------------------------------------------------------------------------------------------------------
    
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
    
    