#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 09:52:17 2018

@author: dagoret
"""

import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as cm
import matplotlib as mpl
cmap = cm.jet

import pysynphot as S
from scipy.interpolate import interp1d
import astropy.units as u
from astropy import constants as const





#PATH_LSSTFiltersKG='../LSSTFiltersKG'
#PATH_ATMPARAMSIM='../atmparamsim'
#PATH_ATMTRANSPSIM='../libradtransim'
#PATH_CADENCE='../cadence'
#PATH_LSSTPHOTO='../lsstphotometry'
#PATH_SED='../pysynphotsed'

#sys.path.append(PATH_LSSTFiltersKG)
#sys.path.append(PATH_ATMPARAMSIM)
#sys.path.append(PATH_ATMTRANSPSIM)
#sys.path.append(PATH_CADENCE)
#sys.path.append(PATH_LSSTPHOTO)
#sys.path.append(PATH_SED)


#import libLSSTFiltersKG as lsst
#import libsimulateTranspLSSTScattAbsAer as atm
#import libLSSTPhotometry as photo

top_pysynphot_data_dir=os.environ['PYSYN_CDBS']
sys.path.append('./snlsphotometry')
sys.path.append('../DESCDec2017/pysynphotsed')

import libSNLSPhotometry as photo
import libCFHTFilters as cfht



cfht_transmissionfile="./all_SNLS_transm.csv"


NBBANDS=5
band_to_number={'u':0,'g':1,'r':2,'i':3,'z':4}
number_to_band={0:'u',1:'g',2:'r',3:'i',4:'z'}
filtercolor=['blue','green','red','orange','black']
NBCOLORS=NBBANDS-1
number_to_color={0:'U-G',1:'G-R',2:'R-I',3:'I-Z'}
color_to_number={'U-G':0,'G-R':1,'R-I':2,'I-Z':3}
mpl_colors_col=['b','g','r','k']

WLMIN=3000. # Minimum wavelength : PySynPhot works with Angstrom
WLMAX=11000. # Minimum wavelength : PySynPhot works with Angstrom

NBINS=int(WLMAX-WLMIN) # Number of bins between WLMIN and WLMAX
BinWidth=(WLMAX-WLMIN)/float(NBINS) # Bin width in Angstrom
WL=np.linspace(WLMIN,WLMAX,NBINS)   # Array of wavelength in Angstrom



CFHT_COLL_SURF=np.pi/4.*(3.6*u.m)**2/(u.cm)**2  # LSST collectif surface
S.refs.setref(area=CFHT_COLL_SURF.decompose(), waveset=None)
S.refs.set_default_waveset(minwave=WLMIN, maxwave=WLMAX, num=NBINS, delta=BinWidth, log=False)
S.refs.showref()

EXPOSURE=30.0                      # LSST Exposure time



#------------------------------------------------------------------------------------------
#  Compute the multiplicative factor as calcilated for SpectractorSim to be used for AuxTel
#-------------------------------------------------------------------------------------------
Tel_Surf=CFHT_COLL_SURF*(u.cm)**2            # collection surface of telescope
Time_unit=1*u.s                              # flux for 1 second
SED_unit=1*u.erg/u.s/(u.cm)**2/(u.nanometer) # Units of SEDs in flam (erg/s/cm2/nm)
hc=const.h*const.c                           # h.c product of fontamental constants c and h 
wl_dwl_unit=(u.nanometer)**2                 # lambda.dlambda  in wavelength in nm
g_elec=3.0                                   # electronic gain : elec/ADU
g_disperser_ronchi=0.2                       # theoretical gain for order+1 : 20%
#Factor=2.1350444e11
Factor=(Tel_Surf*SED_unit*Time_unit*wl_dwl_unit/hc/g_elec*g_disperser_ronchi).decompose()
#-------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
    
    # power law
    #-----------
    #plt.figure()
    #pl = S.PowerLaw(10000, -2,waveunits='Angstrom',fluxunits='flam')
    #plt.loglog(pl.wave, pl.flux)
    #plt.axvline(10000, ls='--', color='k')
    #plt.axhline(1, ls='--', color='k')
    #plt.xlabel(pl.waveunits)
    #plt.ylabel(pl.fluxunits)
    #plt.title(pl.name)
    #plt.show()
    
    # get a typical spectrum in flam unit (erg/sec/cm)
    plt.figure()
    flatsp = S.FlatSpectrum(1, fluxunits='flam')
    plt.plot(flatsp.wave, flatsp.flux)
    plt.xlabel(flatsp.waveunits)
    plt.ylabel(flatsp.fluxunits)
    plt.title(flatsp.name)
    plt.show()
    
    
    #  CFHT passband in Angstrom
    bp_u,bp_g,bp_r,bp_i,bp_z=cfht.GetAllCFHTBands(cfht_transmissionfile)
    cfht.PlotAllCFHTBands(bp_u,bp_g,bp_r,bp_i,bp_z)
    
    

    
    # atmosphere simulation 
    #am=1.
    #pwv=5.
    #ozone=300.
    #lambda0_aerosol=500.
    #tau_aerosol=0.05
    
    ##### call to libradtran
    #path,thefile=atm.ProcessSimulationaer(am,pwv,ozone,lambda0_aerosol,tau_aerosol) 
    #fullfilename=os.path.join(path,thefile)
    #atm_data=np.loadtxt(fullfilename)
    #wl_atm=atm_data[:,0]
    #tr_atm=atm_data[:,1] 
    
    #index=0
    #photo_atm=photo.Atmosphere('libradtran')
    #if index==0:
    #    photo_atm.fill_array(wl_atm)        
    #photo_atm.fill_array(tr_atm)        
    #index+=1
    #
    #theatmosph=photo_atm.get_array()
    #photo_atm.plot_pys_bp()
    # retrieve the atmosphere passband
    #all_bp_atm=photo_atm.get_pys_pb()
    
    photo_atm=photo.Atmosphere('SNLS atmosphere')
    df=pd.read_csv(cfht_transmissionfile)
    df.sort_index(axis=0,ascending=True,inplace=True)     
    wl_atm=df["lambda"]
    tr_atm=df["atm"]
    wl_atm=np.array(wl_atm)   
    tr_atm=np.array(tr_atm)
    photo_atm.fill_array(wl_atm)        
    photo_atm.fill_array(tr_atm)  
    
    theatmosph=photo_atm.get_array()
    photo_atm.plot_pys_bp()
    all_bp_atm=photo_atm.get_pys_pb()
    
    
    
    # CFHT detector and fill its passbands
    cfhtdetector=photo.SNLSTransmission('cfhtel')
    
    
    cfhtdetector.fill_det_allbands([bp_u,bp_g,bp_r,bp_i,bp_z])
    cfhtdetector.fill_atm_allevents(all_bp_atm)
    
    all_transmissions=cfhtdetector.make_transmissions()
    cfhtdetector.plot_transmissions()
    
    # Compute the observations    

    #  OBSERVATION
    #
    #  So apprently, pysynphot do not more than a multiplication of SED by the transmission
    #  The counts in the band must be done by us
    #
    the_observation=photo.SNLSObservation('SNLSObs')               # create a set of observation 
    the_observation.fill_sed([flatsp])                             # get the SED from the SED model model
    the_observation.fill_transmission(all_transmissions)           # provide LSST Trroughput transmission
    the_observation.make_observations()                            # start calculations 
    the_observation.make_samplobservations()
    the_observation.compute_counts()
    the_observation.compute_magnitude()
    the_observation.compute_colors()
    
    the_observation.plot_observations(0)
    the_observation.plot_samplobservations(0)
    
    
    
    # for debugging
    if 0:
    
        # calculate by hand
        final_pb_u=bp_u*all_bp_atm[0]
        final_pb_g=bp_g*all_bp_atm[0]
        final_pb_r=bp_r*all_bp_atm[0]
        final_pb_i=bp_i*all_bp_atm[0]
        final_pb_z=bp_z*all_bp_atm[0]
   
    
        flatsp.convert('photlam')
        final_flux_u=S.Observation(flatsp,final_pb_u)
        final_flux_g=S.Observation(flatsp,final_pb_g)
        final_flux_r=S.Observation(flatsp,final_pb_r)
        final_flux_i=S.Observation(flatsp,final_pb_i)
        final_flux_z=S.Observation(flatsp,final_pb_z)
 
    
    
        plt.figure()
        plt.plot(final_pb_u.wave,final_pb_u.throughput,'b:')
        plt.xlabel(final_pb_u.waveunits.name)
        #plt.ylabel(final_pb_u.throughputunits.name)
        plt.plot(final_pb_g.wave,final_pb_g.throughput,'g:')
        plt.plot(final_pb_r.wave,final_pb_r.throughput,'r:')
        plt.plot(final_pb_i.wave,final_pb_i.throughput,'y:')
        plt.plot(final_pb_z.wave,final_pb_z.throughput,'k:')
    
        plt.title('Pure passband')
        plt.grid()
        plt.show()
    
        plt.figure()
        plt.plot(final_pb_u.wave,final_pb_u.throughput*final_pb_u.wave/(S.units.C*S.units.H),'b:')
        plt.xlabel(final_pb_u.waveunits.name)
        #plt.ylabel(final_pb_u.throughputunits.name)
        plt.plot(final_pb_g.wave,final_pb_g.throughput*final_pb_g.wave/(S.units.C*S.units.H),'g:')
        plt.plot(final_pb_r.wave,final_pb_r.throughput*final_pb_r.wave/(S.units.C*S.units.H),'r:')
        plt.plot(final_pb_i.wave,final_pb_i.throughput*final_pb_i.wave/(S.units.C*S.units.H),'y:')
        plt.plot(final_pb_z.wave,final_pb_z.throughput*final_pb_z.wave/(S.units.C*S.units.H),'k:')
        plt.title('passband * WL')
        plt.grid()
        plt.show()
    
    
        plt.figure()
        plt.plot(final_flux_u.wave,final_flux_u.flux,'b-.')
        #plt.xlabel(final_flux_u.waveunits.name)
        #plt.ylabel(final_flux_u.fluxunits.name)
        plt.plot(final_flux_g.wave,final_flux_g.flux,'g-.')
        plt.plot(final_flux_r.wave,final_flux_r.flux,'r-.')
        plt.plot(final_flux_i.wave,final_flux_i.flux,'y-.')
        plt.plot(final_flux_z.wave,final_flux_z.flux,'k-.')
        plt.grid()
        plt.show()
    
    all_magnitudes=the_observation.get_magnitudes()
    
    
    print all_magnitudes[0]
    