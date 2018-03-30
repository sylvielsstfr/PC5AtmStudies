#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 14:37:29 2017

@author: dagoret
"""

import numpy as np
import matplotlib.pyplot as plt

import os
import re

from astropy.io import fits
from scipy.interpolate import interp1d


# to enlarge the sizes
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 6),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)

import pysynphot as S
S.primary_area=6*1e4
S.binning=10.



# directories
top_pysynphot_data_dir=os.environ['PYSYN_CDBS']
dir_star='calspec'
dir_nostar='grid'
dir_submodels=['ags','bpgs','extinction','jacobi','phoenix','bc95','bz77','galactic','k93models','pickles','bkmodels','ck04models','gunnstryker','kc96']


#star models
sed_colors_mpl=['c','b', 'g', 'r', 'm', 'y', 'k', 'w']
TypeStar=["O","B","A","F","G","K"]

TMAX=50000.
TMIN=4000.
TSTEP=100.

Temperature_range=np.arange(TMIN,TMAX,TSTEP)
TypeStar_to_Temperature = {"O":25000., "B":10000., "A":7500.,"F":6000.,"G":5000.,"K":3500 }
TypeStar_to_color= {"O":'c', "B":'b', "A":'g',"F":'r',"G":'m',"K":"k" }
TypeStar_to_number= {"O":0, "B":1, "A":2,"F":3,"G":4,"K":5 }
Set_Log_G=np.array([0.,1.,2.,3.,4.,5.])
Set_Log_Z=np.array([-2.5,-2.,-1.5,-1.,-.5,0.,0.2,0.5])



# wavelength definitions
WLMIN=3000. # Minimum wavelength : PySynPhot works with Angstrom
WLMAX=11000. # Minimum wavelength : PySynPhot works with Angstrom
NBWLBINS=7000 # Number of bins between WLMIN and WLMAX
BinWidth=(WLMAX-WLMIN)/float(NBWLBINS) # Bin width in Angstrom
WL=np.linspace(WLMIN,WLMAX,NBWLBINS)   # Array of wavelength in Angstrom


NBROW=len(Temperature_range)*len(Set_Log_Z)*len(Set_Log_G)
NBCOL=NBWLBINS+5


# SED container
# index, validity, Temperature , logg, logz, 
index_num=0
index_val=1
index_temp=2
index_logg=3
index_logz=4
index_spec=5
index_mag=index_spec+NBWLBINS

data=np.zeros((NBROW+1,NBCOL))
data[0,index_spec:]=WL


#------------------------------------------------------------------------------------------
def get_grid_phoenixmodels():
    
    index=1
    for temp in Temperature_range:
        for logg in Set_Log_G:
            for logz in Set_Log_Z:
                data[index,index_temp]=temp
                data[index,index_logg]=logg
                data[index,index_logz]=logz
                
                
                sed = S.Icat('phoenix', temp, logz, logg) 
                sed.convert('flam') # to be sure every spectrum is in flam unit
                if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                    data[index,index_val]=1
                    func=interp1d(sed.wave,sed.flux,kind='cubic')
                    flux=func(WL)
                    data[index,index_spec:]=flux
                index+=1
#------------------------------------------------------------------------------------------
def get_grid_phoenixmodels_extinct(extvalue):
    
    index=1
    for temp in Temperature_range:
        for logg in Set_Log_G:
            for logz in Set_Log_Z:
                data[index,index_temp]=temp
                data[index,index_logg]=logg
                data[index,index_logz]=logz
                
                
                sed = S.Icat('phoenix', temp, logz, logg) * S.Extinction(extvalue,'mwavg')
                sed.convert('flam') # to be sure every spectrum is in flam unit
                if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                    data[index,index_val]=1
                    func=interp1d(sed.wave,sed.flux,kind='cubic')
                    flux=func(WL)
                    data[index,index_spec:]=flux
                index+=1
  
    
#------------------------------------------------------------------------------------------

def FitsToPySynphotSED(file_fits):
    
    all_sed=[]
    all_indexes_inFits=[]
    
    hdul = fits.open(file_fits)
    data = hdul[0].data
    
    
    wl=data[0,index_spec:]
    
    good_indexes=np.where(data[0:,index_val]>0)[0]
    
    
    # loop on good spectra only
    for index in good_indexes:
        flux=data[index,index_spec:]
        
        sp = S.ArraySpectrum(wave=wl, flux=flux, waveunits='angstrom', fluxunits='flam')
        all_sed.append(sp)
        all_indexes_inFits.append(index)
        
    all_indexes_inFits=np.array(all_indexes_inFits)
    return all_sed,all_indexes_inFits   

#---------------------------------------------------------------------------------
def plot_allsed():
    plt.figure()   
    img=plt.imshow(data[1:,index_spec:],origin='lower',cmap='jet')
    plt.colorbar(img)
    plt.grid(True)
    plt.show()
    
#------------------------------------------------------------------------------------------------

  #---------------------------------------------------------------------------------
def plot_allsed2():
    plt.figure()   
    img=plt.imshow(data[1:,index_spec:],origin='lower',cmap='jet')
    plt.colorbar(img)
    plt.grid(True)
    plt.show()
    
#------------------------------------------------------------------------------------------------
      

#------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    print 'hello'

    #all_sed=get_all_ck04models()
    
    #plot_allsed_starmodels(all_sed,"ck04models","ck04models.png")
    
    Flag_CALSPEC_HD=False
    Flag_BC95_Z3=False
    Flag_BC95=False
    Flag_KC93_Z0=False
    Flag_THERMALBB=False
    Flag_PHOENIX=True
    Flag_CK04=False
    Flag_PICKLE=False
    Flag_K93=False
    Flag_BK=False
    Flag_BPGS=False
    
  
    
    if Flag_PHOENIX:
        #get_grid_phoenixmodels()
        get_grid_phoenixmodels_extinct(5.0)
        plot_allsed()
        
    hdr = fits.Header()
    
    hdr['NBSED'] = NBROW
    hdr['NBWLBIN']=NBWLBINS
    hdr['WLMIN']=WLMIN
    hdr['WLMAX']=WLMAX
    hdr['SEDMODEL'] = 'phoenix'
    hdr['TMIN'] = TMIN
    hdr['TMAX'] = TMAX
    hdr['TSTEP'] = TSTEP
    hdr['LOGZ'] =np.array_str(Set_Log_Z)
    hdr['LOGG']=np.array_str(Set_Log_G)
    
    print hdr
    
    hdu = fits.PrimaryHDU(data,header=hdr)
    hdu.writeto('sedgrid_phoenixmodels_extinct_50.fits',overwrite=True)