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
NBWLBINS=int(WLMAX-WLMIN) # Number of bins between WLMIN and WLMAX
WLBinWidth=(WLMAX-WLMIN)/float(NBWLBINS) # Bin width in Angstrom
WL=np.linspace(WLMIN,WLMAX,NBWLBINS)   # Array of wavelength in Angstrom


NBROW=len(Temperature_range)*len(Set_Log_Z)*len(Set_Log_G)
NBCOL=NBWLBINS+5


# SED container
# index, validity, Temperature , logg, logz, 
index_num=0  # index for the SED numbers
index_val=1  # flag to indicate if the calculation of the SED is valid
index_temp=2 # temperature of the star in the SED model
index_logg=3 # logg gravity parameter in the SED model
index_logz=4 # logz metallicity parameter in the SED model
index_spec=5 # index where the sed-flux start
index_mag=index_spec+NBWLBINS   # index where the magnitudes would start





SELECTION_MODEL_PHOENIX=False
SELECTION_MODEL_PICKLE_UVI=False
SELECTION_MODEL_PICKLE_UVK=True


 
if SELECTION_MODEL_PHOENIX:
    outputfile_fits='sedgrid_phoenixmodels_all.fits'
    # data container for the SED:
    data=np.zeros((NBROW+1,NBCOL))
    # first row is the wavelengthes
    data[0,index_spec:]=WL
   
elif SELECTION_MODEL_PICKLE_UVI:
    outputfile_fits='sedgrid_pickle_uvi_all.fits'
    SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[9],'dat_uvi')
    fits_files = [f for f in os.listdir(SEDfile_dir) if f.endswith('.fits')]
    if 'pickles.fits' in fits_files:
        fits_files.remove('pickles.fits')    # remove one file
    NBROW=len(fits_files)   # redefine NBROW
    # redefine data
    data=np.zeros((NBROW+1,NBCOL))
    data[0,index_spec:]=WL
    
elif SELECTION_MODEL_PICKLE_UVK:
    outputfile_fits='sedgrid_pickle_uvk_all.fits'
    SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[9],'dat_uvk')
    fits_files = [f for f in os.listdir(SEDfile_dir) if f.endswith('.fits')]
    if 'pickles_uk.fits' in fits_files:
        fits_files.remove('pickles_uk.fits')    # remove one file
    NBROW=len(fits_files)   # redefine NBROW
    # redefine data
    data=np.zeros((NBROW+1,NBCOL))
    data[0,index_spec:]=WL
     
    
else:
    outputfile_fits='sedgrid_unknown_all.fits'
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
def get_grid_pickle_uvi():
    
    fits_files = [f for f in os.listdir(SEDfile_dir) if f.endswith('.fits')]
    
    if 'pickles.fits' in fits_files:
        fits_files.remove('pickles.fits')    # remove one file
    
    
    obj_names = []
    obj_nums = []
    obj_files = []
    
    index=0
    # scan files in the directory
    for thefile in fits_files:
    #thenames=re.findall('^bk_([a-z][0-9]+).fits$',thefile)
        thenames=re.findall('^(.*).fits$',thefile) 
        thenumstr=re.findall('^pickles_(.*).fits$',thefile)
        thenum=int(thenumstr[0])
        thenam=thenames[0]
        if(len(thenames)>0):
            obj_names.append(thenam)
            obj_nums.append(thenum)
            obj_files.append(thefile)
        else:
            print 'bad file ',thefile
        index+=1
    
    
    objames_and_objfiles = zip(obj_names, obj_files)
    objnums_and_objfiles = zip(obj_nums, obj_files)
    
    # make a dictionnary
    OBJDict= {}
    for obj,thefile in objnums_and_objfiles:
        OBJDict[obj]=thefile
      
        
    index=1   
    
    # loop
    for keyobj in OBJDict:
      
        the_file=OBJDict[keyobj]
        
        selected_file=the_file
        selected_fullfile=os.path.join(SEDfile_dir,selected_file)
        
        sed=S.FileSpectrum(selected_fullfile)
        
        
        sed.convert('flam') # to be sure every spectrum is in flam unit
        if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
            data[index,index_val]=1
            
            wave=sed.wave
            flux=sed.flux
            
            #extrapolate
            if(wave[0]>WL[0]):
                wave=np.insert(wave,0,WL[0])
                flux=np.insert(flux,0,0)
            if(wave[0]<WL[-1]):
                wave=np.append(wave,WL[-1])
                flux=np.append(flux,0)
           
            func=interp1d(wave,flux,kind='linear')
            theflux=func(WL)
            data[index,index_spec:]=theflux
        else:
            print 'flux is empty'
            
            
        index+=1
#------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------    
#------------------------------------------------------------------------------------------
def get_grid_pickle_uvk():
    
    fits_files = [f for f in os.listdir(SEDfile_dir) if f.endswith('.fits')]
    
    if 'pickles_uk.fits' in fits_files:
        fits_files.remove('pickles_uk.fits')    # remove one file
        
        
    obj_names = []
    obj_nums = []
    obj_files = []
    
    index=0
    # scan files in the directory
    for thefile in fits_files:
    #thenames=re.findall('^bk_([a-z][0-9]+).fits$',thefile)
        thenames=re.findall('^(.*).fits$',thefile) 
        thenumstr=re.findall('^pickles_uk_(.*).fits$',thefile)
        thenum=int(thenumstr[0])
        thenam=thenames[0]
        if(len(thenames)>0):
            obj_names.append(thenam)
            obj_nums.append(thenum)
            obj_files.append(thefile)
        else:
            print 'bad file ',thefile
        index+=1
    
    
    objames_and_objfiles = zip(obj_names, obj_files)
    objnums_and_objfiles = zip(obj_nums, obj_files)
    
    # make a dictionnary
    OBJDict= {}
    for obj,thefile in objnums_and_objfiles:
        OBJDict[obj]=thefile
      
        
    index=1   
    
    # loop
    for keyobj in OBJDict:
      
        the_file=OBJDict[keyobj]
        
        selected_file=the_file
        selected_fullfile=os.path.join(SEDfile_dir,selected_file)
        
        sed=S.FileSpectrum(selected_fullfile)
        
        
        sed.convert('flam') # to be sure every spectrum is in flam unit
        if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
            data[index,index_val]=1
            
            wave=sed.wave
            flux=sed.flux
            
            #extrapolate
            if(wave[0]>WL[0]):
                wave=np.insert(wave,0,WL[0])
                flux=np.insert(flux,0,0)
            if(wave[0]<WL[-1]):
                wave=np.append(wave,WL[-1])
                flux=np.append(flux,0)
         
            
            func=interp1d(wave,flux,kind='linear')
            theflux=func(WL)
            data[index,index_spec:]=theflux
            data[index,index_temp]=keyobj
        else:
            print 'flux is empty'
            
            
        index+=1
#------------------------------------------------------------------------------------------


                
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
def plot_sedimg():
    plt.figure(figsize=(15,10))   
    
    fluxmin=data[1:,index_spec:].min()
    fluxmax=data[1:,index_spec:].max()
    
    
    
    img=plt.imshow(data[1:,index_spec:],origin='lower',vmin=fluxmin,vmax=fluxmax,cmap='jet')
    plt.colorbar(img)
    plt.grid(True)
    plt.title('sed grid')
    plt.xlabel('bin number of wavelengths')
    plt.ylabel('sed number')
    plt.show()
    
#------------------------------------------------------------------------------------------------

  #---------------------------------------------------------------------------------
def plot_allsed(figsize=(15,8)):
    plt.figure()   
       
    fluxmin=data[1:,index_spec:].min()
    fluxmax=data[1:,index_spec:].max()
    
    
    for idx in np.arange(NBROW-1):
        plt.semilogy(WL,data[1+idx,index_spec:],'-')
        
        
    plt.grid(True)
    plt.title('sed')
    plt.xlabel('wavelength (A)')
    plt.ylabel('sed (flam)')
    plt.ylim(fluxmin,fluxmax)
    plt.show()
    
#------------------------------------------------------------------------------------------------
      

#------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    
    print '*************************************************************************************************************************'
    print 'Parameters for the SIMULATION OF SED'
    print '\t - wavelength : WLMIN= ',WLMIN,' WLMAX=',WLMAX,' WLBinWidth=',WLBinWidth,' NBWLBINS=',NBWLBINS,' Len-WL=',len(WL)
    print '\t - NBROW=',NBROW
    print '\t - NBCOL',NBCOL
    print '\t - index_num=',index_num
    print '\t - index_val=',index_val
    print '\t - index_temp=',index_temp
    print '\t - index_logg=',index_logg
    print '\t - index_logz=',index_logz
    print '\t - index_spec=',index_spec
    print '\t - index_mag=',index_mag
    print '\t - data shape=',data.shape
    print '\t - output-file =',outputfile_fits
    print '*************************************************************************************************************************'
    


    
  
    # simulate SED
    if SELECTION_MODEL_PHOENIX:
        get_grid_phoenixmodels()
        #plot_allsed()
        #plot_sedimg()
    elif SELECTION_MODEL_PICKLE_UVI:
        get_grid_pickle_uvi()
        plot_sedimg()
        plot_allsed()
    elif SELECTION_MODEL_PICKLE_UVK:
        get_grid_pickle_uvk()
        plot_sedimg()
        plot_allsed()
        
        
    # save SED in a fits file    
    hdr = fits.Header()
    hdr['NBSED'] = NBROW
    hdr['NBWLBIN']=NBWLBINS
    hdr['WLMIN']=WLMIN
    hdr['WLMAX']=WLMAX
    hdr['WLBINWDT']=WLBinWidth
    
    
    
    
    if SELECTION_MODEL_PHOENIX:
        hdr['SEDMODEL'] = 'phoenix'
        hdr['TMIN'] = TMIN
        hdr['TMAX'] = TMAX
        hdr['TSTEP'] = TSTEP
        hdr['LOGZ'] =np.array_str(Set_Log_Z)
        hdr['LOGG']=np.array_str(Set_Log_G)
        hdr['IDX_NUM']=index_num
        hdr['IDX_VAL']=index_val
        hdr['IDX_TEMP']=index_temp
        hdr['IDX_LOGG']=index_logg
        hdr['IDX_LOGZ']=index_logz
        hdr['IDX_SPEC']=index_spec
        hdr['IDX_MAG']=index_mag
    elif SELECTION_MODEL_PICKLE_UVI:
        hdr['SEDMODEL'] = 'pickle_uvi'
        hdr['IDX_SPEC']=index_spec
        hdr['IDX_MAG']=index_mag
        hdr['IDX_PKNU']=index_temp
    elif SELECTION_MODEL_PICKLE_UVK:
        hdr['SEDMODEL'] = 'pickle_uvk'
        hdr['IDX_SPEC']=index_spec
        hdr['IDX_MAG']=index_mag
        hdr['IDX_PKNU']=index_temp
    else:
        hdr['SEDMODEL'] = 'unknown'
        hdr['IDX_SPEC']=index_spec
        hdr['IDX_MAG']=index_mag

    
    print hdr
    
    hdu = fits.PrimaryHDU(data,header=hdr)
    hdu.writeto(outputfile_fits,overwrite=True)