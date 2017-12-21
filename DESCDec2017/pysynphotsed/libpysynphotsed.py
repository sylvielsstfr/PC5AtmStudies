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

# to enlarge the sizes
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 8),
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


TypeStar_to_Temperature = {"O":25000., "B":10000., "A":7500.,"F":6000.,"G":5000.,"K":3500 }
TypeStar_to_color= {"O":'c', "B":'b', "A":'g',"F":'r',"G":'m',"K":"k" }
TypeStar_to_number= {"O":0, "B":1, "A":2,"F":3,"G":4,"K":5 }
Set_Log_g=[0.,1.,2.,3.,4.,5.]
Set_Log_Z=[-2.5,-2.,-1.5,-1.,-.5,0.,0.2,0.5]




#-------------------------------------------------------------------------------------------
def get_ck04models(temp):
    all_sed=[]
    for log_g in Set_Log_g:
        for logz in Set_Log_Z:
            sed = S.Icat('ck04models', temp, logz, log_g) 
            if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                all_sed.append(sed)
    return all_sed
#------------------------------------------------------------------------------------------
def get_all_ck04models():
    
    all_sed=[]   # common SED container

    #for key, temp in TypeStar_to_Temperature.iteritems(): 
    for typestar in TypeStar:
        temp=TypeStar_to_Temperature[typestar] 
        
        all_sub_sed=[]
        for log_g in Set_Log_g:
            for logz in Set_Log_Z:
                sed = S.Icat('ck04models', temp, logz, log_g) 
                if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                    all_sub_sed.append(sed)
        
        all_sed.append(all_sub_sed)
    return all_sed
#---------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
def get_k93models(temp):
    all_sed=[]
    for log_g in Set_Log_g:
        for logz in Set_Log_Z:
            sed = S.Icat('k93models', temp, logz, log_g) 
            if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                all_sed.append(sed)
    return all_sed
#------------------------------------------------------------------------------------------
def get_all_k93models():
    
    all_sed=[]   # common SED container

    #for key, temp in TypeStar_to_Temperature.iteritems(): 
    for typestar in TypeStar:
        temp=TypeStar_to_Temperature[typestar]       

        all_sub_sed=[]
        for log_g in Set_Log_g:
            for logz in Set_Log_Z:             
                sed = S.Icat('k93models', temp, logz, log_g) 
                if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                    all_sub_sed.append(sed)
        
        all_sed.append(all_sub_sed)
    return all_sed
#---------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
def get_phoenixmodels(temp):
    all_sed=[]
    for log_g in Set_Log_g:
        for logz in Set_Log_Z:
            sed = S.Icat('phoenix', temp, logz, log_g) 
            if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                all_sed.append(sed)
    return all_sed
#------------------------------------------------------------------------------------------
def get_all_phoenixmodels():
    
    all_sed=[]   # common SED container

    #for key, temp in TypeStar_to_Temperature.iteritems(): 
    for typestar in TypeStar:
        temp=TypeStar_to_Temperature[typestar]   

        all_sub_sed=[]
        for log_g in Set_Log_g:
            for logz in Set_Log_Z:
                sed = S.Icat('phoenix', temp, logz, log_g) 
                if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                    all_sub_sed.append(sed)
        
        all_sed.append(all_sub_sed)
    return all_sed
#---------------------------------------------------------------------------------------
def plot_allsed_starmodels(all_sed,thetitle,figfilename,yscale='lin',XMIN=3200.,XMAX=10000.):
    
    NBSEDCOLORS=len(all_sed)
    
    for icol in np.arange(NBSEDCOLORS):
        sed_coll=all_sed[icol]
        
        index=0
        for sed in sed_coll:
            
            if index==0:
                if yscale=='log' :
                    plt.semilogy(sed.wave,sed.flux,sed_colors_mpl[icol],label=TypeStar[icol])
                else:
                    plt.plot(sed.wave,sed.flux,sed_colors_mpl[icol],label=TypeStar[icol])
            else:
                if yscale=='log' :
                    plt.semilogy(sed.wave,sed.flux,sed_colors_mpl[icol])
                else:
                    plt.plot(sed.wave,sed.flux,sed_colors_mpl[icol])

            index+=1
    plt.xlim(XMIN, XMAX)
    plt.xlabel(sed.waveunits)
    plt.ylabel(sed.fluxunits)
    plt.grid(True)
    plt.legend()
    plt.title(thetitle)
    plt.savefig(figfilename)
    
    
#-----------------------------------------------------------------------------------------------

def get_all_calspec_hd():
    
    SEDfile_dir=os.path.join(top_pysynphot_data_dir,'calspec')
    filelist=os.listdir(SEDfile_dir) 
    fits_files = [f for f in os.listdir(SEDfile_dir) if f.endswith('.fits')]
    
    # extract header and filename
    star_header = []
    star_file_calspec = []   
    for filename in filelist:       
        index=0
        if re.search('fits',filename) and re.search('hd',filename) and re.search('stis',filename):  #example of filename filter
            index+=1
            fullfilename = os.path.join(SEDfile_dir,filename)
            hdr = fits.getheader(fullfilename)
            star_header.append(hdr)
            star_file_calspec.append(filename)
            
    # extract starname    
    star_names = []
    index=0
    for hdr in star_header: 
#    print index
        if index!=433:
            star_name=star_header[index]['TARGETID']
            star_names.append(star_name)
            index+=1
        else:
            print '>>>>>> skip file # ',index, 'BAD HEADER'
            print '>>>>>> filename = ', filelist[index]
            print hdr
            index+=1
        
    # sort filename    
    star_names_sorted=sorted(star_names,key=star_names.count,reverse=True)
    star_names_sorted_upper = map(lambda s: s.upper(), star_names_sorted)
    star_names_set=set(star_names_sorted_upper)

    # build dictionnary of filename (several files per star name)
    StarDict= {}
    for star in star_names_set:
        print star,': \n'
        star_set_of_file= []
        tag_upper='^'+star.upper()+'*'
        tag_lower='^'+star.lower()+'*'
        
        for thefile in fits_files:
            if re.search(tag_upper,thefile) or re.search(tag_lower,thefile):         
                star_set_of_file.append(thefile)
        #StarDict[star]=sorted(star_set_of_file,key=star_names.count,reverse=True)
        StarDict[star]=sorted(star_set_of_file,key=star_names.count)
        print StarDict[star] ,'\n'
      
    # SED 
    all_sed=[]    
    for keystar in StarDict:
        the_files=StarDict[keystar]
        if(len(the_files))>0 and keystar != 'SUN_REFERENCE':
        #print sorted(the_files,reverse=True)
        
            selected_file=the_files[0]
            selected_fullfile=os.path.join(SEDfile_dir,selected_file)
     
            sed=S.FileSpectrum(selected_fullfile)
            if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                all_sed.append(sed)
    return all_sed
        
#-----------------------------------------------------------------------------------------------
   
def plot_allsed(all_sed,thetitle,figfilename,yscale='lin',XMIN=3200.,XMAX=10000.):
      
    for sed in all_sed:     
        if yscale=='log':
            plt.semilogy(sed.wave,sed.flux)
        else:
            plt.plot(sed.wave,sed.flux)

    plt.xlim(XMIN,XMAX)
    plt.xlabel(sed.waveunits)
    plt.ylabel(sed.fluxunits)
    plt.grid(True)
    plt.title(thetitle)
    plt.savefig(figfilename)
    


#------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    print 'hello'

    #all_sed=get_all_ck04models()
    
    #plot_allsed_starmodels(all_sed,"ck04models","ck04models.png")
    
    all_sed=get_all_calspec_hd()
    plot_allsed(all_sed,'SED of CALSPEC stars','calspec_hd_lin.png',yscale='lin')
    plot_allsed(all_sed,'SED of CALSPEC stars','calspec_hd_log.png',yscale='log')

    
    
    
    
    
    
    
    
    
    
    
    