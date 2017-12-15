#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 14:37:29 2017

@author: dagoret
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import os

# to enlarge the sizes
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)

top_pysynphot_data_dir=os.environ['PYSYN_CDBS']
import pysynphot as S
S.primary_area=6*1e4
S.binning=10.


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
def plot_allsed(all_sed,thetitle,figfilename):
    
    NBSEDCOLORS=len(all_sed)
    
    for icol in np.arange(NBSEDCOLORS):
        sed_coll=all_sed[icol]
        
        index=0
        for sed in sed_coll:
            
            if index==0:
                plt.semilogy(sed.wave,sed.flux,sed_colors_mpl[icol],label=TypeStar[icol])
            else:
                plt.semilogy(sed.wave,sed.flux,sed_colors_mpl[icol])

            index+=1
    plt.xlim(0, 11000)
    plt.xlabel(sed.waveunits)
    plt.ylabel(sed.fluxunits)
    plt.grid(True)
    plt.legend()
    plt.title(thetitle)
    plt.savefig(figfilename)
    
    
#-----------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    print 'hello'

    all_sed=get_all_ck04models()
    
    plot_allsed(all_sed,"ck04models","ck04models.png")

    
    
    
    
    
    
    
    
    
    
    
    