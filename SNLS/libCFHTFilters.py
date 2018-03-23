#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:20:28 2017

@author: dagoret
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysynphot as S
from scipy.interpolate import interp1d

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 6),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)

files = ['CFHT_MegaPrime_Transmission.dat', \
         'CFHT_Primary_Transmission.dat', \
         'QE_camera_high_res_model.dat', \
         'SNIFS_extinction_buton2012_with_tl_X1_25.dat']

tag_u=["u0","u1","u2","u3","u4","u5","u6","u7","u8","u9"]
tag_g=["g0","g1","g2","g3","g4","g5","g6","g7","g8","g9"]
tag_r=["r0","r1","r2","r3","r4","r5","r6","r7","r8","r9"]
tag_i=["i0","i1","i2","i3","i4","i5","i6","i7","i8","i9"]
tag_z=["z0","z1","z2","z3","z4","z5","z6","z7","z8","z9"]

file_u =["u0.list","u1.list","u2.list","u3.list","u4.list","u5.list","u6.list","u7.list","u8.list","u9.list"]
file_g =["g0.list","g1.list","g2.list","g3.list","g4.list","g5.list","g6.list","g7.list","g8.list","g9.list"]
file_r =["r0.list","r1.list","r2.list","r3.list","r4.list","r5.list","r6.list","r7.list","r8.list","r9.list"]
file_i =["i0.list","i1.list","i2.list","i3.list","i4.list","i5.list","i6.list","i7.list","i8.list","i9.list"]
file_z =["z0.list","z1.list","z2.list","z3.list","z4.list","z5.list","z6.list","z7.list","z8.list","z9.list"]
file_y =["y0.list","y1.list","y2.list","y3.list","y4.list","y5.list","y6.list","y7.list","y8.list","y9.list"]

path_transmissions="all_products_v3_2/MegaCam_v3.2/"




WLMIN=3000. # Minimum wavelength : PySynPhot works with Angstrom
WLMAX=11000. # Minimum wavelength : PySynPhot works with Angstrom

NBINS=7690 # Number of bins between WLMIN and WLMAX
BinWidth=(WLMAX-WLMIN)/float(NBINS) # Bin width in Angstrom
WL=np.linspace(WLMIN,WLMAX,NBINS)   # Array of wavelength in Angstrom


def insert_row(idx, df, df_insert):
    return df.iloc[:idx, ].append(df_insert).append(df.iloc[idx:, ]).reset_index(drop = True)



def insert_boundary(df,tagname,value=0):
    dfnew=insert_row(value,df,pd.DataFrame({'lambda':[WLMIN],tagname:[value]}))
    dfnew=insert_row(len(dfnew),dfnew,pd.DataFrame({'lambda':[WLMAX],tagname:[value]}))
    return dfnew

def Read_all_filter_files(all_files,all_tag):
    all_df=[]
    index=0
    for file in all_files:
        df=pd.read_table(os.path.join(path_transmissions,file),sep=' ' ,skiprows=4,names=["lambda",all_tag[index]],index_col=False)
        df=insert_row(0,df,pd.DataFrame({'lambda':[WLMIN],all_tag[index]:[0.]}))
        df=insert_row(len(df),df,pd.DataFrame({'lambda':[WLMAX],all_tag[index]:[0.]}))
        all_df.append(df) 
        index+=1
        
    return all_df

def Interpolate(all_df,all_tag):
    
    all_filt_arr=[]
    index=0
    for df in all_df:
        func=interp1d(df["lambda"],df[all_tag[index]],kind='linear')
        sf=func(WL)
        all_filt_arr.append(sf)
        index+=1
    return all_filt_arr


#------------------------------------------------------------------------------------
def GetFiltersTransmissions(path="./all_SNLS_transm.csv"):
      
    df=pd.read_csv(path)
    
    df.sort_index(axis=0,ascending=True,inplace=True)
      
    wl_u=df["lambda"]
    u=df["u"]
    
    wl_g=df["lambda"]
    g=df["g"]
    
    wl_r=df["lambda"]
    r=df["r"]
    
    wl_i=df["lambda"]
    i=df["i"]
    
    wl_z=df["lambda"]
    z=df["z"]
  
    return wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z

#---------------------------------------------------------------------------------
def PlotFiltersTransmissions(wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z):
    plt.figure()
    plt.plot(wl_u,u,'b-',lw=2)
    plt.plot(wl_g,g,'g-',lw=2)
    plt.plot(wl_r,r,'r-',lw=2)
    plt.plot(wl_i,i,'y-',lw=2)
    plt.plot(wl_z,z,color='grey',lw=2)
    plt.grid()
    plt.title("Ideal Filters CFHT",weight='bold')
    plt.xlabel("wavelength (nm)",weight='bold')
    plt.ylabel("transmission",weight='bold')
#---------------------------------------------------------------------------------


#--------------------------------------------------------------------------
def GetThroughputAndCCDQE(path="./all_SNLS_transm.csv"):
    
    df=pd.read_csv(path)
     
    df.sort_index(axis=0,ascending=True,inplace=True)
    
    wl=df["lambda"]
    throughput=df["thrpt"]
    mirror=df["mirror"]
    ccdqe=df["qe"]
    
    trans_opt_elec=np.array(throughput*mirror*ccdqe)
    
    #index_to_zero=np.where(np.logical_or(wl<WLMIN,wl>WLMAX))
    #ccdqe[index_to_zero]=0
    
    return wl,throughput,ccdqe,trans_opt_elec
#----------------------------------------------------------------------------------
def PlotThroughputAndCCDQE(wl,throughput,ccdqe,trans_opt_elec):
    plt.figure()
    plt.plot(wl,throughput,'b-',label='throughput',lw=2)
    plt.plot(wl,ccdqe,'r-',label='CCD-QE',lw=2)
    plt.plot(wl,trans_opt_elec,'k-',label='Comb',lw=2)
    plt.grid()
    plt.title("CFHT Throughput and CCD QE",weight='bold')
    plt.xlabel("wavelength",weight='bold')
    plt.ylabel("transmission",weight='bold')
    plt.legend()
#---------------------------------------------------------------------------------

def GetAllCFHTTransmissions(path="./all_SNLS_transm.csv"):
    
    wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z=GetFiltersTransmissions(path)
    wl2,throughput,ccdqe,trans_opt_elec=GetThroughputAndCCDQE(path)
    
    tot_u=np.array(u*trans_opt_elec)
    tot_g=np.array(g*trans_opt_elec)
    tot_r=np.array(r*trans_opt_elec)
    tot_i=np.array(i*trans_opt_elec)
    tot_z=np.array(z*trans_opt_elec)

    return wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z

#---------------------------------------------------------------------------------

def PlotAllCFHTTransmissions(wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z):
    plt.figure()
    plt.plot(wl_u,tot_u,'b-',lw=2)
    plt.plot(wl_g,tot_g,'g-',lw=2)
    plt.plot(wl_r,tot_r,'r-',lw=2)
    plt.plot(wl_i,tot_i,'y-',lw=2)
    plt.plot(wl_z,tot_z,color='black',lw=2)
    plt.grid()
    plt.title("Total transmission CFHT (no atm)",weight='bold')
    plt.xlabel("wavelength (nm)",weight='bold')
    plt.ylabel("transmission",weight='bold')
    plt.savefig("CFHT-total-transm.png")
#----------------------------------------------------------------------------------------------

def GetAllCFHTBands(path):
       
    wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z=GetAllCFHTTransmissions(path)
    
    # Pysynphot requires to convert wl in a true numpy array
    # Pysynphot is not compatile with pandas
    wl_u=np.array(wl_u)
    wl_g=np.array(wl_g)
    wl_r=np.array(wl_r)
    wl_i=np.array(wl_i)
    wl_z=np.array(wl_z)
      
    
    tot_u=np.array(tot_u)
    tot_g=np.array(tot_g)
    tot_r=np.array(tot_r)
    tot_i=np.array(tot_i)
    tot_z=np.array(tot_z)
    
    bp_u = S.ArrayBandpass(wl_u,tot_u, name='CFHT_U')
    bp_g = S.ArrayBandpass(wl_g,tot_g, name='CFHT_G')
    bp_r = S.ArrayBandpass(wl_r,tot_r, name='CFHT_R')
    bp_i = S.ArrayBandpass(wl_i,tot_i, name='CHFT_I')
    bp_z = S.ArrayBandpass(wl_z,tot_z, name='CFHT_Z')


    return bp_u,bp_g,bp_r,bp_i,bp_z

#--------------------------------------------------------------------------


def PlotAllCFHTBands(bp_u,bp_g,bp_r,bp_i,bp_z):
    plt.figure()
    plt.plot(bp_u.wave, bp_u.throughput, 'b',lw=2)
    plt.plot(bp_g.wave, bp_g.throughput, 'g',lw=2)
    plt.plot(bp_r.wave, bp_r.throughput, 'r',lw=2)
    plt.plot(bp_i.wave, bp_i.throughput, 'y',lw=2)
    plt.plot(bp_z.wave, bp_z.throughput, color='grey',lw=2)

    #plt.ylim(0, 1.)
    plt.grid()
    plt.title("Total transmission LSST (no atm) with PySynPhot",weight='bold')
    plt.xlabel(bp_u.waveunits,weight='bold')
    plt.ylabel('throughput',weight='bold')
    plt.legend([bp_u.name, bp_g.name,bp_r.name, bp_i.name,bp_z.name], loc=1,fontsize=9)



#---------------------------------------------------------------------------------
if __name__ == "__main__":
    print 'hello'
    
    wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z=GetFiltersTransmissions("./all_SNLS_transm.csv")
    PlotFiltersTransmissions(wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z)
      
    wl,throughput,ccdqe,trans_opt_elec=GetThroughputAndCCDQE("./all_SNLS_transm.csv")
    PlotThroughputAndCCDQE(wl,throughput,ccdqe,trans_opt_elec)
    
    
    wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z=GetAllCFHTTransmissions("./all_SNLS_transm.csv")
    
    
    
    PlotAllCFHTTransmissions(wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z)
    
    
    bp_u,bp_g,bp_r,bp_i,bp_z=GetAllCFHTBands("./all_SNLS_transm.csv")
    
    
    PlotAllCFHTBands(bp_u,bp_g,bp_r,bp_i,bp_z)
    
    
    
    
    
    
    
