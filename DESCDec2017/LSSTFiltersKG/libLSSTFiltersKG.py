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


files_idealfilters=['LSSTFiltersKG/fdata/ideal_u.txt','LSSTFiltersKG/fdata/ideal_g.txt','LSSTFiltersKG/fdata/ideal_r.txt','LSSTFiltersKG/fdata/ideal_i.txt','LSSTFiltersKG/fdata/ideal_z.txt','LSSTFiltersKG/fdata/ideal_y4.txt']
file_lsstoptccd='LSSTFiltersKG/fdata/LSST-ThroughputCCD.xlsx'

#------------------------------------------------------------------------------------
def GetFiltersTransmissions(path):
    data_u=np.loadtxt(os.path.join(path,files_idealfilters[0]),skiprows=2)
    data_g=np.loadtxt(os.path.join(path,files_idealfilters[1]),skiprows=2)
    data_r=np.loadtxt(os.path.join(path,files_idealfilters[2]),skiprows=2)
    data_i=np.loadtxt(os.path.join(path,files_idealfilters[3]),skiprows=2)
    data_z=np.loadtxt(os.path.join(path,files_idealfilters[4]),skiprows=2)
    data_y4=np.loadtxt(os.path.join(path,files_idealfilters[5]),skiprows=2)
    
    wl_u=data_u[:,0]
    u=data_u[:,1]
#
    wl_g=data_g[:,0]
    g=data_g[:,1]
#
    wl_r=data_r[:,0]
    r=data_r[:,1]
#
    wl_i=data_i[:,0]
    i=data_i[:,1]
#
    wl_z=data_z[:,0]
    z=data_z[:,1]
#
    wl_y4=data_y4[:,0]
    y4=data_y4[:,1]
    
    return wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z,wl_y4,y4
#---------------------------------------------------------------------------------
def PlotFiltersTransmissions(wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z,wl_y4,y4):
    plt.figure()
    plt.plot(wl_u,u,'b-')
    plt.plot(wl_g,g,'g-')
    plt.plot(wl_r,r,'r-')
    plt.plot(wl_i,i,'y-')
    plt.plot(wl_z,z,'k-')
    plt.plot(wl_y4,y4,'-',color='grey')
    plt.grid()
    plt.title("Ideal Filters LSST")
    plt.xlabel("wavelength (nm)")
    plt.ylabel("transmission")
#---------------------------------------------------------------------------------


#--------------------------------------------------------------------------
def GetThroughputAndCCDQE(path):
    data_throuthput=pd.read_excel(os.path.join(path,file_lsstoptccd),skiprow=1)
    
    wl=np.array(data_throuthput["WL"])
    throughput=np.array(data_throuthput["THROUGHPUT"])
    ccdqe=np.array(data_throuthput["CCD2"])
    trans_opt_elec=np.array(data_throuthput["THROUGHPUT"]*data_throuthput["CCD2"])
    return wl,throughput,ccdqe,trans_opt_elec
#----------------------------------------------------------------------------------
def PlotThroughputAndCCDQE(wl,throughput,ccdqe,trans_opt_elec):
    plt.figure()
    plt.plot(wl,throughput,'b-',label='throughput')
    plt.plot(wl,ccdqe,'r-',label='CCD-QE')
    plt.plot(wl,trans_opt_elec,'k-',label='Comb')
    plt.grid()
    plt.title("LSST Throughput and CCD QE")
    plt.xlabel("wavelength")
    plt.ylabel("transmission")
    plt.legend()
#---------------------------------------------------------------------------------

def GetAllLSSTTransmissions(path):
    wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z,wl_y4,y4=GetFiltersTransmissions(path)
    wl2,throughput,ccdqe,trans_opt_elec=GetThroughputAndCCDQE(path)
    
    tot_u=np.array(u*trans_opt_elec)
    tot_g=np.array(g*trans_opt_elec)
    tot_r=np.array(r*trans_opt_elec)
    tot_i=np.array(i*trans_opt_elec)
    tot_z=np.array(z*trans_opt_elec)
    tot_y4=np.array(y4*trans_opt_elec)   
    return wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z,wl_y4,tot_y4

#---------------------------------------------------------------------------------

def PlotAllLSSTTransmissions(wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z,wl_y4,tot_y4):
    plt.figure()
    plt.plot(wl_u,tot_u,'b-')
    plt.plot(wl_g,tot_g,'g-')
    plt.plot(wl_r,tot_r,'r-')
    plt.plot(wl_i,tot_i,'y-')
    plt.plot(wl_z,tot_z,'k-')
    plt.plot(wl_y4,tot_y4,'-',color='grey')
    plt.grid()
    plt.title("Total transmission LSST (no atm)")
    plt.xlabel("wavelength (nm)")
    plt.ylabel("transmission")
    plt.savefig("lsst-total-transm.png")
#----------------------------------------------------------------------------------------------

def GetAllLSSTBands(path):
    
    wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z,wl_y4,tot_y4=GetAllLSSTTransmissions(path)
    
    bp_u = S.ArrayBandpass(wl_u*10.,tot_u, name='LSST_U')
    bp_g = S.ArrayBandpass(wl_g*10.,tot_g, name='LSST_G')
    bp_r = S.ArrayBandpass(wl_r*10,tot_r, name='LSST_R')
    bp_i = S.ArrayBandpass(wl_i*10.,tot_i, name='LSST_I')
    bp_z = S.ArrayBandpass(wl_z*10.,tot_z, name='LSST_Z')
    bp_y4 = S.ArrayBandpass(wl_y4*10,tot_y4, name='LSST_Y4')
    return bp_u,bp_g,bp_r,bp_i,bp_z,bp_y4

#----------------------------------------------------------------------------------------------
def PlotAllLSSTBands(bp_u,bp_g,bp_r,bp_i,bp_z,bp_y4):
    plt.figure()
    plt.plot(bp_u.wave, bp_u.throughput, 'b')
    plt.plot(bp_g.wave, bp_g.throughput, 'g')
    plt.plot(bp_r.wave, bp_r.throughput, 'r')
    plt.plot(bp_i.wave, bp_i.throughput, 'y')
    plt.plot(bp_z.wave, bp_z.throughput, color='grey')
    plt.plot(bp_y4.wave, bp_y4.throughput,'k')
    plt.ylim(0, 1.)
    plt.grid()
    plt.title("Total transmission LSST (no atm) with PySynPhot")
    plt.xlabel(bp_u.waveunits)
    plt.ylabel('throughput')
    plt.legend([bp_u.name, bp_g.name,bp_r.name, bp_i.name,bp_z.name, bp_y4.name], loc=1,fontsize=9)

#---------------------------------------------------------------------------------
if __name__ == "__main__":
    print 'hello'
    
    wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z,wl_y4,y4=GetFiltersTransmissions('..')
    PlotFiltersTransmissions(wl_u,u,wl_g,g,wl_r,r,wl_i,i,wl_z,z,wl_y4,y4)
      
    wl,throughput,ccdqe,trans_opt_elec=GetThroughputAndCCDQE('..')
    PlotThroughputAndCCDQE(wl,throughput,ccdqe,trans_opt_elec)
    
    
    wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z,wl_y4,tot_y4=GetAllLSSTTransmissions('..')
    PlotAllLSSTTransmissions(wl_u,tot_u,wl_g,tot_g,wl_r,tot_r,wl_i,tot_i,wl_z,tot_z,wl_y4,tot_y4)
    
    
    bp_u,bp_g,bp_r,bp_i,bp_z,bp_y4=GetAllLSSTBands('..')
    PlotAllLSSTBands(bp_u,bp_g,bp_r,bp_i,bp_z,bp_y4)
    
    
    
    
    
    
    
