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


#-------------------------------------------------------------------------------------
class Atmosphere(object):
    def __init__(self,name):
        self.atmcodename = name
        self.nbevent=0
        self.NBWL=0
        self.array = np.array([])
        
        
    def fill_array(self,newarray):
        
        nbwl=len(newarray)
        if self.NBWL==0:
            self.NBWL=nbwl
            self.array=np.array([]).reshape(0,self.NBWL)
            
        if self.NBWL==nbwl:
            self.array = np.r_[self.array, [newarray]]
        else:
            print 'error in the number of wavelength' 
       
    def print_array(self):  
        print self.array
        print type(self.array)
        
    def get_array(self):
        return self.array
#------------------------------------------------------------------------------------
        
    
class LSSTTransmission(object):
    def __init__(self,name,nbbands):
        self.name = name
        self.NBBANDS = nbbands
        self.nbfilled=0
        
        self.passband= []
        
    def fill_passband(self,band):
        self.passband.append(band)
        self.nbfilled+=1
        print 'one band added ',self.nbfilled,' / ',self.NBBANDS
        
    def get_bands(self):
        return self.passband
 
    
        
        
    


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
    
    