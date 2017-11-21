
# coding: utf-8

# # Test to read good filter data
# 
# - author Sylvie Dagoret-Campagne
# - affiliation : LAL/IN2P3/CNRS
# - date May 24th 2017
# - update May 27th 2017
#-------------------------------------------------------------------

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#get_ipython().magic(u'matplotlib inline')


# other paclage

import pandas as pd
import os
import re


# config matplotlib
#-----------------------

# to enlarge the sizes
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)


#file path definition
#----------------------

files_idealfilters=['fdata/ideal_u.txt','fdata/ideal_g.txt','fdata/ideal_r.txt','fdata/ideal_i.txt','fdata/ideal_z.txt','fdata/ideal_y4.txt']

if __name__ == '__main__':
    #Number of files
    #----------------

    NBFILES=len(files_idealfilters)

    #plt.style.use('ggplot') # gnuplot style
    print(plt.style.available)
    print matplotlib.matplotlib_fname() # : where is the matplotlib config file

    # load files
    #-------------

    data_u=np.loadtxt(files_idealfilters[0],skiprows=2)
    data_g=np.loadtxt(files_idealfilters[1],skiprows=2)
    data_r=np.loadtxt(files_idealfilters[2],skiprows=2)
    data_i=np.loadtxt(files_idealfilters[3],skiprows=2)
    data_z=np.loadtxt(files_idealfilters[4],skiprows=2)
    data_y4=np.loadtxt(files_idealfilters[5],skiprows=2)


    # extract data
    #-------------

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


    # Plot figures
    #--------------
    plt.figure()
    plt.plot(wl_u,u,'b-')
    plt.plot(wl_g,g,'g-')
    plt.plot(wl_r,r,'r-')
    plt.plot(wl_i,i,'y-')
    plt.plot(wl_z,z,'k-')
    plt.plot(wl_y4,y4,'-',color='grey')
    plt.grid(True)
    plt.title("Ideal Filters LSST")
    plt.xlabel("wavelength (nm)")
    plt.ylabel("transmission")
    plt.savefig("idealfilterstransm.png")
    plt.show()






