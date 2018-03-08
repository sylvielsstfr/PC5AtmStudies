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
            sed.convert('flam') # to be sure every spectrum is in flam unit
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
                sed.convert('flam') # to be sure every spectrum is in flam unit
                if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                    all_sub_sed.append(sed)
        
        all_sed.append(all_sub_sed)
    return all_sed
#---------------------------------------------------------------------------------------  
def get_many_ck04models():
    
    all_Temp=np.linspace(3500.,50000.,50)
    all_sed=[]
    
    for temp in all_Temp:
        set_of_sed=get_ck04models(temp)
        all_sed=np.concatenate((all_sed,set_of_sed))
        
    return all_sed    

#-------------------------------------------------------------------------------------------
def get_k93models(temp):
    all_sed=[]
    for log_g in Set_Log_g:
        for logz in Set_Log_Z:
            sed = S.Icat('k93models', temp, logz, log_g) 
            sed.convert('flam') # to be sure every spectrum is in flam unit
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
                sed.convert('flam') # to be sure every spectrum is in flam unit
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
            sed.convert('flam') # to be sure every spectrum is in flam unit
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
                sed.convert('flam') # to be sure every spectrum is in flam unit
                if(max(sed.flux)>0): # remove empty fluxes because of bad parameters
                    all_sub_sed.append(sed)
        
        all_sed.append(all_sub_sed)
    return all_sed
#------------------------------------------------------------------------------------------
def get_all_pickle():
    
    all_sed=[]   # common SED container
    
    SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[9],'dat_uvi')
    fits_files = [f for f in os.listdir(SEDfile_dir) if f.endswith('.fits')]
    fits_files.remove('pickles.fits')
    
    obj_headers = []
    obj_files = []
    for filename in fits_files:
        index=0
        if re.search('fits',filename):  #example of filename filter
            index+=1
            fullfilename = os.path.join(SEDfile_dir,filename)
            hdr = fits.getheader(fullfilename)
            obj_headers.append(hdr)
            obj_files.append(filename)
            
    obj_names2 = []
    index=0
    for thefile in fits_files:
        #thenames=re.findall('^bk_([a-z][0-9]+).fits$',thefile)
        thenames=re.findall('^(.*).fits$',thefile) 
        if(len(thenames)>0):
            obj_names2.append(thenames[0])
        else:
            print 'bad file ',thefile
        index+=1
        
    obj_names=obj_names2
    objames_and_objfiles = zip(obj_names, obj_files)
    OBJDict= {}
    for obj,thefile in objames_and_objfiles:
        #print obj,': '
        OBJDict[obj]=thefile
        #print OBJDict[obj] 
        
    for keyobj in OBJDict:
        the_file=OBJDict[keyobj]
        
        selected_file=the_file
        selected_fullfile=os.path.join(SEDfile_dir,selected_file)
        
        sed=S.FileSpectrum(selected_fullfile)
        sed.convert('flam') # to be sure every spectrum is in flam unit
        all_sed.append(sed)
        
    return all_sed
      
    
#---------------------------------------------------------------------------------------
def plot_allsed_starmodels(all_sed,thetitle,figfilename,yscale='lin',XMIN=3200.,XMAX=10000.,YMIN=0,YMAX=0):
    
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
    
    if(YMAX>0):
        plt.ylim(YMIN,YMAX)
        
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
 
def get_all_thermalbb(N=10,T0=6000.,sigT=100.):
    
    all_sed=[]   
    
    all_random_values=np.random.standard_normal(N)
    for index in np.arange(N):
        T=T0+all_random_values[index]*sigT
        
        print T
        
        sed=S.BlackBody(T)
        sed.convert('flam') # to be sure every spectrum is in flam unit
        all_sed.append(sed)
        
    return all_sed

#------------------------------------------------------------------------------------------
def get_all_phoenix():
    
     all_sed=[]
     
     # define the top directory
     SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[4])
     
     filelist1=os.listdir(SEDfile_dir+'/phoenixm00') 
     filelist2=os.listdir(SEDfile_dir+'/phoenixm05') 
     filelist3=os.listdir(SEDfile_dir+'/phoenixm10') 
     filelist4=os.listdir(SEDfile_dir+'/phoenixm15') 
     filelist5=os.listdir(SEDfile_dir+'/phoenixm20') 
     filelist6=os.listdir(SEDfile_dir+'/phoenixm25') 
     filelist7=os.listdir(SEDfile_dir+'/phoenixm30') 
     filelist8=os.listdir(SEDfile_dir+'/phoenixm35') 
     filelist9=os.listdir(SEDfile_dir+'/phoenixm40') 
     filelist10=os.listdir(SEDfile_dir+'/phoenixp03') 
     filelist11=os.listdir(SEDfile_dir+'/phoenixp05') 
     
     filelist1_group = [os.path.join('phoenixm00',f) for f in filelist1 if f.endswith('.fits')]
     filelist2_group = [os.path.join('phoenixm05',f) for f in filelist2 if f.endswith('.fits')]
     filelist3_group = [os.path.join('phoenixm10',f) for f in filelist3 if f.endswith('.fits')]
     filelist4_group = [os.path.join('phoenixm15',f) for f in filelist4 if f.endswith('.fits')]
     filelist5_group = [os.path.join('phoenixm20',f) for f in filelist5 if f.endswith('.fits')]
     filelist6_group = [os.path.join('phoenixm25',f) for f in filelist6 if f.endswith('.fits')]
     filelist7_group = [os.path.join('phoenixm30',f) for f in filelist7 if f.endswith('.fits')]
     filelist8_group = [os.path.join('phoenixm35',f) for f in filelist8 if f.endswith('.fits')]
     filelist9_group = [os.path.join('phoenixm40',f) for f in filelist9 if f.endswith('.fits')]
     filelist10_group = [os.path.join('phoenixp03',f) for f in filelist10 if f.endswith('.fits')]
     filelist11_group = [os.path.join('phoenixp05',f) for f in filelist11 if f.endswith('.fits')]
     
     filelist_group=filelist1_group + filelist2_group + filelist3_group + filelist4_group + filelist5_group+ \
filelist6_group + filelist7_group + filelist8_group + filelist9_group +filelist10_group+filelist11_group
     
     fits_files=filelist_group
     
     
     obj_headers = []
     obj_files = []
     for filename in fits_files:
         index=0
         if re.search('fits',filename):  #example of filename filter
             index+=1
             fullfilename = os.path.join(SEDfile_dir,filename)
             hdr = fits.getheader(fullfilename)
             obj_headers.append(hdr)
             obj_files.append(filename)
     
        
     obj_temperatures = []
     obj_log_z_all = []
     index=0
     for hdr in obj_headers: 
        obj_temp=float(obj_headers[index]['TEFF'])
        obj_logz=float(obj_headers[index]['LOG_Z'])
        obj_temperatures.append(obj_temp)
        obj_log_z_all.append(obj_logz)
        index+=1   
      
        
     obj_names2 = []
     index=0
     for thefile in fits_files:
         #thenames=re.findall('^bk_([a-z][0-9]+).fits$',thefile)
         thenames=re.findall('([a-z].+_[0-9].+).fits$',thefile) 
         if(len(thenames)>0):
             obj_names2.append(thenames[0])
         else:
             print 'bad file ',thefile
         index+=1
        
     obj_names=obj_names2
     obj_files=filelist_group    
        
     #objames_and_objfiles = zip(obj_names, obj_files)
     #objames_and_objtemp = zip(obj_names, obj_temperatures)
     objtemp_and_objlogz = zip(obj_temperatures,obj_log_z_all)  
           
     #all_logg=np.array([0.0,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5])
     #all_logg=np.array([0.0,1.,2.,3.,4.])
     all_logg=np.array([0.0])
     
     index=0
     for temp,logz in objtemp_and_objlogz:
         if index%100==0:
             print 'phoenix star : T=',temp,' metal=', logz
         for logg in all_logg:      
             #Icat(model,temp,logz,logg)
             sed = S.Icat('phoenix', temp,logz,logg)    
             sed.convert('flam') # to be sure every spectrum is in flam unit
             all_sed.append(sed)
             index+=1
     
        
     return all_sed

#---------------------------------------------------------------------------------------------
def get_many_k93model():

    SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[8])
    all_sed=[]

    filelist1=os.listdir(SEDfile_dir+'/km01') 
    filelist2=os.listdir(SEDfile_dir+'/km02') 
    filelist3=os.listdir(SEDfile_dir+'/km03') 
    filelist4=os.listdir(SEDfile_dir+'/km05') 
    filelist5=os.listdir(SEDfile_dir+'/km10') 
    filelist6=os.listdir(SEDfile_dir+'/km20') 
    filelist7=os.listdir(SEDfile_dir+'/km25') 
    filelist8=os.listdir(SEDfile_dir+'/km30') 
    filelist9=os.listdir(SEDfile_dir+'/km35') 
    filelist10=os.listdir(SEDfile_dir+'/km40') 
    filelist11=os.listdir(SEDfile_dir+'/km45') 
    filelist12=os.listdir(SEDfile_dir+'/km50') 
    filelist13=os.listdir(SEDfile_dir+'/kp00') 
    filelist14=os.listdir(SEDfile_dir+'/kp01') 
    filelist15=os.listdir(SEDfile_dir+'/kp02') 
    filelist16=os.listdir(SEDfile_dir+'/kp03') 
    filelist17=os.listdir(SEDfile_dir+'/kp05') 
    filelist18=os.listdir(SEDfile_dir+'/kp10') 


    filelist1.remove('AA_README')
    filelist2.remove('AA_README')
    filelist3.remove('AA_README')
    filelist4.remove('AA_README')
    filelist5.remove('AA_README')
    filelist6.remove('AA_README')
    filelist7.remove('AA_README')
    filelist8.remove('AA_README')
    filelist9.remove('AA_README')
    filelist10.remove('AA_README')
    filelist11.remove('AA_README')
    filelist12.remove('AA_README')
    filelist13.remove('AA_README')
    filelist14.remove('AA_README')
    filelist15.remove('AA_README')
    filelist16.remove('AA_README')
    filelist17.remove('AA_README')
    filelist18.remove('AA_README')
    
    filelist=filelist1 + filelist2 + filelist3 + filelist4 + filelist5+ filelist6 + filelist7 + filelist8 + filelist9 + \
filelist10 + filelist11 + filelist12 + filelist13 + filelist14 + filelist15+ filelist16 + filelist17 + filelist18 

    filelist1_group = [os.path.join('km01',f) for f in filelist1 if f.endswith('.fits')]
    filelist2_group = [os.path.join('km02',f) for f in filelist2 if f.endswith('.fits')]
    filelist3_group = [os.path.join('km03',f) for f in filelist3 if f.endswith('.fits')]
    filelist4_group = [os.path.join('km05',f) for f in filelist4 if f.endswith('.fits')]
    filelist5_group = [os.path.join('km10',f) for f in filelist5 if f.endswith('.fits')]
    filelist6_group = [os.path.join('km20',f) for f in filelist6 if f.endswith('.fits')]
    filelist7_group = [os.path.join('km25',f) for f in filelist7 if f.endswith('.fits')]
    filelist8_group = [os.path.join('km30',f) for f in filelist8 if f.endswith('.fits')]
    filelist9_group = [os.path.join('km35',f) for f in filelist9 if f.endswith('.fits')]
    filelist10_group = [os.path.join('km40',f) for f in filelist10 if f.endswith('.fits')]
    filelist11_group = [os.path.join('km45',f) for f in filelist11 if f.endswith('.fits')]
    filelist12_group = [os.path.join('km50',f) for f in filelist12 if f.endswith('.fits')]
    filelist13_group = [os.path.join('kp00',f) for f in filelist13 if f.endswith('.fits')]
    filelist14_group = [os.path.join('kp01',f) for f in filelist14 if f.endswith('.fits')]
    filelist15_group = [os.path.join('kp02',f) for f in filelist15 if f.endswith('.fits')]
    filelist16_group = [os.path.join('kp03',f) for f in filelist16 if f.endswith('.fits')]
    filelist17_group = [os.path.join('kp05',f) for f in filelist17 if f.endswith('.fits')]
    filelist18_group = [os.path.join('kp10',f) for f in filelist18 if f.endswith('.fits')]
    
    filelist_group=filelist1_group + filelist2_group + filelist3_group + filelist4_group + filelist5_group+ \
filelist6_group + filelist7_group + filelist8_group + filelist9_group + filelist10_group + filelist11_group + filelist12_group + filelist5_group+ \
filelist13_group + filelist14_group + filelist15_group + filelist16_group + filelist17_group + filelist18_group
    
    fits_files=filelist_group
    
    obj_headers = []
    obj_files = []
    for filename in fits_files:
        index=0
        if re.search('fits',filename):  #example of filename filter
            index+=1
            fullfilename = os.path.join(SEDfile_dir,filename)
            hdr = fits.getheader(fullfilename)
            obj_headers.append(hdr)
            obj_files.append(filename)
    
    obj_temperatures = []
    obj_log_z_all = []
    index=0
    for hdr in obj_headers: 
        obj_temp=obj_headers[index]['TEFF']
        obj_logz=obj_headers[index]['LOG_Z']
        obj_temperatures.append(obj_temp)
        obj_log_z_all.append(obj_logz)
        index+=1
 
    obj_names2 = []
    index=0
    obj_temp = []

    for thefile in fits_files:
        #thenames=re.findall('^bk_([a-z][0-9]+).fits$',thefile)
        thenames=re.findall('([a-z].+_[0-9].+).fits$',thefile) 
        temp=re.findall('([a-z].+_[0-9].+).fits$',thefile) 
        if(len(thenames)>0):
            obj_names2.append(thenames[0])
        else:
            print 'bad file ',thefile
        index+=1
    
    
    obj_names=obj_names2
    obj_files=filelist_group
    
    objames_and_objfiles = zip(obj_names, obj_files)
    objames_and_objtemp = zip(obj_names, obj_temperatures)
    objtemp_and_objlogz = zip(obj_temperatures,obj_log_z_all)
    #all_logg=np.array([0.0,1.,2.,3.,4.])  # gravity
    all_logg=np.array([0.0])  # gravity
    
    for temp,logz in objtemp_and_objlogz:
        #Icat(model,temp,logz,logg)
        for logg in all_logg:
            sed = S.Icat('k93models', temp,logz,logg) 
            sed.convert('flam') # to be sure every spectrum is in flam unit
            all_sed.append(sed)
    
    return all_sed    
#------------------------------------------------------------------------------------------------------    
def get_many_bkmodels():
    
    all_sed=[]  
    
    SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[10])
    filelist=os.listdir(SEDfile_dir) 
    fits_files = [f for f in os.listdir(SEDfile_dir) if f.endswith('.fits')]
    
    obj_headers = []
    obj_files = []
    for filename in fits_files:
        index=0
        if re.search('fits',filename):  #example of filename filter
            index+=1
            fullfilename = os.path.join(SEDfile_dir,filename)
            hdr = fits.getheader(fullfilename)
            obj_headers.append(hdr)
            obj_files.append(filename)
            
            
    obj_names = []
    index=0
    for hdr in obj_headers: 
        obj_name=obj_headers[index]['TARGETID']
        obj_names.append(obj_name)
        index+=1
        
    obj_names2 = []
    index=0
    for thefile in fits_files:
        thenames=re.findall('^bk_([a-z][0-9]+).fits$',thefile)
        if(len(thenames)>0):
            obj_names2.append(thenames[0])
        else:
            print 'bad file ',thefile
        index+=1
        
    obj_names=obj_names2
    
    
    objames_and_objfiles = zip(obj_names, obj_files)
    
    OBJDict= {}
    for obj,thefile in objames_and_objfiles:
        #print obj,': '
        OBJDict[obj]=thefile
    
    
    for keyobj in OBJDict:
        the_file=OBJDict[keyobj]
        
        selected_file=the_file
        selected_fullfile=os.path.join(SEDfile_dir,selected_file)
        
        sed=S.FileSpectrum(selected_fullfile)
        sed.convert('flam') # to be sure every spectrum is in flam unit
        all_sed.append(sed)
    
    return all_sed    
#-------------------------------------------------------------------------------------------------
def get_many_bpgs():

     all_sed=[]  
     SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[1]) 
     fits_files = [f for f in os.listdir(SEDfile_dir) if f.endswith('.fits')]
     
     obj_headers = []
     obj_files = []
     for filename in fits_files:
         index=0
         if re.search('fits',filename):  #example of filename filter
             index+=1
             fullfilename = os.path.join(SEDfile_dir,filename)
             hdr = fits.getheader(fullfilename)
             obj_headers.append(hdr)
             obj_files.append(filename)
             
     obj_names = []
     
     index=0
     for hdr in obj_headers: 
         obj_name=obj_headers[index]['TARGETID']
         obj_names.append(obj_name)
         index+=1
        
     obj_names2 = []
     index=0
     for thefile in fits_files:
        #thenames=re.findall('^bk_([a-z][0-9]+).fits$',thefile)
        thenames=re.findall('^bpgs_([0-9].*).fits$',thefile)
        if(len(thenames)>0):
            obj_names2.append('bpgs_'+thenames[0])
        else:
            print 'bad file ',thefile
        index+=1
     
     obj_names=obj_names2
     
     objames_and_objfiles = zip(obj_names, obj_files)
     
     OBJDict= {}
     for obj,thefile in objames_and_objfiles:
         OBJDict[obj]=thefile
   
     for keyobj in OBJDict:
         the_file=OBJDict[keyobj]
        
         selected_file=the_file
         selected_fullfile=os.path.join(SEDfile_dir,selected_file)
        
         sed=S.FileSpectrum(selected_fullfile)
     
         sed.convert('flam') # to be sure every spectrum is in flam unit
         all_sed.append(sed)
         
         
     return all_sed    
#-----------------------------------------------------------------------------------------------    
    
    

#----------------------------------------------------------------------------------------------    
def get_all_thermalbb_flatT(N=100,TMIN=3000.,TMAX=50000.):
    
    all_sed=[]   
    
    all_random_values=np.random.random_sample((N,))
    all_random_Temp=TMIN+(TMAX-TMIN)*all_random_values
    for index in np.arange(N):
        T=all_random_Temp[index]
        sed=S.BlackBody(T)
        sed.convert('flam') # to be sure every spectrum is in flam unit
        all_sed.append(sed)
        
    return all_sed


#---------------------------------------------------------------------------------
def plot_allsed(all_sed,thetitle,figfilename,yscale='lin',XMIN=3000.,XMAX=10000.,YMIN=0,YMAX=0):
    plt.figure()
    
    for sed in all_sed:  
        
        sed.convert('flam') # to be sure every spectrum is in flam unit
        
        if yscale=='log':
            plt.semilogy(sed.wave,sed.flux,lw=2)
        else:
            plt.plot(sed.wave,sed.flux,lw=2)

    plt.xlim(XMIN,XMAX)
    if(YMAX>0):
        plt.ylim(YMIN,YMAX)
        
    plt.xlabel(sed.waveunits)
    plt.ylabel(sed.fluxunits)
    plt.grid(True)
    plt.title(thetitle)
    plt.savefig(figfilename)
    
#------------------------------------------------------------------------------------------------
def get_all_bc95(z=0):
    SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[5],'templates')
    
    filelist=os.listdir(SEDfile_dir)
   
    
    obj_headers = []
    obj_files = []
    index=0
    for filename in filelist:  
        if re.search('fits',filename):  #example of filename filter
            index+=1
            fullfilename = os.path.join(SEDfile_dir,filename)
            hdr = fits.getheader(fullfilename)
            obj_headers.append(hdr)
            obj_files.append(filename)
            
    obj_names = []
    index=0
    for hdr in obj_headers: 
        obj_name=obj_headers[index]['TARGETID']
        obj_names.append(obj_name)
        index+=1
        
    objames_and_objfiles = zip(obj_names, obj_files)
    
    OBJDict= {}
    for obj,thefile in objames_and_objfiles:
        print obj,': '
        OBJDict[obj]=thefile
        print OBJDict[obj] 
    
    all_sed=[]        
    for keyobj in OBJDict:
        the_file=OBJDict[keyobj]
        
        selected_file=the_file
        selected_fullfile=os.path.join(SEDfile_dir,selected_file)
        
        sed=S.FileSpectrum(selected_fullfile)
        sed.convert('flam') # to be sure every spectrum is in flam unit
        if (z>0):
            sed_z=sed.redshift(z)
            all_sed.append(sed_z)
        else:
            all_sed.append(sed)
    return all_sed
    

#------------------------------------------------------------------------------------------------
def get_many_bc95():
    '''
    
    Get many Galaxy profiles
    '''
    
    
    SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[5],'templates')
    
    filelist=os.listdir(SEDfile_dir)
   
    
    obj_headers = []
    obj_files = []
    index=0
    for filename in filelist:  
        if re.search('fits',filename):  #example of filename filter
            index+=1
            fullfilename = os.path.join(SEDfile_dir,filename)
            hdr = fits.getheader(fullfilename)
            obj_headers.append(hdr)
            obj_files.append(filename)
            
    obj_names = []
    index=0
    for hdr in obj_headers: 
        obj_name=obj_headers[index]['TARGETID']
        obj_names.append(obj_name)
        index+=1
        
    objames_and_objfiles = zip(obj_names, obj_files)
    
    OBJDict= {}
    for obj,thefile in objames_and_objfiles:
        #print obj,': '
        OBJDict[obj]=thefile
        #print OBJDict[obj] 
    
    all_sed=[]  
    all_z_rs=np.linspace(0,2.,10)
      
    #for keyobj in OBJDict:
    for obj,thefile in objames_and_objfiles:
        #the_file=OBJDict[keyobj]
        #selected_file=the_file
        selected_file=thefile
        selected_fullfile=os.path.join(SEDfile_dir,selected_file)
        
        sed=S.FileSpectrum(selected_fullfile)
        sed.convert('flam') # to be sure every spectrum is in flam unit

        # loop on Redshifts
        for z_rs in all_z_rs:
            if (z_rs>0):
                sed_z=sed.redshift(z_rs)
                all_sed.append(sed_z)
            else:
                all_sed.append(sed)
                
                
                
    return all_sed

 #------------------------------------------------------------------------------------------------   
def get_all_kc96(z=0):
    SEDfile_dir=os.path.join(top_pysynphot_data_dir,dir_nostar,dir_submodels[13])
    #filelist=os.listdir(SEDfile_dir) 
    
    fits_files = [f for f in os.listdir(SEDfile_dir) if f.endswith('.fits')]
    
    obj_headers = []
    obj_files = []
    for filename in fits_files:
        index=0
        if re.search('fits',filename):  #example of filename filter
            index+=1
            fullfilename = os.path.join(SEDfile_dir,filename)
            hdr = fits.getheader(fullfilename)
            obj_headers.append(hdr)
            obj_files.append(filename)
      
    
    obj_names = []
    index=0
    for hdr in obj_headers: 
        obj_name=obj_headers[index]['TARGETID']
        obj_names.append(obj_name)
        index+=1
        
    objames_and_objfiles = zip(obj_names, obj_files)
    
    OBJDict= {}
    for obj,thefile in objames_and_objfiles:
        print obj,': '
        OBJDict[obj]=thefile
        print OBJDict[obj] 
    
    all_sed=[]       
    for keyobj in OBJDict:
        the_file=OBJDict[keyobj]      
        selected_file=the_file
        selected_fullfile=os.path.join(SEDfile_dir,selected_file)      
        sed=S.FileSpectrum(selected_fullfile)
        sed.convert('flam') # to be sure every spectrum is in flam unit
        
        
        if (z>0):
            sed_z=sed.redshift(z)
            all_sed.append(sed_z)
        else:
            all_sed.append(sed)
        
    return all_sed
        

#------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    print 'hello'

    #all_sed=get_all_ck04models()
    
    #plot_allsed_starmodels(all_sed,"ck04models","ck04models.png")
    
    Flag_CALSPEC_HD=False
    Flag_BC95_Z3=False
    Flag_BC95=True
    Flag_KC93_Z0=False
    Flag_THERMALBB=False
    Flag_PHOENIX=False
    Flag_CK04=False
    Flag_PICKLE=False
    Flag_K93=False
    Flag_BK=False
    Flag_BPGS=False
    
    if Flag_CALSPEC_HD:
        all_sed=get_all_calspec_hd()
        plot_allsed(all_sed,'SED of CALSPEC stars','calspec_hd_lin.png',yscale='lin',YMIN=0,YMAX=0.4e-10)
        plot_allsed(all_sed,'SED of CALSPEC stars','calspec_hd_log.png',yscale='log',YMIN=1e-14,YMAX=1e-10)

    
    if Flag_BC95_Z3:
        all_sed=get_all_bc95(z=0)
        plot_allsed(all_sed,'SED of Bruzaual-Charlot Atlas (bc95 galaxies)','gal_bc95_lin.png',XMIN=0,XMAX=11000.,yscale='lin')
        plot_allsed(all_sed,'SED of Bruzaual-Charlot Atlas (bc95 galaxies)','gal_bc95_log.png',XMIN=0,XMAX=11000.,yscale='log')
        
    if Flag_BC95:
        all_sed=get_many_bc95()
        plot_allsed(all_sed,'SED of Bruzaual-Charlot Atlas (bc95 galaxies)','gal_bc95_lin.png',XMIN=0,XMAX=11000.,yscale='lin')
        plot_allsed(all_sed,'SED of Bruzaual-Charlot Atlas (bc95 galaxies)','gal_bc95_log.png',XMIN=0,XMAX=11000.,yscale='log')
    
    
    if Flag_KC93_Z0:
        all_sed=get_all_kc96(z=0)
        plot_allsed(all_sed,'SED of Kinney-Calzetti Atlas (kc96 galaxies)','gal_kc96_lin.png',yscale='lin')
        plot_allsed(all_sed,'SED of Kinney-Calzetti Atlas (kc96 galaxies)','gal_kc96_log.png',yscale='log')
    
    if Flag_THERMALBB:
        all_sed=get_all_thermalbb_flatT()
        plot_allsed(all_sed,'SED of thermal BB','star_thermalbb_lin.png',yscale='lin')
        plot_allsed(all_sed,'SED of thermal BB','star_thermalbb_log.png',yscale='log')
    
    if Flag_PHOENIX:
        all_sed=get_all_phoenix()
        plot_allsed(all_sed,'SED of Phoenix stars','star_phoenix_lin.png',yscale='lin')
        plot_allsed(all_sed,'SED of Phoenix stars','star_phoenix_log.png',yscale='log')
        
    if Flag_CK04:
        all_sed=get_many_ck04models()
        plot_allsed(all_sed,'SED of ck04 stars','star_ck04_lin.png',yscale='lin')
        plot_allsed(all_sed,'SED of ck04 stars','star_ck04_log.png',yscale='log')
        
    if Flag_PICKLE:  
        all_sed=get_all_pickle()
        plot_allsed(all_sed,'SED of pickle stars','star_pickle_lin.png',yscale='lin')
        plot_allsed(all_sed,'SED of pickle stars','star_pickle_log.png',yscale='log')
        
    if Flag_K93:
        all_sed=get_many_k93model()
        plot_allsed(all_sed,'SED of k93 stars','star_k93_lin.png',yscale='lin')
        plot_allsed(all_sed,'SED of k93 stars','star_k93_log.png',yscale='log')
        
    if Flag_BK:
        all_sed=get_many_bkmodels()
        plot_allsed(all_sed,'SED of bk stars','star_bk_lin.png',yscale='lin')
        plot_allsed(all_sed,'SED of bk stars','star_bk_log.png',yscale='log')
        
    if Flag_BPGS:    
        all_sed=get_many_bpgs()
        plot_allsed(all_sed,'SED of bpgs stars','star_bpgs_lin.png',yscale='lin')
        plot_allsed(all_sed,'SED of bpgs stars','star_bpgs_log.png',yscale='log')
        
    