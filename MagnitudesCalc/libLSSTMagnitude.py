"""
libLSSTFilter
=============


author : Sylvie Dagoret-Campagne
affiliation : LAL/CNRS/IN2P3/FRANCE
Collaboration : DESC-LSST


Purpose : study of atmosphere on LSST Magnitudes

"""


import os
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
import os
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

#---------------------------------------------------------------------------------
# DATA BLOCK 
detdatafilename='data/transmissions-LSST.dat'
#various path to be used
path_atm_rt_us_sa_rt_oz='/Users/dagoret-campagnesylvie/MacOsX/LSST/MyWork/GitHub/PC5AtmosphericExtinction/LibRadTran/simulations/RT/2.0/LS/pp/us/sa/rt/oz/out'
path_atm_rt_us_sa_lt_oz='/Users/dagoret-campagnesylvie/MacOsX/LSST/MyWork/GitHub/PC5AtmosphericExtinction/LibRadTran/simulations/RT/2.0/LS/pp/us/sa/lt/oz/out'
#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------
class LSST_Magnitude:
    """
    class LSST_Magnitude
    ====================
    
    Class to compute relative magnitude
    
    """
    def ComputeSEDxAtm(self,wl_sed,sed,wl_atm,atm):
        interpol_atm=interp1d(wl_atm,atm)
        return sed*interpol_atm(wl_sed)*wl_sed
        
    def ComputeSEDxAtmxFilt(self,wl_sed,sed,wl_atm,atm,wl_filt,filt):
        interpol_atm=interp1d(wl_atm,atm)
        interpol_filt=interp1d(wl_filt,filt)
        return sed*interpol_atm(wl_sed)*interpol_filt(wl_sed)*wl_sed 
        
    def ComputeRelativeMagSEDxAtm(self,wl_sed,sed,wl_atm,atm):
        interpol_atm=interp1d(wl_atm,atm)
        theproduct=sed*interpol_atm(wl_sed)*wl_sed
        thesum=theproduct.sum()
        themag=2.5*np.log10(thesum)
        return themag
    def ComputeRelativeMagSEDxAtmxFilt(self,wl_sed,sed,wl_atm,atm,wl_filt,filt):
        interpol_atm=interp1d(wl_atm,atm)
        interpol_filt=interp1d(wl_filt,filt)
        theproduct=sed*interpol_atm(wl_sed)*interpol_filt(wl_sed)*wl_sed 
        thesum=theproduct.sum()
        themag=2.5*np.log10(thesum)
        return themag



#----------------------------------------------------------------------------
class RT_Atmosphere:
    """
    class RT_Atmosphere
    ================
    
    class to read libradtran atmosphere
      
    """
    pathname=""
    anyfilelist = []
    atmfilelist = []
    def __init__(self,pathname):
        self.pathname=pathname
        for file in os.listdir(self.pathname):
            self.anyfilelist.append(file)
            if re.search('.OUT',file) or re.search('.out',file) :
                self.atmfilelist.append(file)
    def list_of_anyfiles(self):
        return self.anyfilelist
    def list_of_atmfiles(self):
        return self.atmfilelist
    def get_air_transparency(self,filename):
        if self.atmfilelist.count(filename) == 1:
             fullpath=os.path.join(self.pathname,filename)
             data = np.loadtxt(fullpath)
             x=data[:,0]
             y=data[:,1]
             return x,y
        else:
            print "No atmospheric file ", filename,"found"
            return 0,0
#-----------------------------------------------------------------------------------       
    
    
    
    
#---------------------------------------------------------------------------
class Filter:
    """
    class Filter
    =============
    Variables :
    ----------
     wl : wavelength 1d array
     u  : u filter 1d array
     g  : g filter 1d array
     r  : r filter 1d array
     i  : i filter 1d array
     z  : z filter 1d array
     y4 : y4 filter 1d array
     topt
     tccd
     atm
    Methods :
    ---------
    
    """
    wl=0
    u=0
    g=0
    r=0
    i=0
    z=0
    y4=0
    topt=0
    tccd=0
    atm=0
    
    def __init__(self):
        df=pd.read_csv(detdatafilename,names=['wl','Topt','Tccd','U','G','R','I','Z','Y4','atm'],sep='\t')
        self.wl=np.asarray(df['wl'])
        self.u=np.asarray(df['U'])*0.01
        self.g=np.asarray(df['G'])*0.01
        self.r=np.asarray(df['R'])*0.01
        self.i=np.asarray(df['I'])*0.01
        self.z=np.asarray(df['Z'])*0.01
        self.y4=np.asarray(df['Y4'])*0.01
        self.atm=np.asarray(df['atm'])
        print 'init Filter size=',self.wl.shape
    def wavelength_to_u_spl(self):
        #return UnivariateSpline(self.wl,self.u)
        return interp1d(self.wl,self.u,kind='cubic')
    def wavelength_to_g_spl(self):
        #return UnivariateSpline(self.wl,self.g)   
        return interp1d(self.wl,self.g,kind='cubic') 
    def wavelength_to_r_spl(self):
        return interp1d(self.wl,self.r,kind='cubic')  
    def wavelength_to_i_spl(self):
        return interp1d(self.wl,self.i,kind='cubic')          
    def wavelength_to_z_spl(self):
        return interp1d(self.wl,self.z,kind='cubic')
    def wavelength_to_y4_spl(self):
        return interp1d(self.wl,self.y4,kind='cubic')
    def get_u_tr(self):
        return self.wl,self.u
    def get_g_tr(self):
        return self.wl,self.g
    def get_r_tr(self):
        return self.wl,self.r
    def get_i_tr(self):
        return self.wl,self.i
    def get_z_tr(self):
        return self.wl,self.z
    def get_y4_tr(self):
        return self.wl,self.y4
        
#-------------------------------------------------------------------------------------    


#---------------------------------------------------------------------------------
def MakeSED(lambda_min=300.,lambda_max=1200.,dlambda=1.,slope=0):
    """
     MakeSED(lambda_min=300,lambda_max=1200,dlambda=1,power=0)
     =============================================
     
     input :
     -----
         lambda_min : minimum of SED spectra in nm
         lambda_max : maximum of SED spectra in nm
         dlambda    : bin width of SED in nm
         power      : power in wavelength
     
     output :
     --------
         wl         : wavelength 1D array
         sed        : SED 1D array
    
    """
    
    #NBINS=(lambda_max-lambda_min=300.)/dlambda
    wl=np.arange(lambda_min,lambda_max+dlambda,dlambda)
    nbins=wl.shape[0]
    bincenter=nbins/2
    
    
    sed=np.power(wl,slope)
    sedcenter=sed[bincenter]
    sed=sed/sedcenter
        
    return wl,sed  
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
def MakeFilters():
    df=pd.read_csv(detdatafilename,names=['wl','Topt','Tccd','U','G','R','I','Z','Y4','atm'],sep='\t')

    return df['wl'],df['U']*0.01,df['G']*.01,df['R']*0.01,df['I']*0.01,df['Z']*0.01,df['Y4']*0.01
#---------------------------------------------------------------------------------


def PlotFilter():
    df=pd.read_csv(detdatafilename,names=['wl','Topt','Tccd','U','G','R','I','Z','Y4','atm'],sep='\t')
    df.head()
    plt.figure()
    colors = ['blue','green','red', 'orange','grey','black'] 
    df.plot(x='wl', y=['U','G','R','I','Z','Y4'],color=colors)
    plt.ylim([0,100])
    plt.xlabel("$\lambda$")
    plt.ylabel("Filter transmission")
    plt.title("Filters")
#----------------------------------------------------------------------------------    
def PlotSED():
     # create a SED
    wl,sed=MakeSED(lambda_min=300.,lambda_max=1200.,dlambda=1.,slope=-3)
    
    print wl.shape
    plt.figure()
    plt.plot(wl,sed,'-')
    
#----------------------------------------------------------------------------    
def PlotSESxFilter():
    
    # 1) SED
    wl,sed=MakeSED(lambda_min=300.,lambda_max=1200.,dlambda=10.,slope=-1)
    
    # 2) Filter
    flt=Filter()
    
    # get interpolation function
    wl_to_u= flt.wavelength_to_u_spl()
    wl_to_g= flt.wavelength_to_g_spl()
    wl_to_r= flt.wavelength_to_r_spl()
    wl_to_i= flt.wavelength_to_i_spl()
    wl_to_z= flt.wavelength_to_z_spl()
    wl_to_y4= flt.wavelength_to_y4_spl()
    
    tu=wl_to_u(wl)*sed
    tg=wl_to_g(wl)*sed
    tr=wl_to_r(wl)*sed
    ti=wl_to_i(wl)*sed
    tz=wl_to_z(wl)*sed
    ty4=wl_to_y4(wl)*sed
    
    plt.figure()
    plt.plot(wl,tu)
    plt.plot(wl,tg)
    plt.plot(wl,tr)
    plt.plot(wl,ti)
    plt.plot(wl,tz)
    plt.plot(wl,ty4)
#-----------------------------------------------------------------------------

def PlotAtmosphere():
    atm=RT_Atmosphere(path_atm_rt_us_sa_rt_oz)
    #print atm.list_of_atmfiles()
    wl,tr=atm.get_air_transparency("RT_LS_pp_us_sa_rt_z30_oz22.OUT")
    plt.figure()
    plt.plot(wl,tr)

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
if __name__ == "__main__":


    PlotSED()

    # Build the SED
    (wl_sed,sed)=MakeSED(lambda_min=300.,lambda_max=1200.,dlambda=1.,slope=-3)
    
    print 'wl_sed = ',wl_sed.shape
    print 'sed = ',sed.shape    
    
    # Build the atmosphere
    atm=RT_Atmosphere(path_atm_rt_us_sa_rt_oz)
    #print atm.list_of_atmfiles()
    (wl_atm,tr_atm)=atm.get_air_transparency("RT_LS_pp_us_sa_rt_z30_oz22.OUT")

    print 'wl_atm = ',wl_atm.shape
    print 'tr_atm = ',tr_atm.shape
    

    mag=LSST_Magnitude()    
    fl=mag.ComputeSEDxAtm(wl_sed,sed,wl_atm,tr_atm)
    
    #interpol_atm=interp1d(wl_atm,atm)
    #sed*interpol_atm(wl_sed)    
    
    plt.figure()
    plt.plot(wl_sed,fl)  
    
    flt=Filter()
    wl_i,i=flt.get_i_tr()
    fl2=mag.ComputeSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_i,i)
    
    plt.figure()
    plt.plot(wl_sed,fl2)  
    
    M1_i=mag.ComputeRelativeMagSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_i,i)
    print 'magnitude(I) = ' , M1
    