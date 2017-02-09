"""
libLSSTMagnitude.py
=============----


author : Sylvie Dagoret-Campagne
affiliation : LAL/CNRS/IN2P3/FRANCE
Collaboration : DESC-LSST

Update 2017/02/09


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
path_atm_rt_us_sa_rt_pwv='/Users/dagoret-campagnesylvie/MacOsX/LSST/MyWork/GitHub/PC5AtmosphericExtinction/LibRadTran/simulations/RT/2.0/LS/pp/us/sa/rt/pwv/out'
path_atm_rt_us_sa_lt_pwv='/Users/dagoret-campagnesylvie/MacOsX/LSST/MyWork/GitHub/PC5AtmosphericExtinction/LibRadTran/simulations/RT/2.0/LS/pp/us/sa/lt/pwv/out'
#
path_atm_rt_sw_sa_rt_oz='/Users/dagoret-campagnesylvie/MacOsX/LSST/MyWork/GitHub/PC5AtmosphericExtinction/LibRadTran/simulations/RT/2.0/LS/pp/sw/sa/rt/oz/out'
path_atm_rt_sw_sa_lt_oz='/Users/dagoret-campagnesylvie/MacOsX/LSST/MyWork/GitHub/PC5AtmosphericExtinction/LibRadTran/simulations/RT/2.0/LS/pp/sw/sa/lt/oz/out'
path_atm_rt_sw_sa_rt_pwv='/Users/dagoret-campagnesylvie/MacOsX/LSST/MyWork/GitHub/PC5AtmosphericExtinction/LibRadTran/simulations/RT/2.0/LS/pp/sw/sa/rt/pwv/out'
path_atm_rt_sw_sa_lt_pwv='/Users/dagoret-campagnesylvie/MacOsX/LSST/MyWork/GitHub/PC5AtmosphericExtinction/LibRadTran/simulations/RT/2.0/LS/pp/sw/sa/lt/pwv/out'
#

modtran_path='modtran_samples/MT_FirstSamples'
modtran_atmfile="Pachon_MODTRAN.1.5.kg.1.6.17.xlsx"

path_sed='./SED'
filename_star=['Star.B1.No.3.xlsx']
filename_gal=['Gal.S0.template.xlsx','Gal.BC.95.No.1.xlsx','Gal.GS.39.No.2.xlsx']
filename_pick_uk_xcl=['Pick.UK.No.2.22.xlsx','Pick.UK.50.No.4.xlsx']
filename_pick_uk_fits=['pickles_uk_22.fits','pickles_uk_50.fits']    # fits are the same file as above xcl file
filename_pick_110=['Pickles.No.1.110.xlsx','pickles_110.fits','pickles_110.ascii.txt'] # same files

#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------
class LSST_Magnitude:
    """
    class LSST_Magnitude
    ====================
    
    Class to compute relative magnitude
    
    """
    thesum=0
    themag=0
    
    def ComputeSEDxAtm(self,wl_sed,sed,wl_atm,atm):
        interpol_atm=interp1d(wl_atm,atm)
        return sed*interpol_atm(wl_sed)*wl_sed
        
    def ComputeSEDxAtmxFilt(self,wl_sed,sed,wl_atm,atm,wl_filt,filt):
        interpol_atm=interp1d(wl_atm,atm)
        interpol_filt=interp1d(wl_filt,filt)
        return sed*interpol_atm(wl_sed)*interpol_filt(wl_sed)*wl_sed 
        
    def ComputeRelativeMagSEDxAtm(self,wl_sed,sed,wl_atm,atm):
        interpol_atm=interp1d(wl_atm,atm)
        dwl=self.ComputeWlbinSize(wl_sed)
        theproduct=sed*interpol_atm(wl_sed)*wl_sed*dwl
        thesum=theproduct.sum()
        themag=-2.5*np.log10(thesum)
        return -themag
    def ComputeRelativeMagSEDxAtmxFilt(self,wl_sed,sed,wl_atm,atm,wl_filt,filt):
        dwl=self.ComputeWlbinSize(wl_sed)
        interpol_atm=interp1d(wl_atm,atm)
        interpol_filt=interp1d(wl_filt,filt)
        theproduct=sed*interpol_atm(wl_sed)*interpol_filt(wl_sed)*wl_sed*dwl 
        self.thesum=theproduct.sum()
        self.themag=-2.5*np.log10(self.thesum)
        return self.themag
    def ComputeWlbinSize(self,wl):
        len=wl.shape[0]
        wl_shift_right=np.roll(wl,1)
        wl_shift_left=np.roll(wl,-1)
        wl_bin_size=(wl_shift_left-wl_shift_right)/2. # size of each bin      
        wl_bin_size[0]=wl_bin_size[1]  # erase bin width
        wl_bin_size[len-1]=wl_bin_size[len-2]
        return wl_bin_size



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
    filename_ =''
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
    
class MT_Atmosphere:
    """
    class MT Atmosphere
    ===================
    
    class to read Modtran file
    """
    pathname=""
    anyfilelist = []
    atmfilelist = []
    filename_ =''
    def __init__(self,pathname):
        self.pathname=pathname
        for file in os.listdir(self.pathname):
            self.anyfilelist.append(file)
            if re.search('.xlsx',file) or re.search('.XLSX',file) :
                self.atmfilelist.append(file)
    def list_of_anyfiles(self):
        return self.anyfilelist
    def list_of_atmfiles(self):
        return self.atmfilelist
    def get_air_transparency(self,filename):
        if self.atmfilelist.count(filename) == 1:
             fullpath=os.path.join(self.pathname,filename)
             mtfile = pd.ExcelFile(fullpath)
             sheet_name=mtfile.sheet_names[0]
             df_colname = mtfile.parse(sheet_name,index_row=14,usecols=range(0,6))
             df = mtfile.parse(sheet_name,header=16,usecols=range(0,6))
             df.columns = ["wl", "comb","h2o","o2", "o3","scat"]
             MT_X=df["wl"]
             MT_Y1=df["h2o"]
             MT_Y2=df["o3"]
             MT_Y3=df["scat"]
             MT_Y4=df["o2"]
             MT_Y5=df["comb"]
             MT_Y6=MT_Y1*MT_Y2*MT_Y4  # 
             MT_Y=MT_Y1*MT_Y2*MT_Y3*MT_Y4
             
             x=MT_X
             y=MT_Y
             return x,y
        else:
             print "No modtran atmospheric file ", filename,"found"
             return 0,0				
    
    
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

def ReadSED(filename,filetype='xlsx',lambda_min=300.,lambda_max=1199.,dlambda=1.,headerstop=0,cmin=0,cmax=2):
    """
    ReadSED in a file
    ===================
    
    """
    
    x=np.zeros(1)
    y=np.zeros(1)
    
    if filetype=="xlsx":
        file_xlsx= pd.ExcelFile(filename)
        sheet_name = file_xlsx.sheet_names[0]
        df = file_xlsx.parse(sheet_name,headerstop,usecols=range(cmin,cmax))
        df.columns = ["wl","flux"]
        x=np.array(df["wl"])/10.
        y=np.array(df["flux"])
    elif filetype == 'fits' :
        hdulist = fits.open(filename)
        tabledata = hdulist[1].data
        x=np.array(tabledata.field('WAVELENGTH'))/10.
        y=np.array(tabledata.field('FLUX'))
    elif filetype == 'ascii' :
        arr=np.loadtxt(filename,skiprows=headerstop)
        x=arr[:,0]/10.
        y=arr[:,1]
    else:
        print "Unknown file"
    
    sel =np.where(np.logical_and(x>=lambda_min,x<lambda_max))
    
    xcut=x[sel]
    ycut=y[sel]

    return xcut,ycut

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
    plt.title('filters')
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
    plt.title('SED')
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('sed')
    
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
    
    plt.title('filter')
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('filter transmission')
#-----------------------------------------------------------------------------

def PlotAtmosphere():
    atm=RT_Atmosphere(path_atm_rt_us_sa_rt_oz)
    #print atm.list_of_atmfiles()
    wl,tr=atm.get_air_transparency("RT_LS_pp_us_sa_rt_z30_oz22.OUT")
    plt.figure()
    plt.plot(wl,tr)
    plt.title('atmosphere')
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('atmosphere transmission')

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
if __name__ == "__main__":

    # 1) Example of SED
    # -------------------
    PlotSED()

    # Build the SED
    # ---------------
    (wl_sed,sed)=MakeSED(lambda_min=300.,lambda_max=1100.,dlambda=1.,slope=-3)
    
    
    # 2) Build the atmosphere
    # -----------------------
    atm=RT_Atmosphere(path_atm_rt_us_sa_rt_oz)
    #print atm.list_of_atmfiles()
    (wl_atm,tr_atm)=atm.get_air_transparency("RT_LS_pp_us_sa_rt_z30_oz22.OUT")

    
    # 3) Filter
    # --------------
    flt=Filter()
    
    wl_u,u=flt.get_u_tr()
    wl_g,g=flt.get_g_tr()
    wl_r,r=flt.get_r_tr()
    wl_i,i=flt.get_i_tr()
    wl_z,z=flt.get_z_tr()
    wl_y4,y4=flt.get_y4_tr()
    
    # 4) Magnitude
    #--------------
    
    mag=LSST_Magnitude()    
    fl=mag.ComputeSEDxAtm(wl_sed,sed,wl_atm,tr_atm)
    
    #interpol_atm=interp1d(wl_atm,atm)
    #sed*interpol_atm(wl_sed)    
    
    plt.figure()
    plt.plot(wl_sed,fl)
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('spectrum')
    plt.title('Source spectrum through LibRadtran atmosphere')
    
    fl_u=mag.ComputeSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_u,u)
    fl_g=mag.ComputeSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_g,g)
    fl_r=mag.ComputeSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_r,r)
    fl_i=mag.ComputeSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_i,i)
    fl_z=mag.ComputeSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_z,z)
    fl_y4=mag.ComputeSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_y4,y4)
    
    plt.figure()
    plt.plot(wl_sed,fl_u) 
    plt.plot(wl_sed,fl_g)
    plt.plot(wl_sed,fl_r)
    plt.plot(wl_sed,fl_i) 
    plt.plot(wl_sed,fl_z)
    plt.plot(wl_sed,fl_y4) 
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('spectrum')
    plt.title('Source spectrum through Radtran atmosphere and Filters')
    
    M_u=mag.ComputeRelativeMagSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_u,u)
    print 'magnitude(U) = ' , M_u
    M_g=mag.ComputeRelativeMagSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_g,g)
    print 'magnitude(G) = ' , M_g
    M_r=mag.ComputeRelativeMagSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_r,r)
    print 'magnitude(R) = ' , M_r
    M_i=mag.ComputeRelativeMagSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_i,i)
    print 'magnitude(I) = ' , M_i
    M_z=mag.ComputeRelativeMagSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_z,z)
    print 'magnitude(Z) = ' , M_z
    M_y4=mag.ComputeRelativeMagSEDxAtmxFilt(wl_sed,sed,wl_atm,tr_atm,wl_y4,y4)
    print 'magnitude(Y4) = ' , M_y4
   
   
   #-------------------------------------------------------------------------
    
    # 5) Modtran
    mt_atm=MT_Atmosphere(modtran_path)
    print mt_atm.list_of_atmfiles()
    (mt_wl_atm,mt_tr_atm)=mt_atm.get_air_transparency(modtran_atmfile)
  
    
     
    mt_fl=mag.ComputeSEDxAtm(wl_sed,sed,mt_wl_atm,mt_tr_atm)
    
    plt.figure()
    plt.plot(wl_sed,mt_fl)
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('spectrum')
    plt.title('Source spectrum through Modtran atmosphere')
    
    
    
    
    mt_fl_u=mag.ComputeSEDxAtmxFilt(wl_sed,sed,mt_wl_atm,mt_tr_atm,wl_u,u)
    mt_fl_g=mag.ComputeSEDxAtmxFilt(wl_sed,sed,mt_wl_atm,mt_tr_atm,wl_g,g)
    mt_fl_r=mag.ComputeSEDxAtmxFilt(wl_sed,sed,mt_wl_atm,mt_tr_atm,wl_r,r)
    mt_fl_i=mag.ComputeSEDxAtmxFilt(wl_sed,sed,mt_wl_atm,mt_tr_atm,wl_i,i)
    mt_fl_z=mag.ComputeSEDxAtmxFilt(wl_sed,sed,mt_wl_atm,mt_tr_atm,wl_z,z)
    mt_fl_y4=mag.ComputeSEDxAtmxFilt(wl_sed,sed,mt_wl_atm,mt_tr_atm,wl_y4,y4)
    
    plt.figure()
    plt.plot(wl_sed,mt_fl_u) 
    plt.plot(wl_sed,mt_fl_g)
    plt.plot(wl_sed,mt_fl_r)
    plt.plot(wl_sed,mt_fl_i) 
    plt.plot(wl_sed,mt_fl_z)
    plt.plot(wl_sed,mt_fl_y4) 
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('spectrum')
    plt.title('Source spectrum through Modtran atmosphere and Filters')
    
    
    #-----------------------------------------------------------------------
    #------------------------------------------------------------------------
    
    
    ## Star
    ## ------
    fullfilename= os.path.join(path_sed,filename_star[0])   
    x_star,y_star=ReadSED(fullfilename)
    plt.figure()
    plt.title(filename_star[0])
    plt.plot(x_star,y_star)
    plt.xlabel("$\lambda$ (nm) ")
    plt.ylabel("flux")
    
    ## Galaxy
    ## -------
    fullfilename_gal= [  os.path.join(path_sed,filename_gal[j])  for j in range(0,3) ]

    x_gal_1,y_gal_1=ReadSED(fullfilename_gal[0])
    plt.figure()
    plt.title(filename_gal[0])
    plt.plot(x_gal_1,y_gal_1)
    plt.xlabel("$\lambda$ (nm) ")
    plt.ylabel("flux")
    
    x_gal_3,y_gal_3=ReadSED(fullfilename_gal[2])
    plt.figure()
    plt.title(filename_gal[2])
    plt.plot(x_gal_3,y_gal_3)
    plt.xlabel("$\lambda$ (nm) ")
    plt.ylabel("flux")
    
    x_gal_2,y_gal_2=ReadSED(fullfilename_gal[1],headerstop=22)
    plt.figure()
    plt.title(filename_gal[1])
    plt.plot(x_gal_2,y_gal_2)
    plt.xlabel("$\lambda$ (nm) ")
    plt.ylabel("flux")
    
    ## Pickle excl
    ## -----------
    fullfilename_pick_uk_xcl= [os.path.join(path_sed,filename_pick_uk_xcl[j])  for j in range(0,2) ]
    print fullfilename_pick_uk_xcl
    x_gal_1,y_gal_1=ReadSED(fullfilename_pick_uk_xcl[0],headerstop=38)
    x_gal_2,y_gal_2=ReadSED(fullfilename_pick_uk_xcl[1],headerstop=38)
    plt.figure()
    title = filename_pick_uk_xcl[0] + ' , ' + filename_pick_uk_xcl[1]
    plt.title(title)
    plt.plot(x_gal_1,y_gal_1)
    plt.plot(x_gal_2,y_gal_2)
    plt.xlabel("$\lambda$ (nm) ")
    plt.ylabel("flux")
    
    
    ## Pickle fits
    ## -----------
    fullfilename_pick_uk_fits= [os.path.join(path_sed,filename_pick_uk_fits[j])  for j in range(0,2) ]
    x_gal_3,y_gal_3=ReadSED(fullfilename_pick_uk_fits[0],filetype='fits')
    x_gal_4,y_gal_4=ReadSED(fullfilename_pick_uk_fits[1],filetype='fits')
    plt.figure()
    title = filename_pick_uk_fits[0] + ' , ' + filename_pick_uk_fits[1]
    plt.title(title)
    plt.plot(x_gal_3,y_gal_3)
    plt.plot(x_gal_4,y_gal_4)
    plt.xlabel("$\lambda$ (nm) ")
    plt.ylabel("flux")
    
    
    ## Final obj
    ## ---------
    fullfilename_pick_110= [  os.path.join(path_sed,filename_pick_110[j])  for j in range(0,3) ]
    x_pick110_1,y_pick110_1=ReadSED(fullfilename_pick_110[0],headerstop=38)
    x_pick110_2,y_pick110_2=ReadSED(fullfilename_pick_110[1],filetype='fits')
    x_pick110_3,y_pick110_3=ReadSED(fullfilename_pick_110[2],filetype='ascii')
    plt.figure()
    title = filename_pick_110[0]
    plt.title(title)
    plt.plot(x_pick110_1,y_pick110_1)
    plt.plot(x_pick110_2,y_pick110_2)
    plt.plot(x_pick110_3,y_pick110_3)
    plt.xlabel("$\lambda$ (nm) ")
    plt.ylabel("flux")
    