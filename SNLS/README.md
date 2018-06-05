# README.md:
# ===================

- author : Sylvie Dagoret-Campagne
- date :  March 28th 2018 
- update : Jun 5th 2018




To produce a catalog of SED from Pysynphot according distribution of stars
in SNLS field.


## SNLS:

### information from SNLS retrieved from :

- http://supernovae.in2p3.fr/snls_sdss/

comment from LPNHE :

We release 4 products:

Catalogs of natural AB magnitudes for stars selected as tertiary standards.
MegaCam instrument transmission functions at 1.25 airmass.
    Covariance matrices for the estimated calibration error.


### CFHT & MEGACAM :
http://www.cfht.hawaii.edu/fr/projets/


## Notebooks

### Read CFHT/Megacam transmissions
- **MegaCAMTransm.ipynb** : Extract the average transmission for **Megacam** and save it in files
*"all_SNLS_transm.csv"* and *"tot_SNLS_transm.csv"*


### Read catalog of SNLS and plot the sky

- **ReadSNLSCatalog.ipynb** : Read SNLS catalog and extract relevant magnitudes into a csv file *"SNLS_catD1D2D4.csv"*
- **PlotSNLSSky.ipynb** : Skymaps SNLS catalog in the sky from file *"SNLS_catD1D2D4.csv"*
- **SNLSColors.ipynb** : View the SNLS catalog colors and their correlations from file *"SNLS_catD1D2D4.csv"*



### About Pysynphot models:

- **CheckGridSED.ipynb** : check the SED catalog produced by **libpysynphotgridsed.py**: this is valid for pickles and phoenix stars

### Calculations of SNLS magnitudes from SED

Two csv file are produced either with magnitudes or colors.

- **ShowPlot\_libcolors.ipynb** : Simple example to show how one can compute SNLS magnitudes and colors from SED (deprecated)

- **ShowColorPlot\_libcolors2.ipynb** : This one calculate the  AB magnitudes with Phoenix model : it is the one to be used

- **ShowColorPlot\_libcolors2\_pickles.ipynb** : This one calculate the  AB magnitudes for pickles : it is the one to be used

- **ShowColorPlot_libcolors2_regenerated.ipynb** : Do the same thing on regenerated SED

### Compare the model and the SNLS magnitudes


The goal is to compare colors of SNLS catalog wrt SNLS colors calculated on SED.
This is a validation of the SED catalog used and the way to calculate the colors or magnitudes on SED catalog.

- **CompareSNLSandModel.ipynb** : compare the magnitudes 

- **CompareSNLSandModel\_phoenix.ipynb** : compare the magnitudes for phoenix stars

- **CompareSNLSandModel_pickles.ipynb**	 : compare the magnitudes for pickles stars

From this result, one can validate pickles catalog.


-  **CompareSNLSandModel_regenerated.ipynb** : compare the colors (AB) 
-  **add CompareSNLSandModel_regenerated.ipynb** : Cmpare the magnitudes (AB)

### Select the Nearest Neighbourg SED

The goal is to associate an SED to each obj of SNLS magnitude catalog.
This association is done in color space, where the magnitudes are the AB system.
As the color comparison has been validated previously with pickles, only pickles star match is performed.

- **SelectKNearestNeighbors.ipynb** : it is implemented for pickles star matching using K Nearest Neighbor implemented in scikit learn.


### Regenerate the SED

- **Regenerate_SED.ipynb** : based on picles

### Other Tools
- **CalibrationSpectra.ipynb** : Check how one can use calibrated magnitude system in pynsynphot : it is a pedagogical system.

- **MEGACAM\_FromPySynPhot.ipynb** : Check what physynphot gives as MEGACAM transmission. These transmissions are bad. Better use the one calculated here.

- **CheckGridSED.ipynb** : Simple check of SED produced by libpysynphot	

- **MergeSEDandMAG.ipynb**. : Study color plots for phoenix, including the impact of star temperature, metallicity and gravity.

- **ShowColorPlot\_libcolors_abmag.ipynb** : Tool to play with a signle SED to compute magnitudes and colors in AB magnitudes.

- **Extract\_FewSEDSampl.ipynb** : Deprecated

- **SelectSED_pickle.ipynb** and **SelectSED_phoenix.ipynb**   select the SED produced by **libpysynphotgridsed.py**

## libraries
- **libCFHTFilters.py** : Show all transmissions of CFHT

- **libSNLSPhotometry.py** : Calculate the SNLS magnitudes from SED

- **libpysynphotgridsed.py** : Generate the SED from a model.
Only phoenix and pickles stars are implemented for the moment. A fits file is produced.

