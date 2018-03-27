# README.md:

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


### Read cathalog of SNLS and plot the sky

- **ReadSNLSCatalog.ipynb** : Read SNLS catalog and extract relevant magnitudes into a csv file *"SNLS_catD1D2D4.csv"*
- **PlotSNLSSky.ipynb** : Skymaps SNLS catalog in the sky from file *"SNLS_catD1D2D4.csv"*
- **SNLSColors.ipynb** : View the SNLS catalog colors and their correlations from file *"SNLS_catD1D2D4.csv"*




### About Pysynphot models:

- **CheckGridSED.ipynb** : check the SED catalog produced by **libpysynphotgridsed.py**

### Calculations of magnitudes from SED
- **ShowPlot\_libcolors.ipynb** : Simple example to show how one can compute SNLS magnitudes and colors from SED



## libraries
- **libCFHTFilters.py** : Show all transmissions of CFHT

- **libSNLSPhotometry.py** : Calculate the SNLS magnitudes from SED

- **libpysynphotgridsed.py** : Generate the SED from a model.
Only phoenix stars are implemented for the moment. A fits file is produced.

