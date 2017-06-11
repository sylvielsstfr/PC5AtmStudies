# LSSTFiltersKG
===============

- Author :Sylvie Dagoret-Campagne
- affilication : LAL/IN2P3/CNRS
- date : May 27th 2017

Official LSST transmission including
- Filters
- Throughput
- CCD QE
Courteousy Kirk Gilmore, SLAC 


# Standalone python script
-------------------------
ReadIdealFilters.py

do
python ReadIdealFilters.py
or
ipython ReadIdealFilters.py

It should show now the plot on your screen.

notice the matplotlib config file "matplotlibrc" is inside  LSSTFiltersKG
directory


# ipython notebooks
-------------------

do :
ipython notebook

and you shoud see the notebook list

## ReadIdealFilters.ipynb
-------------------------

Example of notebook to read only filters (ascii file)


## ReadThroughputCCD.ipynb
--------------------------
Example to read throughput and CCD (excel file)

## ReadFilterThroughputAndCCDQE.ipynb
-------------------------------------
Example to read  both filters and CCD and Throughput
