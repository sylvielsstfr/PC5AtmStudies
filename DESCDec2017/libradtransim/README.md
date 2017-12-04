README.md
========

to have verbose info in a file do:

> (python simulate_transparency_CTIO_ScattAbs.py -z 1 -w 4 -o 300) >& verbose_mw.txt


# library to simulate atmospheric transparency with libradtran

- libsimulateTranspLSSTScattAbs.py :
Simulate only scattering and absorption

- libsimulateTranspLSSTScattAbsAer.py :
Simulate scattering , absorption and also aerosols


# python program to simulate atmosphere
- simulate_transparency_LSST_ScattAbs.py

# notebooks to check how aerosols are simulated in libradtran

- ScanAerosols.ipynb  : default aerosol config in libradtran				 
- ScanAerosols2.ipynb : tunable formula config in libdradtran , tau = beta*(wl0/wl)^alpha
- ScanAerosols3.ipynb : user analytic formula for aerosols defined by me


# python tool for libradtran
- UVspec.py											