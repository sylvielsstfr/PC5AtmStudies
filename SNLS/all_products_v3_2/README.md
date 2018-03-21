           SNLS and SDSS-SN joined Photometric Calibration
           ===============================================

Author: Marc Betoule
Date: 2012-12-13 17:43:23 CET


Table of Contents
=================
1 Overview
2 Description of calibration products
    2.1 Calibrated catalogs of tertiary standard stars in the science fields
    2.2 MegaCam filters
    2.3 Covariance matrix
    2.4 MegaCam photometric flatfields
3 Updating the calibration


1 Overview 
-----------

This work is part of the SDSS+SNLS joint analysis. It provides
photometric calibration standards for the SNLS and SDSS surpernova
surveys. The paper associated with products presented here can be
found here (see in particular Appendix E).

We release 4 products:
  - Catalogs of *natural AB magnitudes* for stars selected as
    tertiary standards.
  - MegaCam instrument transmission functions at 1.25 airmass.
  - Covariance matrices for the
    estimated calibration error.
All products are delivered as a [single tarball].

In addition, we release the material required to update the
calibration according to new releases of the HST CALSPEC
spectrophotometry.
  - Measured magnitudes of the primary standards
  - Covariance matrices of primary standard measurements


  [single tarball]: ./all_products_v3_2.tgz

2 Description of calibration products 
--------------------------------------

2.1 Calibrated catalogs of tertiary standard stars in the science fields 
=========================================================================

Those catalogs deliver *natural AB magnitudes* for stars selected as
tertiary standards for the [SNLS] and [SDSS supernova surveys].  The
catalogs deliver MegaCam \(u_{M} g_{M} r_{M} i_{M} i2_{M}\) and
$z_{M}$ AB magnitudes. The SNLS cover the 4 deep fields of the [CFHTLS]:
 - =catalogs/D1_release_ab_v3_2.ascii=
 - =catalogs/D2_release_ab_v3_2.ascii=
 - =catalogs/D3_release_ab_v3_2.ascii=
 - =catalogs/D4_release_ab_v3_2.ascii=

The release of the catalogs is associated with a position dependant
model for MegaPrime instrument response [(see below MegaCam filters)].

The SDSS SN survey cover a $2.5^^\times120^^$ degree region
centered on the celestial equator and designated as stripe 82. The
calog deliver AB magnitudes in the 2.5m SDSS photometric system $ugriz$.
 - =catalogs/S82_release_ab_v3_2.ascii=

The associated model for the system response of the SDSS 2.5-meter is
delivered in [Doi et al. (2010)]. We use the average filter model.



[SNLS]: http://www.cfht.hawaii.edu/SNLS/
[SDSS supernova surveys]: http://www.sdss.org/supernova/aboutsupernova.html
[CFHTLS]: http://www.cfht.hawaii.edu/Science/CFHTLS/
[(see below MegaCam filters)]: sec-2-2
[Doi et al. (2010)]: http://adsabs.harvard.edu/abs/2010AJ....139.1628D

2.2 MegaCam filters 
====================
   
The model for the transmission curves incorporate the following contributions:
  - CCD quantum efficiency curves
  - Position-dependent filter transmission curves at several radii
  - Transmission curves for the wide-field corrector optics
  - Reflectivity of the primary mirror
  - Average transmission curve of the atmosphere above Mauna-Kea at an
    airmass of 1.25, as published in [Buton et al. (2012)].

Each component is provided as a tabulated function in the MegaCam_v3.2
directory. Respectively, the corresponding files are:
  - =QE_camera_high_res_model.dat=
  - ={u,g,r,i,y,z}*.list= for the \(u_{M} g_{M} r_{M} i_{M} i2_{M}\)
    and $z_{M}$ filters
  - =CFHT_MegaPrime_Transmission.dat=
  - =CFHT_Primary_Transmission.dat=
  - =SNIFS_extinction_buton2012_with_tl_X1_25.dat=

The effective filter is the product of the five components. Filter
transmission curves must be interpolated a given position.


[Buton et al. (2012)]: http://adsabs.harvard.edu/abs/2012arXiv1210.2619B

2.3 Covariance matrix 
======================

The covariance matrix of the \(u_{M} g_{M} r_{M} i_{M} i2_{M} z_{M} u
g r i z\) calibration offsets is delivered as a fits image in the file:

=covariance/ab_offset_cov_mat_full.fits=

2.4 MegaCam photometric flatfields 
===================================

We do not release the flat-field photometric corrections, as they are
specific to the photometry method used in the calibration paper.

3 Updating the calibration 
---------------------------

Calibration offsets for the tertiary catalogs have been computed by
comparing measured magnitudes of a few CALSPEC primary
spectrophotometric standards to synthetic magnitudes. Offsets computed
with the release 003 of the CALPSEC spectra have been applied to the
catalogs.

Any revision of the HST calibration can be propagated to the SNLS and
SDSS tertiary standards by recomputing calibration offsets. For this
purpose we release the measured magnitudes of the CALPSEC standards in
the tertiary catalog magnitude system, as well as the covariance
matrix of those measurements:
  - =CALSPEC_v3.2/hst_mmags.ascii=
  - =covariance/cov_meas.fits=


