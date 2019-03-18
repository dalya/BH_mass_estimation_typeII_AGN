Online data information 

########################################################################################
########################################################################################

typeI_AGN_spectra.csv

The file contains the continuum-subtracted emission line spectra of the AGN used in this study. The spectra are already ordered according to the detected sequence, such that neighboring spectra in the file are neighbors in the sequence. All spectra are interpolated to a common wavelength range between 4000 and 7000 A with 1 A resolution. Every row in the file corresponds to a single type I AGN spectrum, and contains 3000 flux values. There are 1941 rows in the file, corresponding to the 1941 type I AGN considered in this work.

########################################################################################
########################################################################################

typeI_AGN_metadata.csv

The file contains additional information on the 1941 type I AGN considered in this study. The rows are ordered according to the detected sequence, and thus they match the ordering of the emission line spectra from file ‘typeI_AGN_spectra.csv’. 

The columns in the file:

> plate - the SDSS plate of the spectrum.
> mjd - the SDSS MJD of the spectrum.
> fiber - the SDSS fiber of the spectrum. 
> FWHM_broad_Halpha - the FWHM of the broad Halpha line, in units of km/sec, taken from the Shen et al. (2011) catalog. 
> log_L_OIII_5007 - the log of the narrow [OIII]5007A luminosity, in units of erg/sec, taken from the Shen et al. (2011) catalog.
> log_L_narrow_Hbeta - the log of the narrow Hbeta luminosity, in units of erg/sec, taken from the Shen et al. (2011) catalog.
> EW_broad_Hbeta - the equivalent with of the broad Hbeta line, in units of A, taken from the Shen et al. (2011) catalog.
> EW_broad_FeII - the equivalent width of the broad FeII line, in units of A, taken from the Shen et al. (2011) catalog.
> EW_OIII_5007 - the equivalent width of the [OIII]5007A line, in units of A, taken from the Shen et al. (2011) catalog.
> r_mag - r band magnitude from the SDSS.
> W1_mag - W1 magnitude from WISE.
> W2_mag - W2 magnitude from WISE.
> W3_mag - W3 magnitude from WISE.
> W4_mag - W4 magnitude from WISE.
> log_BH_mass - log of the virial black hole mass, in units of solar masses, taken from the Shen et al. (2011) catalog.
> log_L_bol - log of the AGN bolometric luminosity, in units of erg/sec, taken from the Shen et al. (2011) catalog.
> log_Eddington_ratio - log of the Eddington ratio, defined as Lbol/LEdd, taken from the Shen et al. (2011) catalog.

########################################################################################
########################################################################################

typeI_AGN_stacked_measurements.csv

The file contains the narrow log([OIII]/Hbeta) versus FWHM(broad Halpha) measurements obtained from emission line decomposition of the stacked spectra presented in the paper.

########################################################################################
########################################################################################

typeII_AGN_metadata.csv

The file contains information on the ~10,000 type II AGN considered in this study.
The columns in the file:

> plate - the SDSS plate of the spectrum.
> mjd - the SDSS MJD of the spectrum.
> fiberid - the SDSS fiber of the spectrum.
> z - the redshift of the object.
> h_beta_flux - the flux of the narrow Hbeta emission line, in units of 10^-17 erg/sec/cm^2, taken from the MPA/JHU catalog.
> h_beta_flux_err - the uncertainty of the flux of the narrow Hbeta emission line.
> oiii_5007_flux - the flux of the narrow [OIII]5007A emission line, in units of 10^-17 erg/sec/cm^2, taken from the MPA/JHU catalog.
> oiii_5007_flux_err - the uncertainty of the flux of the narrow [OIII]5007A emission line.
> h_alpha_flux - the flux of the narrow Halpha emission line, in units of 10^-17 erg/sec/cm^2, taken from the MPA/JHU catalog.
> h_alpha_flux_err - the uncertainty of the flux of the narrow Halpha emission line.
> nii_6584_flux - the flux of the narrow [NII]6584A emission line, in units of 10^-17 erg/sec/cm^2, taken from the MPA/JHU catalog.
> nii_6584_flux_err - the uncertainty of the flux of the narrow [NII]6584A emission line.
> log_stellar_sigma - log of stellar velocity dispersion, in units of km/sec, measured using pPXF.
> log_bh_mass - log black hole mass, in units of solar mass.
> objid - SDSS DR7 object ID.
> psfMag_u - SDSS PSF u-band magnitude.
> psfMag_g - SDSS PSF g-band magnitude.
> psfMag_r - SDSS PSF r-band magnitude.
> psfMag_i - SDSS PSF i-band magnitude.
> psfMag_z - SDSS PSF z-band magnitude.
> psfMagErr_u - SDSS PSF u-band magnitude uncertainty.
> psfMagErr_g - SDSS PSF g-band magnitude uncertainty.
> psfMagErr_r - SDSS PSF r-band magnitude uncertainty.
> psfMagErr_i - SDSS PSF i-band magnitude uncertainty.
> psfMagErr_z - SDSS PSF z-band magnitude uncertainty. 
> mendel_logM_p50 - median log total stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logM_p16 - 16th percentile of log total stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logM_p84 - 84th percentile of log total stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logMt_p50 - median Log total stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logMt_p16 - 16th percentile of Log total stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logMt_p84 - 84th percentile of Log total stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logMb_p50 - median log bulge stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logMb_p16 - 16th percentile of log bulge stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logMb_p84 - 84th percentile of log bulge stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logMd_p50 - median log disk stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logMd_p16 - 16th percentile of log disk stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> mendel_logMd_p84 - 84th percentile of log disk stellar mass, in units of solar mass, taken from Mendel et al. (2014).
> simard_b_t_g - g-band bulge fraction, taken from Simard et al. (2011).
> simard_e_b_t_g - uncertainty of g-band bulge fraction.
> simard_b_t_r - r-band bulge fraction, taken from Simard et al. (2011).
> simard_e_b_t_r - uncertainty of r-band bulge fraction.
> simard_Rhlg - g-band galaxy semi-major axis, half-light radius, in units of kpc, taken from Simard et al. (2011).
> simard_Rhlr - r-band galaxy semi-major axis, half-light radius, in units of kpc, taken from Simard et al. (2011).
> simard_Rchl_g - g-band galaxy circular half-light radius, in units of kpc, taken from Simard et al. (2011).
> simard_Rchl_r - r-band galaxy circular half-light radius, in units of kpc, taken from Simard et al. (2011). 
> simard_Re - Bulge semi-major effective radius, in units of kpc, taken from Simard et al. (2011).
> simard_e_Re - uncertainty of  bulge semi-major effective radius.
> simard_e - Bulge ellipticity, taken from Simard et al. (2011).
> simard_e_e - uncertainty of bulge ellipticity.
> simard_nb - Bulge Sersic index, taken from Simard et al. (2011).
> simard_e_nb - uncertainty of bulge Sersic index.
> simard_PpS - F-test probability that a B+D model is not required compared to a pure Sersic model, taken from Simard et al. (2011).
> simard_Pn4 - F-test probability that a free nb B+D model is not required compared to a fixed nb=4 B+D model (where nb is the Sersic index), taken from Simard et al. (2011).

