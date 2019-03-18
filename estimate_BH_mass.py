import numpy
from scipy.interpolate import interp1d
from astropy import units as u
from astropy import constants as c
from astropy import cosmology as cos
cosmo = cos.FlatLambdaCDM(H0=70, Om0=0.3)

def return_fwhm_given_narrow_lines(oiii_hbeta):
    """
    function returns the FWHM(broad Halpha) given the measured narrow [OIII]/Hbeta line ratio. 
    function does not filter out objects which are beyond the range, since it was already done when I constructed the 
    files.
    """
    log_oiii_hbeta = numpy.log10(oiii_hbeta)
    log_fwhm = 1.406 * log_oiii_hbeta + 2.580 #1.23 * log_oiii_hbeta + 2.716    
    fwhm = 10**log_fwhm
    return fwhm

def deredden_spectrum_Fitzpatric99(wl, spec, E_bv):
    """
    function dereddens a spectrum based on the given E(B-V) value and Fitzpatric99 model
    """
    # dust model
    wls = numpy.array([ 2600,  2700,  4110,  4670,  5470,  6000, 12200, 26500])
    a_l = numpy.array([ 6.591,  6.265,  4.315,  3.806,  3.055,  2.688,  0.829,  0.265])
    f_interp = interp1d(wls, a_l, kind="cubic")

    a_l_all = f_interp(wl)
    A_lambda = E_bv * a_l_all
    spec_dered = spec * 10 ** (A_lambda / 2.5)

    return spec_dered


def return_ccm_a_b_model(y):
    """
    function returns the a,b coefficients for the CCM model
    """
    a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 -0.77530*y**6 + 0.32999*y**7
    b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 -0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7

    return a, b

def deredden_spectrum_CCM(wl, spec, E_bv):
    """
    function dereddens the observed spectrum according to the CCM 1989 extinction curve, which is more appropriate!
    """
    lambda_wl = wl / 10000.0
    lambda_x = 1.0 / lambda_wl
    lambda_y = lambda_x - 1.82

    a_arr, b_arr = return_ccm_a_b_model(lambda_y)

    Rv = 3.1 # this is appropriate for MW-type extinction
    Av = E_bv * 3.1

    # compute A(lambda) at each point in the array
    A_lambda = Av * (a_arr + b_arr/Rv)

    spec_dered = spec * 10**(A_lambda/2.5)

    return spec_dered

def return_hbeta_luminosity(hbeta_flux, halpha_flux, z):
    """
    function returns Hbeta luminosity after applying dust correction
    """
    ebv = 2.33 * numpy.log10((halpha_flux/hbeta_flux)/2.85)
    hbeta_flux_dered = deredden_spectrum_CCM(4861., hbeta_flux, ebv)
    dl = cosmo.luminosity_distance(z).cgs
    hbeta_lum = 4 * numpy.pi * (dl.value)**2 * 10**-17 * hbeta_flux_dered
    return hbeta_lum

def return_log_bolometric_luminosity(hbeta_flux, halpha_flux, oiii_flux, z):
    """
    function returns the AGN bolometric luminosity according to Hagai's formula
    """
    hbeta_lum = return_hbeta_luminosity(hbeta_flux, halpha_flux, z)
    log_hbeta_lum = numpy.log10(hbeta_lum)
    
    log_oiii_hbeta = numpy.log10(oiii_flux/hbeta_flux)
    val = 0.31 * (log_oiii_hbeta - 0.6)
    val[val < 0] = 0
    log_lbol = log_hbeta_lum + 3.48 + val
        
    return log_lbol

def return_l5100_from_lbol(hbeta_flux, halpha_flux, oiii_flux, z):
    """
    function returns the L(5100) luminosity, given the bolometric correction
    """
    log_lbol = return_log_bolometric_luminosity(hbeta_flux, halpha_flux, oiii_flux, z)
    
    log_l5100 = 1.09 * log_lbol - 5.23
    l5100 = 10**log_l5100
    return l5100

def return_bh_mass(hbeta_flux, halpha_flux, oiii_flux, z):
    """
    function returns the BH mass of the type II AGN
    the input is assumed to be numpy arrays of the sample of type II AGN
    hbeta_flux: the flux of the narrow hbeta line, in units of 10^-17 erg/sec/cm^2
    halpha_flux: the flux of the narrow halpha line, in units of 10^-17 erg/sec/cm^2
    oiii_flux: the flux of the narrow [OIII] line, in units of 10^-17 erg/sec/cm^2
    z: the redshift of eahc object
    """
    # estimate FWHM(broad Halpha)
    oiii_hbeta = oiii_flux / hbeta_flux
    log_oiii_hbeta = numpy.log10(oiii_hbeta)
    fwhm_broad_halpha = numpy.array(return_fwhm_given_narrow_lines(oiii_hbeta)).astype(numpy.float64)
    # estimate AGN continuum luminosity
    l5100 = numpy.array(return_l5100_from_lbol(hbeta_flux, halpha_flux, oiii_flux, z)).astype(numpy.float64)
    log_l5100 = numpy.log10(l5100)
    # estimate BH mass
    log_mbh = numpy.log10(1.075) + 6.90 + 0.54*(log_l5100 - 44.) + 2.06*numpy.log10(fwhm_broad_halpha/(10**3))
    # remove objects for which the [OIII]/Hbeta is outside the region
    log_mbh[(log_oiii_hbeta < 0.55) | (log_oiii_hbeta > 1.1)] = numpy.nan
    return log_mbh