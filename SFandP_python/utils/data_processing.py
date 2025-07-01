from astropy.io import fits
import numpy as np

def load_and_process_fits(filepath):
    hdul = fits.open(filepath)
    data = hdul[1].data
    cataid = data.field('CATAID')
    mass = data.field('mass!')
    mass_err = data.field('mass_err!')  # (84th percentile + 16th percentile)/2
    sfr_halpha = data.field('SFR')

    log_mass = np.log10(mass)
    log_mass_err = mass_err / (mass * np.log(10))
    log_sfr_halpha = np.log10(sfr_halpha)
    
    return cataid, log_mass, log_mass_err, log_sfr_halpha
    