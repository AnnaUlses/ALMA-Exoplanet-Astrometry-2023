# ALMA Band 7 Astrometry

This repository contains code and data for an analysis of the feasbility of using ALMA for astrometric exoplanet searches. If you found this repository you presumably found my paper first so I won't explain the results here. The following is an explanation of what each of the data files contain. Cheers

ast_pars_149 = minimum planet mass & semi-major axis, astrometric accuracy/signature,
brighter3mag_photometry = fitted photometry parameters for larger sample of 223 (before cuts were made on declinations),
converted_mags = all of Hipparcos DR2 magnitudes converted to Gaia passband,
hip_data_149 = Hipparcos DR2 parameters for final sample of 149 stars,
observable_sample = fitted photometry parameters for observable sample,
resolved_stars = stars in this sample which are resolved by ALMA in band 7,
sample_names = just a list of the Hipparcos ID's for this sample,
sensitivity_alma = calculated sensitivies for all 8 ALMA bands for 15 and 60 minutes integration times,
snr_alma = snr in each ALMA band for each star in sample,
stellar_distances = distance to stars in sample, can be calculated from Hipparcos DR2 parallax values,
unresolved_stars = stars which remain unresolved
