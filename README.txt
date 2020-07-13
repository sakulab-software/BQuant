BQuant v0.1

Wong JH et al. Basic proline-rich protein-mediated microtubules are essential for lobe growth and flattened cell geometry. Plant Physiol, 2019,
https://doi.org/10.1104/pp.19.00811

-----------------------------------------------------
# Data format
- Prepare a csv file in which the fluorescence intensity series is entered vertically.
- If you want to analyze more than one series, arrange the vertical series horizontally.

Example:

 (  Linescaned intensity index ---->)
1141.1, 219, 375.4
949.69, 263.86, 369.16
879.32, 171.33, 512.53
629.82, 118.64, 877.86
506.09, 236.84, 1417.47
305.48, 473.64, 1944.38
240.93, NaN, 1739.57
316.47, NaN, 1541.18
298.35, NaN,  NaN
289.82, NaN,  NaN

- To match the longest series, enter "NaN" where there is no number in the short series.
- Place data files in the "data" folder.

-----------------------------------------------------
# Analysis procedure

- Basically you just edit and run the MATLAB script "main.m".
- Edit "file", "dx", and "pp". "smooth_level" is optional.
- The analysis result of each step is saved in the "data_working" folder.

0.  step0_backgrund_subtraction.m
- Remove gradient spots in microscopic images.
- Save images and feature values.


1.  step1_estimate_MT_location.m
- Smooth intensity data by gaussian filter to reduce small noise
- Estimate MT locations as:
  *positive peaks
  *Inflection points
- Save images and feature values.


2.  step2_gaussian_fitting.m
- Fit Gaussians to intensity series by optimization algorithm without using any library for calculation.
- Determine the parameter that maximizes the posterior distribution by locally minimizing the likelihood with the Gaussian location as the prior distribution.
- Save images and feature values.

NOTE: Calculation time will be required depending on data size.


3. step3_quantification.m
- Compute statistics on height, width, and area of Gaussians.
-  Draw histogram and  Gamma distribution fitted to the data.
-  Save the values in csv format in the 'data_stat' folder.
