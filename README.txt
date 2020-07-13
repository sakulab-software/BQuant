BQuant v0.1

---
0.  Data format
- 1つの輝度系列を縦に記入した csv ファイルを用意してください。
- 複数の系列を解析する場合は、縦の系列を横に並べてください。
  Example:

 (  MT band index ---->)
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

- 最も長い系列に長さを合わせるため、短い系列で数値がないところは "NaN" を記入してください。

---
1.  A1_estimate_MT_location.m
(1) Smooth intensity data by gaussian filter to reduce small noise
   - load data file from the 'data' folder

(2) Find MT locations. Definitions in this script:
   - positive peaks
   - Inflection points
(3) Save images and feature values of the processed data in the 'data_working' folder.
   Variables:
   -  dx      : scale  um/pixel
   -  x        : real scale coordinate of fluorescent intensity
   -  org_I    : original fluorescent intencsity
   -  smooth_I : smoothed fluorescent intencsity
   -  g_idx    : estimated locations where MT exists (not real scale)

---
2. A2_fitting_gaussians.m
(1) load the processed data from the 'data_working' folder.
(2) Fit Gaussians to intensity series by optimization algorithm without using any library for calculation.
(3)  Save images and feature values of the processed data in the 'data_working' folder.
   Variables:
   -  E             : minimum value of - Log-posterior
   -  g_height : List of Gaussian height
   -  g_center : List of Gaussian center
   -  g_stdev  : List of Gaussian standard deviation
   -  base       : baseline

NOTE 1:  Estimated results may differ slightly from run to run due to the use
of random numbers
NOTE 2:  A2_fitting_gaussians_plot.m  shows the optimization process in animation.

---
3. A3_statistics_MT_***.m
3.1  A3_1_statistics_MT_height_width_distrib.m
(1)  load the processed data from the 'data_working' folder.
(2)  Get statistics on height and width of Gaussians and save the data in the 'data_out' folder.
   -  Draw histogram and  Gamma distribution fitted to the data.
   -  Save the values in csv format

3.2  A3_2_statistics_MT_num_area_density.m
(1)  load the processed data from the 'data_working' folder.
(2)  Count Gaussian number and compute area of Gaussians and save the data in the 'data_out' folder.
   -  Draw scatter plot of number and area, and compute statistical test
   -  Save values in csv format



