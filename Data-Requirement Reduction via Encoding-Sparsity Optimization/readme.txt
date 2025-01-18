Please run Run_main_function_surface.m to plot the surface.

alldata.mat: raw data of speckles with different encoding spasities and quantities.
Recon3.6.mat: reconstructions from alldata.mat.

------------------------------------------------------
deconvtv.m                             Deconvlution
deconvtvl2.m                          Deconvlution
mergeimages.m                        Mergeimages
nnnmf.m                                 NMF
nmfdeconv.m                         UseParallel NMF and deconvlution


------------------------------------------------------
mergeimages.m
nmfdeconv.m  
Reference:
Zhu, L. et al. Repository for ’Large field-of-view non-invasive imaging through 
scattering layers using fluctuating random illumination’ (v1.0.0). Zenodo.
https://doi.org/10.5281/zenodo.5850465 (2021).

------------------------------------------------------
nnnmf.m 
Copyright 2007-2016 The MathWorks, Inc.
Reference:
M.W. Berry et al. (2007), "Algorithms and Applications for Approximate
Nonnegative Matrix Factorization," Computational Statistics and Data
Analysis, vol. 52, no. 1, pp. 155-173.


------------------------------------------------------
deconvtv.m
deconvtvl2.m
Reference:
Stanley Chan (2025). deconvtv - fast algorithm for total variation deconvolution 
(https://www.mathworks.com/matlabcentral/fileexchange/43600-deconvtv-fast-algorithm-for-total-variation-deconvolution), 
MATLAB Central File Exchange. Retrieved January 18, 2025.
S.H. Chan, R. Khoshabeh, K.B. Gibson, P.E. Gill, and T.Q. Nguyen, "An augmented Lagrangian method for total variation 
video restoration", IEEE Trans. Image Process., vol. 20, no. 11, p. 3097-3111, 2011.
