# GRACE_Matlab_Toolbox
GRACE_Matlab_Toolbox
GRACE Matlab Toolbox (GRAMAT) contains a set of open-source functions to process GRACE level-2 spherical harmonic coefficient products. Data processing functions in the GMT contain: (1) destriping of SH coefficients to remove “north-to-south” stripes and Gaussian smoothing, (2) spherical harmonic analysis and synthesis, (3) analyzing and reducing leakage effect in GRACE-derived mass variations, (4) analyzing regional mass variations in spatial and temporal domains. Matlab GUIs are presented to facilitate the use of functions in the GRAMAT.

The work has been published in Earth Science Informatics. Please cite this paper if you use this toolbox https://link.springer.com/article/10.1007/s12145-018-0368-0

Our recent papers on groundwater storage variations using this toolbox are:

(1)	Feng, W., C. Shum, M. Zhong, Y. Pan, Groundwater storage changes in China from satellite gravity: An overview, Remote Sensing, 2018, 10(5), 674.

(2)	Feng W., M. Zhong, J.-M. Lemoine, R. Biancale, H.-T. Hsu, J. Xia, Evaluation of groundwater depletion in North China using the Gravity Recovery and Climate Experiment (GRACE) data and ground-based measurements, Water Resources Research, 2013, doi:10.1002/wrcr.20192

We also appreciate it if you refer to these papers in your publications.


-- Change History --

02/25/2022, Yu Zhang, use centering and scaling in curve fitting in gmt_destriping_chen,
            otherwise Warning: Polynomial is badly conditioned for P4M6

02/21/2022, Yu Zhang, fix a bug in gmt_destriping for choosing CHENP4M6 method