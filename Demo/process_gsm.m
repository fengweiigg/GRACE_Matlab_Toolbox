clear all
clc
%% process GSM data from CSR
% Please modify the control file before running
GRACE_Matlab_Toolbox_preprocessing_core('GRAMAT_Control_File_csr_swenson.txt');
GRACE_Matlab_Toolbox_preprocessing_core('GRAMAT_Control_File_csr_none.txt');

%% process GSM data from GFZ
GRACE_Matlab_Toolbox_preprocessing_core('GRAMAT_Control_File_gfz_swenson.txt');
GRACE_Matlab_Toolbox_preprocessing_core('GRAMAT_Control_File_gfz_none.txt');

%% process GSM data from JPL
GRACE_Matlab_Toolbox_preprocessing_core('GRAMAT_Control_File_jpl_swenson.txt');
GRACE_Matlab_Toolbox_preprocessing_core('GRAMAT_Control_File_jpl_none.txt');

%% process GSM data from GRGS
GRACE_Matlab_Toolbox_preprocessing_core('GRAMAT_Control_File_grgs.txt');


%% TEST
load /Users/fengwei/matlab/GRACE_Matlab_Toolbox/GRACE_results/cs_gsm_csr_swenson_2002_2014_60degree.mat

% Calculate mass variation gridded data for 2004/02
grid_1=gmt_cs2grid(cs_grace(20,:,:),300,1);
figure(1)
gmt_grid2map(grid_1*100,30,'cm','Global mass variations in 2004/02 from GRACE')

% Calculate trend map 
grid_2=gmt_cs2grid(cs_grace,300,1); % 
[ ~, ~, ~,~, ~, ~, ~, ~, Trend, ~, ~, ~, ~] = gmt_harmonic(time,[],grid_2);
figure(2)
gmt_grid2map(Trend*100,5,'cm/yr','Trend map of mass variations for 2002-2014 from GRACE')

% Compute the time series of mass variations in the North China Plain
% keep in mind, no leakage reduction and rescaling should be implemented to get unbiased time series
grid_2=gmt_cs2grid(cs_grace,300,1); % 
plot_region_ncp=gmt_grid2series(grid_2,'NCP.bln','line');
figure(3)
plot(time,plot_region_ncp,'-s')
