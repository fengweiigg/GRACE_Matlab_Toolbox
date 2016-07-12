function [cs_leakagefree]=gmt_cs2leakagefreecs(cs,type,destrip_method,radius_filter)

% Remove the leakage from ocean(or land) in spectrul domain
% 
% INPUT:
%   cs                C_lm & S_lm in CS format (spherical harmonic coefficients)
%   radius_filter     Radius of Gaussian smoothing, unit: km
%   type              'ocean', Ocean Leakage Reduction: remove land's leakage effect 
%                     'land', Land Leakage Reduction: remove ocean's leakage effect
%   destrip_method    destriping methods: options: NONE, SWENSON, CHAMBERS2007,CHAMBERS2012, CHENP3M6, CHENP4M6, DUAN
% 
% OUTPUT:
%   cs_leakagefree    Leakage-removed spherical harmonic coefficients
% 
% FENG Wei 05/09/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

if nargin < 3
    destrip_method='NONE';
    radius_filter=0;
end


if ndims(cs)==2 % cs is the one CS matrix
    [rows,cols] = size(cs);
    if rows ~= cols					% field is not in CS-format
        error('Check format of gravity field data (CS-format).')
    end
    degree_max=rows-1; % maximum degree
elseif ndims(cs)==3 % cs is the series of CS matrixes
    [~,rows,cols] = size(cs);
    if rows ~= cols					% field is not in CS-format
        error('Check format of gravity field data (CS-format).')
    end
    degree_max=rows-1; % maximum degree
end
    
if strcmp(type,'ocean')
    load lg_msk_land_025.mat
elseif strcmp(type,'land')
    load lg_msk_ocean_025.mat
else
    error('Wrong type option, (Options: land, ocean)');
end

% STEP1: Transfer original SH coefficients to original grids
grid_original=gmt_cs2grid(cs,0,0.25,'NONE');

% STEP2: Mask out Land/Ocean for Land/Ocean reduction
if ndims(grid_original)==3  % grid is 3-D matrix, grid(lat, lon, time)
    for ii=1:size(grid_original,3)
        grid_masked(:,:,ii)=grid_original(:,:,ii).*lg_msk;
    end
else % grid is 2-D matrix, grid(lat, lon)
    grid_masked=grid_original.*lg_msk;
end

% STEP3: Transfer grids to SH coefficients again (Land/Ocean masked)
cs_masked=gmt_grid2cs(grid_masked,degree_max);

% STEP4: Do destriping and filter to SH coefficients (Land/Ocean masked)
if ndims(cs_masked)==3  % 3-D matrix
    [num_file, ~, ~]=size(cs_masked);
    % Do destriping
    for ii=1:num_file
        cs_tmp(:,:) = cs_masked(ii,:,:);
        cs_destrip(ii,:,:) = gmt_destriping(cs_tmp,destrip_method);
    end
    % Do Gaussian smoothing
    for ii=1:num_file
        cs_tmp(:,:)         = cs_destrip(ii,:,:);
        cs_leakage=gmt_gaussian_filter(cs_tmp,radius_filter);
    end
    % Do destriping
    for ii=1:num_file
        cs_tmp(:,:) = cs(ii,:,:);
        cs_destrip(ii,:,:) = gmt_destriping(cs_tmp,destrip_method);
    end
    % Do Gaussian smoothing
    for ii=1:num_file
        cs_tmp(:,:)         = cs_destrip(ii,:,:);
        cs_before_leakage=gmt_gaussian_filter(cs_tmp,radius_filter);
    end    
    
end
if ndims(cs_masked)==2  % 2-D square matrix
    % Do destriping
    cs_destrip = gmt_destriping(cs_masked,destrip_method);
    % Do Gaussian smoothing
    cs_leakage=gmt_gaussian_filter(cs_destrip,radius_filter);
    % Do destriping
    cs_destrip = gmt_destriping(cs,destrip_method);
    % Do Gaussian smoothing
    cs_before_leakage=gmt_gaussian_filter(cs_destrip,radius_filter);
end

% STEP5: Remove leakage
cs_leakagefree=cs_before_leakage-cs_leakage;





