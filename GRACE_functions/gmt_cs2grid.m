function [grid_filter]=gmt_cs2grid(cs,radius_filter,type,destrip_method)

% Transfer spherical harmonic coefficients to grid
% 
% INPUT:
%   cs                C_lm & S_lm in CS format (spherical harmonic coefficients)
%   radius_filter     Radius of Gaussian smoothing, unit: km
%   type              grid's interval, options: 0.25 degree or 1 degree, equi-angular grid N*2N
%   destrip_method    destriping methods, options: NONE, SWENSON, CHAMBERS2007,CHAMBERS2012, CHENP3M6, CHENP4M6, DUAN
% 
% OUTPUT:
%   grid_filter       equi-angular grid N*2N, N=180 or 720
% 
% FENG Wei 11/06/2016
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

if nargin < 2, radius_filter = 0; end % default Gaussian filter radius is zero
if nargin < 3, type = 1; end % default 1 degree interval

if nargin < 4, destrip_method='NONE'; end% default destriping method is none

if ndims(cs)==3  % cs is the series of CS matrixes
    [size_a, ~, ~]=size(cs);
    if type==1
        grid_filter        = zeros(180,360,size_a);
    elseif type==0.25
        grid_filter        = zeros(721,1440,size_a);
    else
        error('Wrong input type: 1 or 0.25');
    end
    % Do destriping
    for ii=1:size_a
        cs_tmp(:,:) = cs(ii,:,:);
        cs_destrip(ii,:,:) = gmt_destriping(cs_tmp,destrip_method);
    end
    
    for ii=1:size_a
        cs_tmp(:,:)         = cs_destrip(ii,:,:);
        % Do Gaussian smoothing
        cs_fltr=gmt_gaussian_filter(cs_tmp,radius_filter);
        % SH to Grid
        if type==1
            c11cmn=[0.5 89.5 359.5 -89.5];
            lmcosi_fltr=gmt_cs2lmcosi(cs_fltr); % change the CS format
            grid_tmp=plm2xyz(lmcosi_fltr,type,c11cmn);
        elseif type==0.25
            c11cmn=[0 90 359.75 -90];
            lmcosi_fltr=gmt_cs2lmcosi(cs_fltr); % change the CS format
            grid_tmp=plm2xyz(lmcosi_fltr,type,c11cmn);
        end
        grid_filter(:,:,ii) = grid_tmp(:,:);
    end
end

if ndims(cs)==2 % cs is the one CS matrix
    if type==1
        grid_filter        = zeros(180,360);
    elseif type==0.25
        grid_filter        = zeros(721,1440);
    else
        error('Wrong input type: 1 or 0.25');
    end
    % Do destriping
    cs_destrip = gmt_destriping(cs,destrip_method);
    % Do Gaussian smoothing
    cs_fltr=gmt_gaussian_filter(cs_destrip,radius_filter);
    % SH to Grid
    if type==1
        c11cmn=[0.5 89.5 359.5 -89.5];
        lmcosi_fltr=gmt_cs2lmcosi(cs_fltr); % change the CS format
        grid_filter=plm2xyz(lmcosi_fltr,type,c11cmn);
    elseif type==0.25
        c11cmn=[0 90 359.75 -90];
        lmcosi_fltr=gmt_cs2lmcosi(cs_fltr); % change the CS format
        grid_filter=plm2xyz(lmcosi_fltr,type,c11cmn);
        
    end
end

end


