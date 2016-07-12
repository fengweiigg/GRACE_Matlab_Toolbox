function cs_fltr=gmt_gaussian_filter(field,radius_filter)

% Do Gaussian smoothing to spherical harmonic coefficients 
% 
% INPUT:
%   field           C_lm & S_lm in CS format (spherical harmonic coefficients, |C\S|, (L+1)x(L+1) matrix)
%   radius_filter   Radius of Gaussian smoothing, unit: km
%
% OUTPUT:
%   cs_fltr         C_lm & S_lm in SC format (spherical harmonic coefficients, /S|C\, (L+1)x(2L+1))
%
% FENG Wei 18/12/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

if ndims(field)==2
    [rows,cols] = size(field);
    if rows ~= cols					% field is not in CS-format
        error('Check format of gravity field data (CS-format).')
    end
    lmax=rows-1; % maximum degree
    W = gmt_gaussian(lmax,radius_filter);
    sc(:,:) = gmt_cs2sc(field);
    % gauss filter
    for ll = 0:1:lmax
        scnew(ll+1,:) = W(ll+1)*sc(ll+1,:);
    end
    % chang to CS format back
    cs_fltr = gmt_sc2cs(scnew);
    
elseif ndims(field)==3
    [num_file,rows,cols] = size(field);
    if rows ~= cols					% field is not in CS-format
        error('Check format of gravity field data (CS-format).')
    end
    lmax=rows-1; % maximum degree
    W = gmt_gaussian(lmax,radius_filter);
    for i=1:num_file
        % change to SC format
        cs_tmp(:,:)   = field(i,:,:);
        sc(:,:) = gmt_cs2sc(cs_tmp);
        % gauss filter
        for ll = 0:1:lmax
            scnew(ll+1,:) = W(ll+1)*sc(ll+1,:);
        end
        % chang to CS format back
        cs_fltr(i,:,:) = gmt_sc2cs(scnew);
    end
else
    error('Please check the dimension of input CS data in ''gmt_gaussian_filter'' funciton');
end



