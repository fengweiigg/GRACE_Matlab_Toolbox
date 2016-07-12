function [cs]=gmt_grid2cs(grid,degree_max)

% Transfer gridded field to spherical harmonic coefficients
% Spherical harmonic analysis
% 
% INPUT:
%   grid          Regular global gridded field (equi-angular grid N*2N)
%   degree_max    Maximum degree of the expansion
%  
% OUTPUT:
%   cs            C_lm & S_lm in CS format (spherical harmonic coefficients)
%
% xyz2plm function is created by Frederik Simons from http://geoweb.princeton.edu/people/simons/software.html
% lmcosi        Matrix listing l,m,cosine and sine coefficients
% for full edition of spherical harmonic analysis, please use xyz2plm function
%
% FENG Wei 11/06/2016
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


grid(logical(isnan(grid)))=0; % set NANs in grid to ZEROs


if ndims(grid)==2
    [lmcosi,dw]=xyz2plm(grid,degree_max);%do spherical harmonic analysis
    cs=gmt_lmcosi2cs(lmcosi); % change the CS format
end

if ndims(grid)==3
    for k=1:size(grid,3)
        grid_tmp(:,:)=grid(:,:,k);
        [lmcosi,dw]=xyz2plm(grid_tmp,degree_max); %do spherical harmonic analysis
        cs_tmp=gmt_lmcosi2cs(lmcosi); % change the CS format
        cs(k,:,:)=cs_tmp;
    end
end

end