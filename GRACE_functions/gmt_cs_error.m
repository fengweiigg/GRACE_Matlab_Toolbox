function [grid_error]=gmt_cs_error(time,cs,type)

% Estimate error in CS using Wahr et al. 2006 method
% Reference: Wahr J., Swenson S. & Velicogna I. Accuracy of GRACE mass estimates. Geophys. Res. Lett. 33, L06401, doi:06410.01029/02005GL025305 (2006).
% spherical harmonic coefficients to grid
% 
% INPUT:
%   cs                C_lm & S_lm in CS format (spherical harmonic coefficients)
%   type              grid's interval, options: 0.25, 0.5 or 1 degree
%
% OUTPUT:
%   grid_error       
% 
% FENG Wei 07/02/2017
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn
%
[size_a, ~, ~]=size(cs);
for ii=1:size_a
    cs_tmp(:,:,ii) = cs(ii,:,:);
end
    

[ ~, ~, ~,~, ~, ~, ~, ~, ~, ~, ~, cs_Resid, ~] = gmt_harmonic(time,[],cs_tmp);


for ii=1:size_a
    
    cs_res(:,:)=cs_Resid(:,:,ii);
    % SH to Grid
    if type==1
        c11cmn=[0.5 89.5 359.5 -89.5];
        lmcosi_fltr=gmt_cs2lmcosi(cs_res.^2); % change the CS format
        grid_tmp=plm2xyz_error(lmcosi_fltr,type,c11cmn);
    elseif type==0.25
        c11cmn=[0 90 359.75 -90];
        lmcosi_fltr=gmt_cs2lmcosi(cs_res.^2); % change the CS format
        grid_tmp=plm2xyz_error(lmcosi_fltr,type,c11cmn);
    elseif type==0.5
        c11cmn=[0.25 89.75 359.75 -89.75];
        lmcosi_fltr=gmt_cs2lmcosi(cs_res.^2); % change the CS format
        grid_tmp=plm2xyz_error(lmcosi_fltr,type,c11cmn);
    end
    grid_error(:,:,ii) = sqrt(grid_tmp(:,:));
end




end


