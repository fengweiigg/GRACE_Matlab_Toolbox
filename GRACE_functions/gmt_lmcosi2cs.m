function [cs]=gmt_lmcosi2cs(lmcosi)

% Transfer spherical harmonic coefficients from lmcosi format to CS format
% 
% INPUT:
%   lmcosi     Matrix listing l,m,cosine and sine spherical harmonic coefficients
%
% OUTPUT:
%   cs         C_lm & S_lm in CS format (spherical harmonic coefficients)
%
% FENG Wei 04/09/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


for ii=1:size(lmcosi,1)
    index_l=lmcosi(ii,1);
    index_m=lmcosi(ii,2);
    cs(index_l+1,index_m+1)=lmcosi(ii,3);%C_lm
    if index_m~=0
        cs(index_m,index_l+1)=lmcosi(ii,4);%S_lm
    end
    
end
end