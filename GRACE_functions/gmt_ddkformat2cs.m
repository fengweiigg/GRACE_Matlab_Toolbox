function cs = gmt_ddkformat2cs(shc_ddkformat)

% Transfer spherical harmonic coefficients from DDK-format to CS-format
% 
% INPUT:
%   shc_ddkformat   C_lm & S_lm in DDK format (spherical harmonic coefficients)
% 
% OUTPUT:
%   cs              C_lm & S_lm in CS format (spherical harmonic coefficients, |C\S|, (L+1)x(L+1) matrix)
%
% FENG Wei 09/07/2016
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


lmax=size(shc_ddkformat.C,1)-1;
ntime=size(shc_ddkformat.C,3);
for ii=1:ntime
    sc_tmp=zeros(lmax+1,2*lmax+1);
    sc_tmp(:,lmax+1:2*lmax+1)=shc_ddkformat.C(:,:,ii);
    sc_tmp(2:lmax+1,1:lmax)=fliplr(shc_ddkformat.S(2:lmax+1,2:lmax+1,ii));
    cs_tmp=gmt_sc2cs(sc_tmp);
    cs(ii,:,:)=cs_tmp;
end
cs=squeeze(cs);% if only one month CS data
end


