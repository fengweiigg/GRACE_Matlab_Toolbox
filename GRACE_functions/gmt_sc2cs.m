function cs = gmt_sc2cs(sc)

% Transfer spherical harmonic coefficients from SC-format to CS-format
% 
% INPUT:
%   sc       C_lm & S_lm in SC format (spherical harmonic coefficients, /S|C\, (L+1)x(2L+1))
%
% OUTPUT:
%   cs       C_lm & S_lm in CS format (spherical harmonic coefficients, |C\S|, (L+1)x(L+1) matrix)
% 
% FENG Wei 05/09/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

if ndims(sc)==2 % sc is the one SC matrix
    [rows,cols] = size(sc);
    lmax = rows -1;
    if cols ~= 2*lmax+1, error('gmt_sc2cs: Matrix dimensions must be (L+1)x(2L+1).'), end
    c  = sc(:,lmax+1:2*lmax+1);
    s  = [zeros(lmax+1,1) sc(:,1:lmax)];
    cs = tril(c) + triu(rot90(s),1);
elseif ndims(sc)==3 % sc  is the series of SC matrixes
    [ntime,rows,cols] = size(sc);
    lmax = rows -1;
    if cols ~= 2*lmax+1, error('gmt_sc2cs: Matrix dimensions must be (L+1)x(2L+1).'), end
    for ii=1:ntime   
        sc_tmp(:,:)=sc(ii,:,:);
        c  = sc_tmp(:,lmax+1:2*lmax+1);
        s  = [zeros(lmax+1,1) sc_tmp(:,1:lmax)];
        cs_tmp = tril(c) + triu(rot90(s),1);
        cs(ii,:,:)=cs_tmp;
    end
else
    error('gmt_sc2cs: Input format is wrong.');
end
