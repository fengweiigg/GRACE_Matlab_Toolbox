function sc = gmt_cs2sc(cs)

% Transfer spherical harmonic coefficients from CS-format to SC-format
% 
% INPUT:
%   cs       C_lm & S_lm in CS format (spherical harmonic coefficients, |C\S|, (L+1)x(L+1) matrix)
% 
% OUTPUT:
%   sc       C_lm & S_lm in SC format (spherical harmonic coefficients, /S|C\, (L+1)x(2L+1))
%
% FENG Wei 05/09/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn
if ndims(cs)==2 % cs is the one CS matrix
    [rows,cols] = size(cs);
    lmax = rows -1;
    if cols ~= rows, error('gmt_cs2sc: A square matrix is needed.'), end
    c  = tril(cs);
    s  = rot90(triu(cs,1),-1);
    sc = [s(:,2:lmax+1) c];
elseif ndims(cs)==3  % cs is the series of CS matrixes
    [ntime,rows,cols] = size(cs);
    if cols ~= rows, error('gmt_cs2sc: A square matrix is needed.'), end
    for ii=1:ntime
        cs_tmp(:,:)=cs(ii,:,:);
        c  = tril(cs_tmp);
        s  = rot90(triu(cs_tmp,1),-1);
        sc_tmp = [s(:,2:lmax+1) c];
        sc(ii,:,:)=sc_tmp;
    end
else
    error('gmt_cs2sc: Input format is wrong.');
end