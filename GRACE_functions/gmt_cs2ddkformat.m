function shc_ddkformat = gmt_cs2ddkformat(cs)

% Transfer spherical harmonic coefficients from CS-format to DDK-format
%
% INPUT:
%   cs              C_lm & S_lm in CS format (spherical harmonic coefficients, |C\S|, (L+1)x(L+1) matrix)
%
% OUTPUT:
%   shc_ddkformat   C_lm & S_lm in DDK format (spherical harmonic coefficients)
%
% FENG Wei 09/07/2016
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

if ndims(cs) == 3 % cs is the series of CS matrixes
    [ntime,rows,cols] = size(cs);
elseif ndims(cs) == 2 % cs is the one CS matrix
    [rows,cols] = size(cs);
    ntime=1;
else
    error('gmt_cs2ddkformat: Dimension of matrix is wrong.');
end

lmax = rows -1;
if cols ~= rows, error('gmt_cs2ddkformat: A CS-format matrix is needed.'), end
shc_ddkformat.C=zeros(lmax+1,lmax+1,ntime);
shc_ddkformat.S=zeros(lmax+1,lmax+1,ntime);


for ii=1:ntime
    if ndims(cs) == 3
        cs_tmp(:,:)=cs(ii,:,:);
    elseif ndims(cs) == 2
        cs_tmp=cs;
    end
    % transfer to SC format firstly, /S|C\, (L+1)x(2L+1)
    c  = tril(cs_tmp);
    s  = rot90(triu(cs_tmp,1),-1);
    sc_tmp = [s(:,2:lmax+1) c];
    
    shc_ddkformat.C(:,:,ii)=sc_tmp(:,lmax+1:2*lmax+1);
    
    % tmp=fliplr(sc_tmp(:,1:lamx));
    
    shc_ddkformat.S(2:lmax+1,2:lmax+1,ii)=fliplr(sc_tmp(2:lmax+1,1:lmax));
end

end


