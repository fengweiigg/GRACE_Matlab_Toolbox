function [lmcosi]=gmt_cs2lmcosi(cs)
%
% Transfer spherical harmonic coefficients from CS format to lmcosi format
% 
% INPUT:
%   cs         C_lm & S_lm in CS format (spherical harmonic coefficients)
% 
% OUTPUT:
%   lmcosi     Matrix listing l,m,cosine and sine spherical harmonic coefficients
% 

if ndims(cs)~=2 || size(cs,1)~=size(cs,2)
    error('gmt_cs2lmcosi: CS should be 2-D square matrix.');
end

degree_max=size(cs,1)-1;
index_lmcosi=(degree_max+2)*(degree_max+1)/2;
lmcosi=zeros(index_lmcosi,4);

for ii=1:degree_max+1
    for jj=1:ii
        index_tmp=ii*(ii-1)/2+jj;
        lmcosi(index_tmp,1)=ii-1; % l
        lmcosi(index_tmp,2)=jj-1; % m
        lmcosi(index_tmp,3)=cs(ii,jj); % C_lm
        if jj==1
            lmcosi(index_tmp,4)=0; % S_lm
        else
            lmcosi(index_tmp,4)=cs(jj-1,ii); % C_lm
        end
    end
end
end