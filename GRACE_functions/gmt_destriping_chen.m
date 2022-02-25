function scnew = gmt_destriping_chen(sc,type)

% Do destriping to spherical harmonic coefficients using Chen methods
% References:   Chen P3M6 2007 GRACE detects coseismic and postseismic deformation from the Sumatra-Andaman earthquake
%               Chen P4M6 2009 2005 drought event in the Amazon River basin as measured by GRACE and estimated by climate models
% 
% INPUT:
%   sc       C_lm & S_lm in SC format (original spherical harmonic coefficients, /S|C\, (L+1)x(2L+1))
%   type     destriping methods, options: CHENP3M6,CHENP4M6
%
% OUTPUT:
%   scnew    C_lm & S_lm in CS format (destriped spherical harmonic coefficients, /S|C\, (L+1)x(2L+1))
% 
% FENG Wei 22/03/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


% sc represents the original sphere harmonics before destriping
[rows,cols] = size(sc);
if rows == cols					% field is in CS-format
    maxdeg  = rows - 1;
    sc = cs2sc(sc);			% convert to SC-format
elseif cols-2*rows == -1			% field is in SC-format already
    maxdeg  = rows - 1;
else
    error('Check format of gravity field data.')
end

if strcmp(type,'CHENP3M6')
    poly_ord    = 3;
    start_ord   = 6;
elseif strcmp(type,'CHENP4M6')
    poly_ord    = 4;
    start_ord   = 6;
else
    error('The Chen Destriping method is wrong!');
end

if maxdeg==60 % old release of RL05 CSR GSM files
    end_ord     = 50; 
elseif maxdeg==96 % new release of RL05 CSR GSM files
    end_ord     = 86;
elseif maxdeg==90 % RL05 JPL/GFZ GSM files
    end_ord     = 80;
else
    error('The maximum degree of GRACE gravity field should be 60/90/96!');
end


sc_s = zeros(maxdeg+1,2*maxdeg+1);

for mm = start_ord:end_ord
    
    % for mm = start_ord:end_ord
    clm_col_s = zeros(maxdeg+1,1);
    slm_col_s = zeros(maxdeg+1,1);
    
        % even(odd) degree
        n=mm:2:maxdeg;
        clm = sc(n+1,maxdeg+1+mm)';
        slm = sc(n+1,maxdeg+1-mm)';
        % fit a  polynomial
        % 02/25/2022, Yu Zhang, use centering and scaling in fitting,
        % otherwise Warning: Polynomial is badly conditioned for P4M6
        [poly_clm,~,mu] = polyfit(n,clm,poly_ord);
        clm_col_s(n+1,1) = polyval(poly_clm,n,[],mu);
        [poly_slm,~,mu] = polyfit(n,slm,poly_ord);
        slm_col_s(n+1,1) = polyval(poly_slm,n,[],mu);
               
        % odd(even) degree
        n=mm+1:2:maxdeg;
        clm = sc(n+1,maxdeg+1+mm)';
        slm = sc(n+1,maxdeg+1-mm)';
        % fit a  polynomial
        [poly_clm,~,mu] = polyfit(n,clm,poly_ord);
        clm_col_s(n+1,1) = polyval(poly_clm,n,[],mu);
        [poly_slm,~,mu] = polyfit(n,slm,poly_ord);
        slm_col_s(n+1,1) = polyval(poly_slm,n,[],mu);
        
        sc_s(:,maxdeg+1+mm) = clm_col_s;
        sc_s(:,maxdeg+1-mm) = slm_col_s;    

end



scnew = sc-sc_s;

for mm = end_ord+1:maxdeg
    scnew(:,maxdeg+1+mm) = 0.;
    scnew(:,maxdeg+1-mm) = 0.;
end

end


