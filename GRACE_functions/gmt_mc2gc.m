function outcoeff = gmt_mc2gc(incoeff)

% Converts mass coefficients (mc) to geoid coefficients (gc)
% using the factor before summation and the degree-dependent factor in 
% equation (14) of Wahr et al. 1998
% divided by the density of water to convert surface load (in kg/m^2) to 
% equivalent water height in meter.
%
% INPUT:
%   incoeff    spherical harmonic coefficients (water equivalent height) (m)
% 
% OUTPUT:
%   outcoeff   spherical harmonic coefficients (geoid) in CS-format
%
% FENG Wei 18/04/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

[rows,cols] = size(incoeff);
if rows == cols					% field is in CS-format
   maxdeg  = rows - 1;
%    incoeff = cs2sc(incoeff);			
elseif cols-2*rows == -1			% field is in SC-format already
   maxdeg  = rows - 1;
   incoeff = gmt_sc2cs(incoeff);			% convert to CS-format   
else
   error('Check format of field.')
end


ae        =  6378136.3;             % semi-major axis of ellipsoid [m]
% GM        =  3.986004415e14;         % geocentric grav. constant [m^3 / s^2]

rho_w = 1000;       % density of water kg/m^3
rho_ave = 5517;     % average density of the earth kg/m^3
% dens = rho_ave/(3.*rho_w);
dens = ae*rho_ave/(3.*rho_w);
% Han and Wahr Love numbers
load 'HanWahrLoveNumbers.mat'

%Refer to Chamnbers deg1_coef.txt  Fengwei
k(2) = 0.021;

sc = gmt_cs2sc(incoeff);

% loop over degrees to multiply each row of coefficients in sc-format
% with the factor from equation (13) of Wahr et al. 1998
for ll = 0:1:maxdeg
    factor(ll+1)  = dens*( 2*ll+1 ) / ( 1+k(ll+1) );
    scnew(ll+1,:) = sc(ll+1,:)/factor(ll+1);
end
outcoeff = gmt_sc2cs(scnew);