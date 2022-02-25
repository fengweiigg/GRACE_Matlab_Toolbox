function [cs_destrip] = gmt_destriping(field,destrip_method)

% Do destriping to spherical harmonic coefficients using different methods
% 
% INPUT:
%   field               C_lm & S_lm in CS or SC format (spherical harmonic coefficients)
%   destrip_method      destriping methods: options: NONE, SWENSON, CHAMBERS2007,CHAMBERS2012, CHENP3M6, CHENP4M6, DUAN
% 
% OUTPUT:
%   cs_destrip          destriped-removed spherical harmonic coefficients
% 
% FENG Wei 05/09/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


%
% change to SC format, if the input is CS format
%
[rows,cols] = size(field);
if rows == cols					% field is in CS-format
   field = gmt_cs2sc(field);			% convert to SC-format
elseif cols-2*rows ~= -1			% field is not in SC-format
   error('Check format of gravity field data.')
end

%
% destrip the gravity field (GRACE level-2 GSM files)
%
if strcmp(destrip_method,'SWENSON')
    sc_destrip = gmt_destriping_swenson(field);
    cs_destrip=gmt_sc2cs(sc_destrip);
elseif strcmp(destrip_method,'CHAMBERS2007')
    sc_destrip = gmt_destriping_chambers(field,'CHAMBERS2007');
    cs_destrip=gmt_sc2cs(sc_destrip);
elseif strcmp(destrip_method,'CHAMBERS2012')
    sc_destrip = gmt_destriping_chambers(field,'CHAMBERS2012');
    cs_destrip=gmt_sc2cs(sc_destrip);
elseif strcmp(destrip_method,'CHENP3M6')
    sc_destrip = gmt_destriping_chen(field,'CHENP3M6');
    cs_destrip=gmt_sc2cs(sc_destrip);
elseif strcmp(destrip_method,'CHENP4M6')
    sc_destrip = gmt_destriping_chen(field,'CHENP4M6'); % 02/22/2022, Yu Zhang, fix a bug here
    cs_destrip=gmt_sc2cs(sc_destrip);
elseif strcmp(destrip_method,'DUAN')
    pair1=35;
    pair2=15;
    r=3.5;
    gamma=0.1;
    p=2;
    K=15;
    A=30;
    sc_destrip = gmt_destriping_duan(pair1,pair2,r,gamma,p,K,A,field);
    cs_destrip=gmt_sc2cs(sc_destrip);
elseif strcmp(destrip_method,'NONE')
    sc_destrip = field;
    cs_destrip=gmt_sc2cs(sc_destrip);
    else
    error('gmt_destriping: Destrping method is wrong! (Options: NONE,SWENSON,CHAMBERS2007,CHAMBERS2012,CHENP3M4,CHENP4M6 and DUAN)');
end

end