function [ days ] = gmt_get_days(int_year,mon,day )

% Converts the (year, month, day) to the days in the year
%
% INPUT:
%   int_year    Year
%   mon         Month
%   day         Day
% 
% OUTPUT:
%   days        Days in the year
%
% FENG Wei 05/07/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

d1=[0,31,59,90,120,151,181,212,243,273,304,334,365];
d2=[0,31,60,91,121,152,182,213,244,274,305,335,366]; % leap year

if int_year<0 || mon<0 || mon>12 || day<0 || day>31
    error('Wrong input parameter in gmt_get_days Function (Fengwei)!')
end

if ( mod(int_year,4)==0 &&  mod(int_year,400)~=0 ) || mod(int_year,400)==0 % leap year
    
    days = d2(mon)+day;
else
    if mon==2 && day>=29
        error('Wrong input parameter in get_days Function (Fengwei)!')
    end
    days = d1(mon)+day;
end

