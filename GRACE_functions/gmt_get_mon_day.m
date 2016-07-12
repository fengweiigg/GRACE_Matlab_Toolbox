function [ mon,day ] = gmt_get_mon_day(int_year,days )

% Converts the days in the year to (year, month, day)
%
% INPUT:
%   int_year    Year
%   days        Days in the year
% 
% OUTPUT:
%   mon         Month
%   day         Day

% FENG Wei 05/07/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

d1=[0,31,59,90,120,151,181,212,243,273,304,334,365];
d2=[0,31,60,91,121,152,182,213,244,274,305,335,366]; % leap year

if int_year<0 || days<0 || days>366
    error('Wrong input parameter in get_mon_day Function!')
end

if days ==0
    mon =1;
    day =1;
    return
end

if ( mod(int_year,4)==0 &&  mod(int_year,400)~=0 ) || mod(int_year,400)==0 % leap year
    for i=1:12
        if days>d2(i) && days<=d2(i+1) 
            mon =i;
            day = days - d2(i);
            return
        end
    end
else
    for i=1:12
        if days>d1(i) && days<=d1(i+1)
            mon =i;
            day = days - d1(i);
            return
        end
    end
end
end

