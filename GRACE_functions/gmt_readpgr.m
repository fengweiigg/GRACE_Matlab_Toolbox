function [ grid_pgr ] = gmt_readpgr( dir_in )

% Read the GIA files and get the gridded field
%
% INPUT:
%   dir_in    filename with full path
% 
% OUTPUT:
%   grid_pgr  GIA grid
%
% FENG Wei 22/03/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

% read the header and fine the line number of head
h=0;
fid = fopen(dir_in,'r');
tline = fgetl(fid);
while strcmp(tline(1:3),'HDR')
    h = h+1;
    tline = fgetl(fid);
end
fclose(fid);


fid = fopen(dir_in,'r');
% read the head
for i=1:h
    tline = fgetl(fid);
end
% read the data
% lat: 89.5 ~ -89.5
% lon: -180 ~ 180
    
for i=1:360
        
    for j=1:180

        a = fscanf(fid,'%f %f',[1 3]);
        grid(i,j) = a(3);
    end
end

% change unit from mm/yr to m/yr
grid_pgr = flipud(grid')/1000;


end


