function [plot_region]=gmt_grid2series_area(grid,dir_msk,type,bound)

% Calculate the area of specified region with valide values
% 
% INPUT:
%   grid            regular global gridded field
%                   size of grid must be: 720*1440, 721*1440, 180*360, 181*360, 360*720, 361*720
%   dir_msk         boundary/mask file
%   type            options: line, mask
%   bound           the maximum latitude boundary, default is 90 degree
%  
% OUTPUT:
%   plot_region     get the area of specified region
%
% FENG Wei 12/02/2016
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn



if ndims(grid)==2
    [jj,ii] = size(grid); % jj must be lat, ii must be lon
    kk=1;
elseif ndims(grid)==3
    [jj,ii,kk] = size(grid); % jj must be lat, ii must be lon, kk must be time
else
    warndlg('Input grid file should be 3-D/2-D matrix, i.e., grid(lat,lon,time)/grid(lat,lon) in gmt_grid2series_area fucntion','Warning');
    return;
end

plot_region = zeros(kk,1);

if nargin==3
    bound=90; % bound represents the range of latitude, used for type=='mask'
end

%
% create the mask, if the input is boundary file
%
if (strcmp(type,'line'))
    
    [~,~,FILE_TYPE]=fileparts(dir_msk);
    if ~strcmp(FILE_TYPE,'.bln') % if input is .mat mask file
        warndlg('Input boundary line file should be *.bln in gmt_grid2series function!','Warning');
        return;
    end
    % Read the line file
    fid           = fopen(dir_msk,'r');
    num_points    = fscanf(fid,'%d',1);
    xv            = zeros(num_points,1);
    yv            = zeros(num_points,1);
    for i=1:num_points
        a       = fscanf(fid,'%f %f',[1 2]);
        xv(i)   = a(1);
        yv(i)   = a(2);
    end
    fclose(fid);
    lon_min=min(xv);
    lon_max=max(xv);
    lat_min=min(yv);
    lat_max=max(yv);
    % Create the mask grid for different kinds of grid file
    if (jj==721 && ii==1440)
        lon=0:0.25:359.75;
        lat=90:-.25:-90;
        lg_msk = zeros(jj,ii);
        angle_area = zeros(jj,ii);
        for j=1:jj
            for i=1:ii
                if lon(i)>=lon_min && lon(i)<=lon_max && lat(j)>=lat_min && lat(j)<=lat_max
                    in =inpolygon(lon(i),lat(j),xv,yv);
                    if in==1
                        lg_msk(j,i) = 1;
                        angle_area(j,i) = cos(lat(j)*pi/180)*(0.25*pi/180)*(0.25*pi/180);% get the angle_area
                    end
                end
            end
        end
    elseif (jj==720 && ii==1440) % for altimetry data 0.25
        lon = 0.125:0.25:359.875;
        lat = 89.875:-0.25:-89.875;
        lg_msk = zeros(jj,ii);
        angle_area = zeros(jj,ii);
        for j=1:jj
            for i=1:ii
                if lon(i)>=lon_min && lon(i)<=lon_max && lat(j)>=lat_min && lat(j)<=lat_max
                    in =inpolygon(lon(i),lat(j),xv,yv);
                    if in==1
                        lg_msk(j,i) = 1;
                        angle_area(j,i) = cos(lat(j)*pi/180)*(0.25*pi/180)*(0.25*pi/180);% get the angle_area
                    end
                end
            end
        end
    elseif (jj==180 && ii==360)
        lon=0.5:359.5;
        lat=89.5:-1:-89.5;
        lg_msk = zeros(jj,ii);
        angle_area = zeros(jj,ii);
        for j=1:jj
            for i=1:ii
                if lon(i)>=lon_min && lon(i)<=lon_max && lat(j)>=lat_min && lat(j)<=lat_max
                    in =inpolygon(lon(i),lat(j),xv,yv);
                    if in==1
                        lg_msk(j,i) = 1;
                        angle_area(j,i) = cos(lat(j)*pi/180)*(1*pi/180)*(1*pi/180);% get the angle_area
                    end
                end
            end
        end        
    elseif (jj==360 && ii==720)
        lon=0.25:0.5:359.75;
        lat=89.75:-0.5:-89.75;
        lg_msk = zeros(jj,ii);
        angle_area = zeros(jj,ii);
        for j=1:jj
            for i=1:ii
                if lon(i)>=lon_min && lon(i)<=lon_max && lat(j)>=lat_min && lat(j)<=lat_max
                    in =inpolygon(lon(i),lat(j),xv,yv);
                    if in==1
                        lg_msk(j,i) = 1;
                        angle_area(j,i) = cos(lat(j)*pi/180)*(0.5*pi/180)*(0.5*pi/180);% get the angle_area
                    end
                end
            end
        end
    elseif (jj==361 && ii==720)
        lon=0.25:0.5:359.75;
        lat=90:-0.5:-90;
        lg_msk = zeros(jj,ii);
        angle_area = zeros(jj,ii);
        for j=1:jj
            for i=1:ii
                if lon(i)>=lon_min && lon(i)<=lon_max && lat(j)>=lat_min && lat(j)<=lat_max
                    in =inpolygon(lon(i),lat(j),xv,yv);
                    if in==1
                        lg_msk(j,i) = 1;
                        angle_area(j,i) = cos(lat(j)*pi/180)*(0.5*pi/180)*(0.5*pi/180);% get the angle_area
                    end
                end
            end
        end
    elseif (jj==181 && ii==360) % for GIA grid file
        lon=0:360;
        lat=90:-1:-90;
        lg_msk = zeros(jj,ii);
        angle_area = zeros(jj,ii);
        for j=1:jj
            for i=1:ii
                if lon(i)>=lon_min && lon(i)<=lon_max && lat(j)>=lat_min && lat(j)<=lat_max
                    in =inpolygon(lon(i),lat(j),xv,yv);
                    if in==1
                        lg_msk(j,i) = 1;
                        angle_area(j,i) = cos(lat(j)*pi/180)*(1*pi/180)*(1*pi/180);% get the angle_area
                    end
                end
            end
        end
    else
        warndlg('ERROR in Input boundary line file in gmt_grid2series function','Warning');
        return
    end
end


if (strcmp(type,'mask'))    
    [~,~,FILE_TYPE]=fileparts(dir_msk);
    if strcmp(FILE_TYPE,'.mat') % if input is .mat mask file
        load(dir_msk); % load the lg_msk variable in dir_msk file
        if ~exist('lg_msk','var')
            warnlg('''lg_msk'' variable is not in the input mat file. ERROR in gmt_grid2series function','Warning');
            return;
        end
        [jjj,iii]=size(lg_msk);
        if jj~=jjj || ii~=iii
            warnlg('Size of''lg_msk'' variable is wrong. ERROR in gmt_grid2series function','Warning');
            return;
        end
        if (jjj==180 && iii==360)
            lon=0.5:359.5;
            lat=89.5:-1:-89.5;
        end
        if (jjj==721 && iii==1440)
            lon=0:0.25:359.75;
            lat=90:-.25:-90;
        end
    elseif strcmp(FILE_TYPE,'.xyz') % if input is .xyz mask file
        if (jj==180 && ii==360)
            lon=0.5:359.5;
            lat=89.5:-1:-89.5;
            lg_msk = zeros(jj,ii);
            fid    = fopen(dir_msk,'r');
            angle_area = zeros(jj,ii);
            for j=1:jj
                for i=1:ii
                    a       = fscanf(fid,'%f %f',[1 3]);
                    lg_msk(j,i)   = a(3); %%if lg_msk(j,i)==0 % ocean
                    angle_area(j,i) = cos(lat(j)*pi/180)*(1*pi/180)*(1*pi/180);% get the angle_area
                end
            end
            %         lg_msk = ~lg_msk; % set the ocean to 1, for compute sea level change
            % set the area outside the bound to ZERO
            for j=1:jj
                if abs(lat(j))>bound
                    lg_msk(j,:)=1;
                end
            end
        end
        % Be careful HERE! 
        % To calculate global mean sea level, I set ocean to ONES, land
        % ZEROs, FENG Wei
        lg_msk=~logical(lg_msk);
        
%         elseif (jj==721 && ii==1440)
%             lon=0:0.25:359.75;
%             lat=90:-.25:-90;
%             lg_msk = zeros(jj,ii);
%             fid    = fopen(dir_msk,'r');
%             for j=1:jj
%                 for i=1:ii
%                     a       = fscanf(fid,'%f %f',[1 3]);
%                     lg_msk(j,i)   = a(3); %%if lg_msk(j,i)==0 % ocean
%                 end
%             end
%         end
    else
        error('ERROR in gmt_grid2series_area function');
    end
end

% %
% % create the weight
% %
% grid_weight=zeros(jj,ii);
% for j=1:jj
%     for i=1:ii
%         grid_weight(j,i)=cos(lat(j)*pi/180);
%     end
% end
% grid_weight=grid_weight.*lg_msk;

grid_weight=angle_area.*lg_msk;

ae        =  6378136.3;             % semi-major axis of ellipsoid [m]

%
% create the time series
%
tmmp=zeros(jj,ii);
% grid(logical(isnan(grid)))=0; % set possible NaNs to zeros
for k=1:kk
    tmmp(:,:)=grid(:,:,k);
    tmmp(logical(~isnan(tmmp)))=1; % set the point values to zeros

    tmmppp=tmmp.*grid_weight;
    
    tmmppp(logical(isnan(tmmppp)))=0;% there are NaNs outside the basin
    
%     plot_region(k)=sum(sum(tmmppp))* ae* ae;% unit is m^2; % sum all points angle area in the basin
    % be careful, which kind of unit you need!
    plot_region(k)=sum(sum(tmmppp))* ae* ae/10^6;% unit is km^2; % sum all points angle area in the basin
    
end





end