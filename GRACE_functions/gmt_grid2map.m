function gmt_grid2map(grid_data,colorbar_value,colorbar_unit,title_string,filename)

% Plot global grid using M_map toolbox 
% 
% INPUT:
%   grid_data           Regular global gridded field (equi-angular grid N*2N)
%   colorbar_value      maximum value of colorbar
%   colorbar_unit       unit of colorbar
%   title_string        title above the figure
%   filename            output filename of figure
%  
% OUTPUT:
%
%
% FENG Wei 12/02/2016
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


if nargin<3 || nargin>5 
   error('number of input variables is wrong in function gmt_grid2map!');
end
if ndims(grid_data)~=2  % grid should be 2-D matrix
   error('The input grid should be 2-D matrix in function gmt_grid2map!');
end
    
[jj,ii] = size(grid_data); % jj must be lat, ii must be lon, kk must be time
if (jj==721 && ii==1440)
    lon=0:0.25:359.75;
    lat=90:-.25:-90;
elseif (jj==720 && ii==1440) % for altimetry data 0.25
    lon = 0.125:0.25:359.875;
    lat = 89.875:-0.25:-89.875;
elseif (jj==180 && ii==360)
    lon=0.5:359.5;
    lat=89.5:-1:-89.5;
elseif (jj==360 && ii==720)
    lon=0.25:0.5:359.75;
    lat=89.75:-0.5:-89.75;
elseif (jj==361 && ii==720)
    lon=0.25:0.5:359.75;
    lat=90:-0.5:-90;
elseif (jj==181 && ii==360) % for GIA grid file
    lon=0:360;
    lat=90:-1:-90;
else
    error('Wrong input grid in function gmt_grid2map_land!');
end

set(0,'DefaultFigureRenderer','zbuffer');

set(gcf,'Units','centimeters');
set(gcf,'Position',[10 10 40 20]);
set(gca,'position',[.05 .05, .9 .8])

[LON,LAT]=meshgrid(lon,lat);


m_proj('Equidistant Cylindrical','long',[0 360],'lat',[-90 90]);
% m_proj('Robinson','long',[0 360],'lat',[-90 90]);
m_pcolor(LON,LAT,grid_data);
shading flat;
m_coast('linewidth',1,'color','k');
% m_coast('patch',[0.8 0.8 0.8]);
%  m_gshhs_i('color','k');              % Coastline...
m_grid('xtick',6,'ytick',7,'tickdir','in','xlabeldir','middle',...
    'TickLength',0.008,'LineWidth',1.,'FontName', 'Helvetica','FontSize',15,'fontweight','bold');
if nargin>=4
    title(title_string,'fontsize',20,'FontName', 'Helvetica','fontweight','bold');
end
caxis([-colorbar_value,colorbar_value]);
% caxis([0,colorbar_value]);
% cptcmap('GMT_haxby')
% cptcmap('GMT_rainbow')
% cptcmap('GMT_wysiwyg')

colormap('jet')
h=colorbar('v','FontSize',20,'fontweight','bold');
set(get(h,'title'),'string',colorbar_unit,'FontSize',20,'fontweight','bold');

% h=contourcmap('jet',-colorbar_value:0.5:colorbar_value,'colorbar','on','FontSize',15,'fontweight','bold');

h_pos=get(h,'Position');
h_pos(1)=h_pos(1)+0.05;
h_pos(2)=h_pos(2)*1.2;
h_pos(3)=h_pos(3);
h_pos(4)=h_pos(4);
set(h,'Position',h_pos);



if nargin==5
set(gcf, 'Color', 'white'); % white bckgr
export_fig(gcf, ...      % figure handle
    filename,... % name of output file without extension
    '-painters', ...      % renderer
    '-pdf', ...           % file format
    '-r300' );             % resolution in dpi
end
end
