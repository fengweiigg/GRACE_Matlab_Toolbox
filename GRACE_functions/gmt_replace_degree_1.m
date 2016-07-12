function [ cs_replace,tag ] = gmt_replace_degree_1(dir_in,cs,int_year,int_month,num_file )

% Read the degree 1 files and replace the corresponding original Stokes
% coefficients with them
%
% INPUT:
%   dir_in      full path
%   cs          spherical harmonic coefficients in CS-format
%   int_year    year
%   int_month   mont
%   num_file    number of files
% 
% OUTPUT:
%   cs_replace  spherical harmonic coefficients with replaced degree one
%   tag         check
%
% FENG Wei 22/03/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


cs_replace=zeros(size(cs));
[~, FILE_NAME,~]=fileparts(dir_in);

if (strcmp(FILE_NAME,'deg1_coef'))
    tag=1;
    % read the header
    head_index=0;
    fid = fopen(dir_in,'r');
    tline = fgetl(fid);
    while size(tline,2)<32 || ~strcmp(tline(1:32),'''(a6,2i5,2E19.12,2E11.4,2f14.4)''')
        head_index = head_index+1;
        tline = fgetl(fid);
    end
    fclose(fid);

    % re-read the file, skip the comment lines
    [year,mon, l, m, aa, bb, ~, ~] = textread(dir_in,'%4d %2d %5d %5d %19f %19f %11f %11f %*[^\n]','headerlines',head_index+1);
    
    % Replace Degree 1
    for ii=1:num_file
        cs_replace(ii,:,:) = cs(ii,:,:);
        for jj=1:length(l)
            if (int_year(ii)==year(jj) && int_month(ii)==mon(jj) && l(jj)==1 && m(jj)==0)
                cs_replace(ii,2,1) = aa(jj);
            elseif (int_year(ii)==year(jj) && int_month(ii)==mon(jj) && l(jj)==1 && m(jj)==1)
                cs_replace(ii,2,2) = aa(jj);
                cs_replace(ii,1,2) = bb(jj);
            end
        end
    end
else
    tag=0;
    %     warndlg('Format of degree 1 file is wrong!','Warning');
end

