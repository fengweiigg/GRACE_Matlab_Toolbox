function [cs_replace,tag] = gmt_replace_C20(dir_in,cs,int_year,int_month,num_file)

% Read the C20 files and replace the corresponding original Stokes
% coefficients with them
%
% Read C20_RL05.txt file from CSR's FTP site ftp://ftp.csr.utexas.edu/pub/slr/degree_2/
% Or read TN-07_C20_SLR.txt from JPL PODAAC's FTP site
% Replace C20 in GRACE GSM files with those in C20_RL05.txt or TN-07_C20_SLR.txt
% Values in C20_RL05.txt and TN-07_C20_SLR.txt are the same! 
% Restored AOD should be removed from values in C20_RL05.txt
%
% INPUT:
%   dir_in      full path
%   cs          spherical harmonic coefficients in CS-format
%   int_year    year
%   int_month   mont
%   num_file    number of files
% 
% OUTPUT:
%   cs_replace  spherical harmonic coefficients with replaced C20
%   tag         check
%
% FENG Wei 22/03/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


cs_replace=cs;
time=zeros(num_file,1);
[~, FILE_NAME,~]=fileparts(dir_in);

if (strcmp(FILE_NAME,'TN-07_C20_SLR'))
    tag=1;
    % read the header
    head_index=0;
    fid = fopen(dir_in,'r');
    tline = fgetl(fid);
    while size(tline,2)<8 || ~strcmp(tline(1:8),'PRODUCT:')
        head_index = head_index+1;
        tline = fgetl(fid);
    end
    fclose(fid);
    if head_index<1
        warndlg('The C20 file is wrong!','Warning');
        return;
    end
    % re-read the file, skip the comment lines
    [~, time_year, C20, ~, ~] = textread(dir_in,'%f %f %f %f %f %*[^\n]','headerlines',head_index+1);
    C20_size=max(size(C20));
    for ii=1:num_file
        % at the beginning of month, See TN-07_C20_SLR.txt
        time(ii) = int_year(ii) + get_days(int_year(ii),int_month(ii),0)/365;
        for jj=1:C20_size
            if abs(time(ii)-time_year(jj))<=0.02
                cs_replace(ii,3,1) = C20(jj);
            end
        end
    end
    return;
elseif (strcmp(FILE_NAME,'C20_RL05'))
    tag=1;
    % read the header
    head_index=0;
    fid = fopen(dir_in,'r');
    tline = fgetl(fid);
    while strcmp(tline(1:1),'#')
        head_index = head_index+1;
        tline = fgetl(fid);
    end
    fclose(fid);    
    % re-read the file, skip the comment lines
    [time_year, C20, ~, ~, C20_aod] = textread(dir_in,'%f %f %f %f %f %*[^\n]','headerlines',head_index);
    C20_size=max(size(C20));
    for ii=1:num_file
        % Approximate mid-point of monthly solution, See C20_RL05.txt
        time(ii) = int_year(ii) + gmt_get_days(int_year(ii),int_month(ii),15)/365;
        for jj=1:C20_size
            if abs(time(ii)-time_year(jj))<=0.03
                cs_replace(ii,3,1) = C20(jj)-C20_aod(jj)*(1E-10); % restored AOD has to be removed
            end
        end
    end
    return;
else
    %     warndlg('The name of C20 file is wrong!','Warning');
    tag=0;
    return;
end
end
