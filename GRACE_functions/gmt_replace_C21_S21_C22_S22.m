function [cs_replace,tag] = gmt_replace_C21_S21_C22_S22(dir_in,cs,int_year,int_month,num_file)

% Read the C21, S21, C22, S22 files and replace the corresponding original
% Stokes coefficients with them
%
% Read C21_S21_RL05.txt or C22_S22_RL05.txt file provided from CSR's FTP site ftp://ftp.csr.utexas.edu/pub/slr/degree_2/
% Replace C21 & S21 or C22 & S22 in GRACE GSM files with those in
% C21_S21_RL05.txt or C22_S22_RL05.txt
% Restored AOD should be removed from values in C21_S21_RL05.txt and C22_S22_RL05.txt
%
% INPUT:
%   dir_in      full path
%   cs          spherical harmonic coefficients in CS-format
%   int_year    year
%   int_month   mont
%   num_file    number of files
% 
% OUTPUT:
%   cs_replace  spherical harmonic coefficients with replaced C21, S21, C22, S22
%   tag         check
%

%
% FENG Wei 22/03/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


cs_replace=cs;
time=zeros(num_file,1);
[~, FILE_NAME,~]=fileparts(dir_in);

if strcmp(FILE_NAME,'C21_S21_RL05') || strcmp(FILE_NAME,'C22_S22_RL05')
    tag=1;
    % read the header
    head_index=0;
    fid = fopen(dir_in,'r');
    tline = fgetl(fid);
    while strcmp(tline(1:2),'#')
        head_index = head_index+1;
        tline = fgetl(fid);
    end
    fclose(fid);
    % re-read the file, skip the comment lines
    [time_year, C21_C22, S21_S22, ~, ~, C21_C22_aod, S21_S22_aod ] = textread(dir_in,'%f %f %f %f %f %f %f %*[^\n]','headerlines',head_index);
    C21_C22_size=max(size(C21_C22));
    for ii=1:num_file
        % Approximate mid-point of monthly solution, See C21_S21_RL05.txt
        time(ii) = int_year(ii) + get_days(int_year(ii),int_month(ii),15)/365;
        for jj=1:C21_C22_size
            if abs(time(ii)-time_year(jj))<=0.03
                if strcmp(FILE_NAME,'C21_S21_RL05')
                    cs_replace(ii,3,2) = C21_C22(jj)-C21_C22_aod(jj)*(1E-10);% restored AOD has to be removed;
                    cs_replace(ii,1,3) = S21_S22(jj)-S21_S22_aod(jj)*(1E-10);
                elseif strcmp(FILE_NAME,'C22_S22_RL05')
                    cs_replace(ii,3,3) = C21_C22(jj)-C21_C22_aod(jj)*(1E-10);% restored AOD has to be removed
                    cs_replace(ii,2,3) = S21_S22(jj)-S21_S22_aod(jj)*(1E-10);
                end
            end
        end
    end
else
%     warndlg('The name of C21&S21 or C22&S22 file is wrong!','ERROR');
    tag=0;
end
end

