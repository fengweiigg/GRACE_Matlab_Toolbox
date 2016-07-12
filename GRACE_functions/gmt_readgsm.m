function    [cs,cs_sigma,int_year,int_month,time] = gmt_readgsm(dir_in,file_name,lmax,type)

% Read the gravity filed files
%
% INPUT:
%   dir_in      full path
%   file_name   filename 
%   lmax        maximum degree in files
%   type        options: ICGEM, GRACE
% 
% OUTPUT:
%   cs          spherical harmonic coefficients in CS-format
%   cs_sigma    spherical harmonic coefficients in CS-format (formal error)
%   int_year    year
%   int_month   mont
%   time        fraction in the year
%
% FENG Wei 18/04/2016
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn


cs=zeros(lmax+1,lmax+1);
cs_sigma=zeros(lmax+1,lmax+1);

if  strcmp(type,'ICGEM')
    % read the header
    head_index=0;
    fid = fopen([dir_in file_name],'r');
    tline = fgetl(fid);
    while size(tline,2)<11 || ~strcmp(tline(1:11),'end_of_head')
        head_index = head_index+1;
        tline = fgetl(fid);
        if size(tline,2)>10 && strcmp(tline(1:10),'max_degree') % In JPL GSM files, the maximum degree changes, 60 or 90
            degree_max=str2double(tline(11:end));
        end
    end
    fclose(fid);
    % re-read the file, skip the comment lines
%     [name, l, m, Clm, Slm, Clm_sig, Slm_sig] = textread(strcat(dir_in,file_name),'%3s %6u %5u %20f %20f %12f %12f %*[^\n]','headerlines',head_index+1);
    [~, l, m, Clm, Slm, Clm_sigma, Slm_sigma] = textread(strcat(dir_in,file_name),'%s %u %u %f %f %f %f','headerlines',head_index+1);
    for i = 1:length(l)
        sc_tmp(l(i)+1,degree_max+1-m(i)) = Slm(i);
        sc_tmp(l(i)+1,degree_max+1+m(i)) = Clm(i);
        sc_sigma_tmp(l(i)+1,degree_max+1-m(i)) = Slm_sigma(i);
        sc_sigma_tmp(l(i)+1,degree_max+1+m(i)) = Clm_sigma(i);
    end
    cs_tmp=sc2cs(sc_tmp);
    cs_sigma_tmp=sc2cs(sc_sigma_tmp);
    % Get SH coefficients
    cs = cs_tmp(1:lmax+1,1:lmax+1);
    cs_sigma= cs_sigma_tmp(1:lmax+1,1:lmax+1);
    % Get time tag
    if strcmp(file_name(1:3),'GSM')||strcmp(file_name(1:3),'GAD')||strcmp(file_name(1:3),'GAC')||strcmp(file_name(1:3),'GAA')
        %GSM-2_2002123-2002137_0012_UTCSR_0096_0005.gfc
        year1 = str2num(file_name(7:10));
        year2 = str2num(file_name(15:18));
        day1  = str2num(file_name(11:13));
        day2  = str2num(file_name(19:21));
    elseif strcmp(file_name(1:11),'kfilter_DDK') %kfilter_DDK1_GSM-2_2002095-2002120_0021_UTCSR_0096_0005.gfc
        year1 = str2num(file_name(20:23));
        year2 = str2num(file_name(28:30));
        day1  = str2num(file_name(24:26));
        day2  = str2num(file_name(31:33));
    elseif strcmp(file_name(1:3),'ITG') %ITG-Grace2010_2003-06.gfc
        int_year = str2num(file_name(15:18));
        int_month = str2num(file_name(20:21));
        time=int_year+(int_month-1)*(1/12)+1/24;
        return;
    elseif strcmp(file_name(1:6),'IGGCAS') %IGGCAS_2004-01-01to2004-01-31recovered_gravity_model.txt
        int_year = str2num(file_name(8:11));
        int_month = str2num(file_name(13:14));
        time=int_year+(int_month-1)*(1/12)+1/24;
        return;
    elseif strcmp(file_name(1:3),'IGG') %IGG_sst_leo2008-01.gfc
        int_year = str2num(file_name(12:15));
        int_month = str2num(file_name(17:18));
        time=int_year+(int_month-1)*(1/12)+1/24;
        return;
    end
    if year1 == year2
        meanday = (day1+day2)/2;
    else
        if (day1+(366-day1+day2)/2)>365 % in latter year
            year1   = year1 + 1;
            meanday = day2-(366-day1+day2)/2;
        else
            meanday = day1+(366-day1+day2)/2;  % in former year
        end
    end
    time        = year1 + meanday/365.;
    meanday     = round(meanday);
    int_year    = year1;
    int_month   = gmt_get_mon_day(year1,meanday+1);
%     str_year     = num2str(year1);
%     str_month   = num2str(month,'%02u');
end

if  strcmp(type,'GRACE')
    % read the header
    head_index=0;
    fid = fopen([dir_in file_name],'r');
    tline = fgetl(fid);
    while size(tline,2)<6 || ~strcmp(tline(1:6),'GRCOF2')
        head_index = head_index+1;
        tline = fgetl(fid);
        if size(tline,2)>3 && strcmp(tline(1:3),'SHM') % In JPL GSM files, the maximum degree changes, 60 or 90
            degree_max=str2double(tline(9:11));
        end
    end
    fclose(fid);
    % re-read the file, skip the comment lines
    [~, l, m, Clm, Slm, Clm_sigma, Slm_sigma] = textread(strcat(dir_in,file_name),'%s %u %u %f %f %f %f %*[^\n]','headerlines',head_index);
    for i = 1:length(l)
        sc_tmp(l(i)+1,degree_max+1-m(i)) = Slm(i);
        sc_tmp(l(i)+1,degree_max+1+m(i)) = Clm(i);
        sc_sigma_tmp(l(i)+1,degree_max+1-m(i)) = Slm_sigma(i);
        sc_sigma_tmp(l(i)+1,degree_max+1+m(i)) = Clm_sigma(i);
    end
    cs_tmp=gmt_sc2cs(sc_tmp);
    cs_sigma_tmp=gmt_sc2cs(sc_sigma_tmp);
    % Get SH coefficients
    cs = cs_tmp(1:lmax+1,1:lmax+1);
    cs_sigma= cs_sigma_tmp(1:lmax+1,1:lmax+1);
    % Get time tag GSM-2_2003001-2003031_0026_EIGEN_G---_0005
    year1 = str2num(file_name(7:10));
    year2 = str2num(file_name(15:18));
    day1  = str2num(file_name(11:13));
    day2  = str2num(file_name(19:21));
    if year1 == year2
        meanday = (day1+day2)/2;
    else
        if (day1+(366-day1+day2)/2)>365 % in latter year
            year1   = year1 + 1;
            meanday = day2-(366-day1+day2)/2;
        else
            meanday = day1+(366-day1+day2)/2;  % in former year
        end
    end
    time        = year1 + meanday/365.;
    meanday     = round(meanday);
    int_year    = year1;
    int_month   = gmt_get_mon_day(year1,meanday+1);
end
end