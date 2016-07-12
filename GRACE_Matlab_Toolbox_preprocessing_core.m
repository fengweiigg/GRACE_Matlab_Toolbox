function GRACE_Matlab_Toolbox_preprocessing_core(controlfile_path)


% Read the Control File
fid=fopen(controlfile_path,'r');
num_file        = fscanf(fid,'%d',1);
radius_filter   = fscanf(fid,'%d',1);

destrip_method  = fscanf(fid,'%s',1);
% Check destriping method option
if ~strcmp(destrip_method,'NONE') ...
        && ~strcmp(destrip_method,'SWENSON') ...
        && ~strcmp(destrip_method,'CHAMBERS2007') ...
        && ~strcmp(destrip_method,'CHAMBERS2012') ...
        && ~strcmp(destrip_method,'CHENP3M6') ...
        && ~strcmp(destrip_method,'CHENP4M6') ...
        && ~strcmp(destrip_method,'DUAN')
    warndlg('Destrping method in the control file is wrong! (Options: NONE,SWENSON,CHAMBERS2007,CHAMBERS2012,CHENP3M4,CHENP4M6 and DUAN)','Warning');
    return;
end

option_gia      = fscanf(fid,'%s',1);
% Check GIA option
if ~strcmp(option_gia,'GIA_notRemoved') && ~strcmp(option_gia,'GIA_Removed_Geru')
    warndlg('GIA option in the control file is wrong! (Options: GIA_notRemoved and GIA_Removed_Geru)','Warning');
    return;
end

InputFileFormat = fscanf(fid,'%s',1);
% Check Input file format option
if ~strcmp(InputFileFormat,'GRACE') && ~strcmp(InputFileFormat,'ICGEM')
    warndlg('Format of input files in the control file is wrong! (Options: GRACE and ICGEM)','Warning');
    return;
end

OutputFileFormat_1= fscanf(fid,'%s',1);
% Check Output file format option
if strcmp(OutputFileFormat_1,'SH_MAT')
    OutputFileFormat_2= fscanf(fid,'%f',1);% max degree
    OutputFileFormat_3= fscanf(fid,'%s',1); % output name
elseif strcmp(OutputFileFormat_1,'GRID_MAT')
    OutputFileFormat_2= fscanf(fid,'%f',1); % max degree
    OutputFileFormat_3= fscanf(fid,'%s',1); % output name
    OutputFileFormat_4= fscanf(fid,'%f',1); % spatial resolution
    if OutputFileFormat_4~=1 && OutputFileFormat_4~=0.25
        warndlg('Resolution of output grids in the control file is wrong! (Options: 1 and 0.25)','Warning');
        return;
    end
else
    warndlg('Output format in the control file is wrong! (Options: SH_MAT and GRID_MAT)','Warning');
    return;
end

dir_c20         = fscanf(fid,'%s',1);
dir_c21_s21     = fscanf(fid,'%s',1);
dir_c22_s22     = fscanf(fid,'%s',1);
dir_degree_1    = fscanf(fid,'%s',1);
dir_in          = fscanf(fid,'%s',1);
dir_out         = fscanf(fid,'%s',1);
% Check files and directories
if ~exist(dir_c20,'file') && ~strcmp(dir_c20,'NAN')
    warndlg('Degree C20 file does not exist!','Warning');
    return;
end
if ~exist(dir_c21_s21,'file') && ~strcmp(dir_c21_s21,'NAN')
    warndlg('Degree C21 & S21 file does not exist!','Warning');
    return;
end
if ~exist(dir_c22_s22,'file') && ~strcmp(dir_c22_s22,'NAN')
    warndlg('Degree C22 & S22 file does not exist!','Warning');
    return;
end
if ~exist(dir_degree_1,'file') && ~strcmp(dir_degree_1,'NAN')
    warndlg('Degree 1 file does not exist!','Warning');
    return;
end
if ~exist(dir_in,'dir')
    warndlg('Directory of Input GRACE files does not exist!','Warning');
    return;
end
if ~exist(dir_out,'dir')
    warndlg('Directory of Output files does not exist!','Warning');
    return;
end

file_name=cell(num_file,1);% Name list of input GRACE files
for i = 1:num_file % read the name list of files in control file
    file_name{i} = fscanf(fid,'%s',1);
    if ~exist([dir_in,file_name{i}],'file')
        warndlg(strcat('Input GRACE files',[dir_in,file_name{i}],'does not exist!'),'Warning');
        return;
    end
end
fclose(fid);

%
% GSM-2_2003060-2003090_0030_UTCSR_0096_0005.gfc
% GSM-2_2004001-2004013_0013_UTCSR_0096_0005
% GSM-2_2002102-2002120_0018_UTCSR_0060_0004
%
% GSM-2_2002122-2002137_0013_EIGEN_G---_005a.gfc
% GSM-2_2003032-2003059_0028_EIGEN_G---_0005
% GSM-2_2003001-2003031_0029_JPLEM_0000_0005
% GSM-2_2002305-2002334_0028_JPLEM_0000_0005.gfc
%
%  ------------------------------------------------------
%   Initialize
%  ------------------------------------------------------

lmax=OutputFileFormat_2; % maximum degree of output files

cs              = zeros(num_file,lmax+1,lmax+1);
cs_sgi          = zeros(num_file,lmax+1,lmax+1);
int_year        = zeros(num_file,1);
int_month       = zeros(num_file,1);
time            = zeros(num_file,1);


cs_res          = zeros(num_file,lmax+1,lmax+1);
cs_destrip      = zeros(num_file,lmax+1,lmax+1);
cs_mss          = zeros(num_file,lmax+1,lmax+1);
cs_fltr         = zeros(num_file,lmax+1,lmax+1);

if strcmp(OutputFileFormat_1,'GRID_MAT')
    if strcmp(OutputFileFormat_4,'1')
        grid_data_grace = zeros(180,360,num_file);
    elseif strcmp(OutputFileFormat_4,'0.25')
        grid_data_grace = zeros(721,1440,num_file);
    end
end
%
% cs_tmp          = zeros(lmax+1,lmax+1);
% sc_tmp          = zeros(lmax+1,2*lmax+1);
% grid_tmp        = zeros(180,360);

% region_plot = zeros(num_file,1);
% region_plot_csr_error = zeros(num_file,1);


hwait=waitbar(0,'Waiting>>>>>>>>');
% ----------------------------------------
% Read the GRACE Files
% ----------------------------------------
if strcmp(InputFileFormat,'ICGEM')
    for ii=1:num_file
        waitbar(0.1*ii/num_file,hwait,'Read the GRACE Files');
        [cs(ii,:,:),cs_sgi(ii,:,:),int_year(ii),int_month(ii),time(ii)] = gmt_readgsm(dir_in,file_name{ii},lmax,'ICGEM');
    end
elseif strcmp(InputFileFormat,'GRACE')
    for ii=1:num_file
        waitbar(0.1*ii/num_file,hwait,'Read the GRACE Files');
        [cs(ii,:,:),cs_sgi(ii,:,:),int_year(ii),int_month(ii),time(ii)] = gmt_readgsm(dir_in,file_name{ii},lmax,'GRACE');
    end
end

% ----------------------------------------
% Read and Replace C20 and Degree_1
% ----------------------------------------
waitbar(0.2,hwait,'Replace Degree 2 and Degree 1');
pause(0.5);

cs_replace=cs;
if ~strcmp(dir_degree_1,'NAN')
    [cs_replace,tag] = gmt_replace_degree_1(dir_degree_1,cs_replace,int_year,int_month,num_file);
    if tag==0
        warndlg('Format of degree 1 file is wrong!','Warning');
        close(hwait);
        return;
    end
end
if ~strcmp(dir_c20,'NAN')
    [cs_replace,tag]  = gmt_replace_C20(dir_c20,cs_replace,int_year,int_month,num_file);
    if tag==0
        warndlg('Format of C20 file is wrong!','Warning');
        close(hwait);
        return;
    end
end
if ~strcmp(dir_c21_s21,'NAN')
    [cs_replace,tag] = gmt_replace_C21_S21_C22_S22(dir_c21_s21,cs_replace,int_year,int_month,num_file);
    if tag==0
        warndlg('Format of C21 & S21 file is wrong!','Warning');
        close(hwait);
        return;
    end
end
if ~strcmp(dir_c22_s22,'NAN')
    [cs_replace,tag] = gmt_replace_C21_S21_C22_S22(dir_c22_s22,cs_replace,int_year,int_month,num_file);
    if tag==0
        warndlg('Format of C22 & S22 file is wrong!','Warning');
        close(hwait);
        return;
    end
end


% ----------------------------------------
% Remove average
% ----------------------------------------
waitbar(0.3,hwait,'Remove the average gravity field');
pause(0.5);

if num_file>1
    cs_mean = mean(cs_replace);
    for i=1:num_file
        cs_res(i,:,:)  = cs_replace(i,:,:)-cs_mean(1,:,:);
    end
else % only one input files
    cs_res=cs_replace;
end

% ----------------------------------------
% Destriping
% ----------------------------------------
for i=1:num_file
    %     waitbar(0.3+0.4*i/num_file,hwait,'Do destriping');
    waitbar(0.3+0.5*i/num_file,hwait,'Convert geoid to mass & Do destriping');
    cs_tmp(:,:) = cs_res(i,:,:);
    
    % ----------------------------------------
    % Geoid to Mass
    % ----------------------------------------
    cs_tmp(:,:)=gmt_gc2mc(cs_tmp);

%     %     % TEST
%     if i==17
%         ttt(:,:)=gmt_gc2mc(cs_tmp);
%         cs_destrip(i,:,:) = gmt_destriping_test(ttt,destrip_method);
%     end
    
% Do destriping and get the mass coefficients
    cs_mss(i,:,:) = gmt_destriping(cs_tmp,destrip_method);

end
pause(0.5);

% figure(2)
% ttt1(:,:)=cs_res(17,:,:);
% ttt2=gmt_gc2mc(ttt1);
% plot(ttt2(:,2),'-s')
% 
% hold on
% ttt(:,:)=cs_destrip(17,:,:);
% plot(ttt(:,2),'-r')
% ----------------------------------------
% Geoid to Mass
% ----------------------------------------
% for i=1:num_file
%     waitbar(0.7+0.1*i/num_file,hwait,'Convert Gravity Geoid to Equivalent Water Height (unit:meter)');
%     cs_tmp(:,:) = cs_destrip(i,:,:);
%     cs_mss(i,:,:) = gmt_gc2mc(cs_tmp);
% end
% pause(0.5);

% figure(2)
% ttt1(:,:)=cs_res(17,:,:);
% ttt2=gmt_gc2mc(ttt1);
% plot(ttt2(:,2),'-s')
% 
% hold on
% ttt(:,:)=cs_mss(17,:,:);
% plot(ttt(:,2),'-r')


% disp('. remove Uplift PGR');
% grid_pgr  = readpgr_new( 'GIA_n100_uplift_0km.txt');

% ----------------------------------------
% Remove GIA effect
% ----------------------------------------
waitbar(0.8,hwait,'Remove the GIA effect');
if strcmp(option_gia,'GIA_Removed_Geru')
    % GIA file must be in Matlab Path
    grid_pgr  = gmt_readpgr('GIA_n100_mass_0km.txt'); 
    cs_pgr = gsha(grid_pgr,'mean','block',lmax);
    for ii=1:num_file
        cs_tmp(:,:)         = cs_mss(ii,:,:);
        % remove PGR from stokes coefficients
        cs_tmp              = cs_tmp - (time(ii)-mean(time))*cs_pgr;
        cs_mss(ii,:,:)      = cs_tmp(:,:);
    end
end
pause(0.5);

% ----------------------------------------
% Gaussian filtering
% ----------------------------------------        
waitbar(0.9,hwait,'Do Gaussian filtering');
cs_fltr=gmt_gaussian_filter(cs_mss,radius_filter);
pause(0.5);

% ----------------------------------------
% Save the output files
% ----------------------------------------
% from number to string
for ii=1:num_file
    str_year{ii}     = num2str(int_year(ii));
    str_month{ii}   = num2str(int_month(ii),'%02u');
end
% Save SH coefficients
if strcmp(OutputFileFormat_1,'SH_MAT')
    cs_grace=cs_fltr(:,1:OutputFileFormat_2+1,1:OutputFileFormat_2+1);
    save(strcat(dir_out, OutputFileFormat_3,'.mat'),'cs_grace','str_year','str_month','time');
    waitbar(1,hwait,'Save the Files');
    close(hwait);
    return;
end
% Create Grids
for ii=1:num_file
    cs_tmp(:,:)         = cs_fltr(ii,:,:);
    if OutputFileFormat_4==1
        [grid_tmp,~,~]= gshs(cs_tmp,'none','block',180,0,0.,0);
    elseif OutputFileFormat_4==0.25
        [grid_tmp,~,~]= gshs(cs_tmp,'none','pole',720,0,0.,0);
    end
    grid_data_grace(:,:,ii) = grid_tmp(:,:);
end
% Save Grids
if strcmp(OutputFileFormat_1,'GRID_MAT')    
    waitbar(1,hwait,'Save the grid file in MAT format');
    save(strcat(dir_out, OutputFileFormat_3,'.mat'),'grid_data_grace','str_year','str_month','time');
    close(hwait);
    return;
end
% if strcmp(OutputFileFormat_1,'GRID_GRD')
%     for ii=1:num_file
%         waitbar(0.9+0.1*ii/num_file,hwait,'Save the monthly grid files');
%         grid_tmp(:,:)=grid_grace(:,:,ii);
%         if OutputFileFormat_4==1
%             grdwrite2(0.5:359.5,-89.5:89.5,double(flipud(grid_tmp)),strcat(dir_out,'GSM_',str_year{ii},str_month{ii},'.grd'));
%         elseif OutputFileFormat_4==0.25
%             grdwrite2(0:0.25:359.75,-90:0.25:90,double(flipud(grid_tmp)),strcat(dir_out,'GSM_',str_year{ii},str_month{ii},'.grd'));
%         end
%     end
    
    %     write_fid = fopen(strcat(dir_out,'CSR_GSM_',str_year{ii},str_month{ii},'.txt'),'w');
    %     for j=1:180
    %         for i=1:360
    % %             fprintf(write_fid,'%4.1f %4.1f %f',lam(i),90.-theta(j),grid(j,i,ii)-(time(ii)-mean(time))*grid_pgr(j,i) );
    %             fprintf(write_fid,'%4.1f %4.1f %f',lam(i),90.-theta(j),grid(j,i,ii) );
    %             fprintf(write_fid,'\n');
    %         end
    %     end
    %     fclose(write_fid);
    
% end

close(hwait);