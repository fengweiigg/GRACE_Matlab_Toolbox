function varargout = GRACE_Matlab_Toolbox_preprocessing(varargin)
% GRACE_MATLAB_TOOLBOX_PREPROCESSING MATLAB code for GRACE_Matlab_Toolbox_preprocessing.fig
%      GRACE_MATLAB_TOOLBOX_PREPROCESSING, by itself, creates a new GRACE_MATLAB_TOOLBOX_PREPROCESSING or raises the existing
%      singleton*.
%
%      H = GRACE_MATLAB_TOOLBOX_PREPROCESSING returns the handle to a new GRACE_MATLAB_TOOLBOX_PREPROCESSING or the handle to
%      the existing singleton*.
%
%      GRACE_MATLAB_TOOLBOX_PREPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRACE_MATLAB_TOOLBOX_PREPROCESSING.M with the given input arguments.
%
%      GRACE_MATLAB_TOOLBOX_PREPROCESSING('Property','Value',...) creates a new GRACE_MATLAB_TOOLBOX_PREPROCESSING or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GRACE_Matlab_Toolbox_preprocessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GRACE_Matlab_Toolbox_preprocessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GRACE_Matlab_Toolbox_preprocessing

% Last Modified by GUIDE v2.5 23-Mar-2015 18:15:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GRACE_Matlab_Toolbox_preprocessing_OpeningFcn, ...
    'gui_OutputFcn',  @GRACE_Matlab_Toolbox_preprocessing_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before GRACE_Matlab_Toolbox_preprocessing is made visible.
function GRACE_Matlab_Toolbox_preprocessing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GRACE_Matlab_Toolbox_preprocessing (see VARARGIN)

% Choose default command line output for GRACE_Matlab_Toolbox_preprocessing
handles.output = hObject;

guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% --- Outputs from this function are returned to the command line.
function varargout = GRACE_Matlab_Toolbox_preprocessing_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end
% Update handles structure
guidata(handles.figure1, handles);

% --- Main Function in this Toolbox ---
function PushbuttonCalculate_Callback(hObject, eventdata, handles)
if ~isempty(get(handles.EditOpenControlFile,'String'))
    GRACE_Matlab_Toolbox_preprocessing_core(get(handles.EditOpenControlFile,'String'));
    return;
elseif ~isempty(get(handles.EditSaveControlFile,'String'))
    GRACE_Matlab_Toolbox_preprocessing_core(get(handles.EditSaveControlFile,'String'));
    return;
else
    warndlg('Control file should be specified in STEP3/STEP1!','Warning');
end

% Specify the GRACE Level2 Files
function PushButtonOpenFiles_Callback(hObject, eventdata, handles)
% Open GRACE GSM Files
[InputFilename, InputPathname, filterindex] = uigetfile( ...
    {'*','GRACE-files (*)'; ...
    '*.gfc','ICGEM-files (*.gfc)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Pick GRACE Level2 files', ...
    'MultiSelect', 'on');

if filterindex==1
    set(handles.InputFileListbox,'string',InputFilename);% Show the Files in the listbox
    handles.InputFilename=InputFilename;
    handles.InputPathname=InputPathname;
    guidata(hObject,handles);
else
    warndlg('Input files should be specified!!!','Warning');
    return;
end
function InputFileListbox_Callback(hObject, eventdata, handles)
function InputFileListbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Specify the Maximum Degree
function EditMaxDegreeOutput_Callback(hObject, eventdata, handles)
handles.MaxDegreeOutput=str2double(get(hObject,'String'));
guidata(hObject, handles);
function EditMaxDegreeOutput_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Specify the Gaussin Filter Radius
function EditFilterRadius_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
function EditFilterRadius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Specify the path of C20, C21, S21, C22, S22
function CheckboxC20_Callback(hObject, eventdata, handles)
if get(handles.CheckboxC20,'value')==1
    [InputFilename, InputPathname, filterindex] = uigetfile( ...
        {'*.*',  'All Files (*.*)'}, ...
        'Pick C20 file');
    if filterindex==1 % file selected
        % Show the Directory and File in the textbox
        set(handles.EditC20, 'String', [InputPathname,InputFilename]);
        handles.PathC20=[InputPathname,InputFilename];
        guidata(hObject, handles);
    else
        set(handles.CheckboxC20,'value',0);
        set(handles.EditC20, 'String', 'NAN');
    end
else
    set(handles.EditC20, 'String', 'NAN');
end
function CheckboxC21S21_Callback(hObject, eventdata, handles)
if get(handles.CheckboxC21S21,'value')==1
    [InputFilename, InputPathname, filterindex] = uigetfile( ...
        {'*.*',  'All Files (*.*)'}, ...
        'Pick C21 & S21 file');
    if filterindex==1 % file selected
        % Show the Directory and File in the textbox
        set(handles.EditC21S21, 'String', [InputPathname,InputFilename]);
        handles.PathC21S21=[InputPathname,InputFilename];
        guidata(hObject, handles);
    else
        set(handles.CheckboxC21S21,'value',0);
        set(handles.EditC21S21, 'String', 'NAN');
    end
else
    set(handles.EditC21S21, 'String', 'NAN');
end
function CheckboxC22S22_Callback(hObject, eventdata, handles)
if get(handles.CheckboxC22S22,'value')==1
    [InputFilename, InputPathname, filterindex] = uigetfile( ...
        {'*.*',  'All Files (*.*)'}, ...
        'Pick C22 & S22 file');
    if filterindex==1 % file selected
        % Show the Directory and File in the textbox
        set(handles.EditC22S22, 'String', [InputPathname,InputFilename]);
        handles.PathC22S22=[InputPathname,InputFilename];
        guidata(hObject, handles);
    else
        set(handles.CheckboxC22S22,'value',0);
        set(handles.EditC22S22, 'String', 'NAN');
    end
else
    set(handles.EditC22S22, 'String', 'NAN');
end
function CheckboxDegree1_Callback(hObject, eventdata, handles)
if get(handles.CheckboxDegree1,'value')==1
    [InputFilename, InputPathname, filterindex] = uigetfile( ...
        {'*.*',  'All Files (*.*)'}, ...
        'Pick Degree 1 file');
    if filterindex==1 % file selected
        % Show the Directory and File in the textbox
        set(handles.EditDegree1, 'String', [InputPathname,InputFilename]);
        handles.PathDegree1=[InputPathname,InputFilename];
        guidata(hObject, handles);
    else
        set(handles.CheckboxDegree1,'value',0);
        set(handles.EditDegree1, 'String', 'NAN');
    end
else
    set(handles.EditDegree1, 'String', 'NAN');
end

function EditC20_Callback(hObject, eventdata, handles)
function EditC20_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function EditC21S21_Callback(hObject, eventdata, handles)
function EditC21S21_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function EditC22S22_Callback(hObject, eventdata, handles)
function EditC22S22_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function EditDegree1_Callback(hObject, eventdata, handles)
function EditDegree1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Specify the Destriping method
function RadiobuttonNonedestriping_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonNonedestriping,'value',1);
set(handles.RadiobuttonSwenson,'value',0);
set(handles.RadiobuttonChambers2007,'value',0);
set(handles.RadiobuttonChambers2012,'value',0);
set(handles.RadiobuttonChenP3M6,'value',0);
set(handles.RadiobuttonChenP4M6,'value',0);
set(handles.RadiobuttonDuan,'value',0);
function RadiobuttonSwenson_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonNonedestriping,'value',0);
set(handles.RadiobuttonSwenson,'value',1);
set(handles.RadiobuttonChambers2007,'value',0);
set(handles.RadiobuttonChambers2012,'value',0);
set(handles.RadiobuttonChenP3M6,'value',0);
set(handles.RadiobuttonChenP4M6,'value',0);
set(handles.RadiobuttonDuan,'value',0);
function RadiobuttonChambers2007_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonNonedestriping,'value',0);
set(handles.RadiobuttonSwenson,'value',0);
set(handles.RadiobuttonChambers2007,'value',1);
set(handles.RadiobuttonChambers2012,'value',0);
set(handles.RadiobuttonChenP3M6,'value',0);
set(handles.RadiobuttonChenP4M6,'value',0);
set(handles.RadiobuttonDuan,'value',0);
function RadiobuttonChambers2012_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonNonedestriping,'value',0);
set(handles.RadiobuttonSwenson,'value',0);
set(handles.RadiobuttonChambers2007,'value',0);
set(handles.RadiobuttonChambers2012,'value',1);
set(handles.RadiobuttonChenP3M6,'value',0);
set(handles.RadiobuttonChenP4M6,'value',0);
set(handles.RadiobuttonDuan,'value',0);
function RadiobuttonChenP3M6_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonNonedestriping,'value',0);
set(handles.RadiobuttonSwenson,'value',0);
set(handles.RadiobuttonChambers2007,'value',0);
set(handles.RadiobuttonChambers2012,'value',0);
set(handles.RadiobuttonChenP3M6,'value',1);
set(handles.RadiobuttonChenP4M6,'value',0);
set(handles.RadiobuttonDuan,'value',0);
function RadiobuttonChenP4M6_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonNonedestriping,'value',0);
set(handles.RadiobuttonSwenson,'value',0);
set(handles.RadiobuttonChambers2007,'value',0);
set(handles.RadiobuttonChambers2012,'value',0);
set(handles.RadiobuttonChenP3M6,'value',0);
set(handles.RadiobuttonChenP4M6,'value',1);
set(handles.RadiobuttonDuan,'value',0);
function RadiobuttonDuan_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonNonedestriping,'value',0);
set(handles.RadiobuttonSwenson,'value',0);
set(handles.RadiobuttonChambers2007,'value',0);
set(handles.RadiobuttonChambers2012,'value',0);
set(handles.RadiobuttonChenP3M6,'value',0);
set(handles.RadiobuttonChenP4M6,'value',0);
set(handles.RadiobuttonDuan,'value',1);

% Specify the GIA removed or not
function RadiobuttonGIAnotRemoved_Callback(hObject, eventdata, handles)
    set(handles.RadiobuttonGIAnotRemoved, 'value', 1);
    set(handles.RadiobuttonGIARemovedGeru, 'value', 0);
function RadiobuttonGIARemovedGeru_Callback(hObject, eventdata, handles)
    set(handles.RadiobuttonGIAnotRemoved, 'value', 0);
    set(handles.RadiobuttonGIARemovedGeru, 'value', 1);

% Specify the Output Directory
function PushbuttonOutputDirectory_Callback(hObject, eventdata, handles)
OutputPathname= uigetdir('','Pick the output directory');
if OutputPathname==0
    %     set(handles.EditOutputDirectory, 'String', 'Output directry must be specified!!!');
    warndlg('Output directry must be specified!!!','Warning');
else
    set(handles.EditOutputDirectory, 'String', OutputPathname); % show directory
    handles.OutputPathname=[OutputPathname,'/'];
    guidata(hObject, handles);
end
function EditOutputDirectory_Callback(hObject, eventdata, handles)
function EditOutputDirectory_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Specify the Output Format
function EditOutputSHFilename_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonOutputFormatSH,'value',1);
set(handles.RadiobuttonOutputFormatGrid1degree,'value',0);
set(handles.RadiobuttonOutputFormatGrid025degree,'value',0);
function EditOutputSHFilename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EditOutputGridFilename1degree_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonOutputFormatSH,'value',0);
set(handles.RadiobuttonOutputFormatGrid1degree,'value',1);
set(handles.RadiobuttonOutputFormatGrid025degree,'value',0);
function EditOutputGridFilename1degree_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EditOutputGridFilename025degree_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonOutputFormatSH,'value',0);
set(handles.RadiobuttonOutputFormatGrid1degree,'value',0);
set(handles.RadiobuttonOutputFormatGrid025degree,'value',1);
function EditOutputGridFilename025degree_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RadiobuttonOutputFormatSH_Callback(hObject, eventdata, handles)
set(handles.EditOutputGridFilename1degree,'String','');
set(handles.EditOutputGridFilename025degree,'String','');
set(handles.RadiobuttonOutputFormatSH,'value',1);
set(handles.RadiobuttonOutputFormatGrid1degree,'value',0);
set(handles.RadiobuttonOutputFormatGrid025degree,'value',0);
function RadiobuttonOutputFormatGrid1degree_Callback(hObject, eventdata, handles)
set(handles.EditOutputSHFilename,'String','');
set(handles.EditOutputGridFilename025degree,'String','');
set(handles.RadiobuttonOutputFormatSH,'value',0);
set(handles.RadiobuttonOutputFormatGrid1degree,'value',1);
set(handles.RadiobuttonOutputFormatGrid025degree,'value',0);
function RadiobuttonOutputFormatGrid025degree_Callback(hObject, eventdata, handles)
set(handles.EditOutputSHFilename,'String','');
set(handles.EditOutputGridFilename1degree,'String','');
set(handles.RadiobuttonOutputFormatSH,'value',0);
set(handles.RadiobuttonOutputFormatGrid1degree,'value',0);
set(handles.RadiobuttonOutputFormatGrid025degree,'value',1);


% Save the Control File
function PushbuttonSaveControlFile_Callback(hObject, eventdata, handles)
[InputFilename,InputPathname,FilterIndex] = uiputfile('GMT_Control_File.txt','Save control file');
if FilterIndex==1 % file saved
    set(handles.EditSaveControlFile, 'String', [InputPathname,InputFilename]);
    set(handles.EditOpenControlFile, 'String', '');
    handles.PathControlFile=[InputPathname,InputFilename];
    guidata(hObject, handles);
    
    % Check the possible errors in Input Settings, Output Settings
    if isempty(get(handles.InputFileListbox,'String'))
        warndlg('Input files should be specified! Please select GRACE Level2 Files.','Warning');
        return;
    end
    if isempty(get(handles.EditMaxDegreeOutput,'String'))
        warndlg('Maximun Output Degree should be specified!','Warning');
        return;
    end
    if isnan(str2double(get(handles.EditMaxDegreeOutput,'String')))
        warndlg('Maximun Output Degree should be a number!','Warning');
        return;
    end
    if isempty(get(handles.EditFilterRadius,'String'))
        warndlg('Filter Radius should be specified!','Warning');
        return;
    end
    if isnan(str2double(get(handles.EditFilterRadius,'String')))
        warndlg('Gaussian Filter Radius should be a number!','Warning');
        return;
    end
    if isempty(get(handles.EditOutputDirectory,'String'))
        warndlg('Output Directory should be specified!','Warning');
        return;
    end
    if get(handles.RadiobuttonOutputFormatSH,'value') && isempty(get(handles.EditOutputSHFilename,'String'))
        warndlg('Name of SH output file should be specified!','Warning');
        return;
    end
    if get(handles.RadiobuttonOutputFormatGrid1degree,'value') && isempty(get(handles.EditOutputGridFilename1degree,'String'))
        warndlg('Name of Grid output file should be specified!','Warning');
        return;
    end
    if get(handles.RadiobuttonOutputFormatGrid025degree,'value') && isempty(get(handles.EditOutputGridFilename025degree,'String'))
        warndlg('Name of Grid output file should be specified!','Warning');
        return;
    end
    
    
    
    % Get the destriping method option
    if get(handles.RadiobuttonNonedestriping,'value')
        option_destriping='NONE';
    elseif get(handles.RadiobuttonSwenson,'value')
        option_destriping='SWENSON';
    elseif get(handles.RadiobuttonChambers2007,'value')
        option_destriping='CHAMBERS2007';
    elseif get(handles.RadiobuttonChambers2012,'value')
        option_destriping='CHAMBERS2012';
    elseif get(handles.RadiobuttonChenP3M6,'value')
        option_destriping='CHENP3M6';
    elseif get(handles.RadiobuttonChenP4M6,'value')
        option_destriping='CHENP4M6';
    elseif get(handles.RadiobuttonDuan,'value')
        option_destriping='DUAN';
    end
    
    
    % Get the GIA option
    if get(handles.RadiobuttonGIAnotRemoved,'value')
        option_gia='GIA_notRemoved';
    elseif get(handles.RadiobuttonGIARemovedGeru,'value')
        option_gia='GIA_Removed_Geru';
    end
    
    % Get the output format
    if get(handles.RadiobuttonOutputFormatSH,'value')
        OutputFileFormat=['SH_MAT ', get(handles.EditMaxDegreeOutput,'String'),' ', get(handles.EditOutputSHFilename,'String')];
    elseif get(handles.RadiobuttonOutputFormatGrid1degree,'value')
        OutputFileFormat=['GRID_MAT ',get(handles.EditMaxDegreeOutput,'String'),' ',get(handles.EditOutputGridFilename1degree,'String'),' 1'];
    elseif get(handles.RadiobuttonOutputFormatGrid025degree,'value')
        OutputFileFormat=['GRID_MAT ',get(handles.EditMaxDegreeOutput,'String'),' ',get(handles.EditOutputGridFilename025degree,'String'),' 0.25'];
    end

    % Get the file extension and determine the format of input files
    handles.InputFilename=get(handles.InputFileListbox,'String');
    if iscell(handles.InputFilename) % more than one input file
        [~, ~,FILE_TYPE]=fileparts([handles.InputPathname,handles.InputFilename{1}]);
    else % only one input file
        [~, ~,FILE_TYPE]=fileparts([handles.InputPathname,handles.InputFilename]);
    end
    if strcmp(FILE_TYPE,'.gfc') || strcmp(FILE_TYPE,'.GFC')
        InputFileFormat='ICGEM';
    else
        InputFileFormat='GRACE';
    end
  

    % Create Control File
    fid_write=fopen([InputPathname,InputFilename],'w');
    if iscell(handles.InputFilename) % more than one input file
        fprintf(fid_write,'%d\n',length(handles.InputFilename));
    else % only one input file
        fprintf(fid_write,'%d\n',1);
    end
    fprintf(fid_write,'%s\n',get(handles.EditFilterRadius,'String'));
    fprintf(fid_write,'%s\n',option_destriping);
    fprintf(fid_write,'%s\n',option_gia);
    fprintf(fid_write,'%s\n',InputFileFormat);
    fprintf(fid_write,'%s\n',OutputFileFormat);
    fprintf(fid_write,'%s\n',get(handles.EditC20, 'String'));
    fprintf(fid_write,'%s\n',get(handles.EditC21S21, 'String'));
    fprintf(fid_write,'%s\n',get(handles.EditC22S22, 'String'));
    fprintf(fid_write,'%s\n',get(handles.EditDegree1, 'String'));
    fprintf(fid_write,'%s\n',handles.InputPathname);
    fprintf(fid_write,'%s\n',strcat(get(handles.EditOutputDirectory,'String'),'/'));
    if iscell(handles.InputFilename) % more than one input file
        for jj=1:length(handles.InputFilename)-1
            fprintf(fid_write,'%s',handles.InputFilename{jj});
            fprintf(fid_write,'\n');
        end
        fprintf(fid_write,'%s',handles.InputFilename{length(handles.InputFilename)});
        fclose(fid_write);
    else % only one input file
        fprintf(fid_write,'%s',handles.InputFilename);
        fclose(fid_write);
    end
    
end
function EditSaveControlFile_Callback(hObject, eventdata, handles)
function EditSaveControlFile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Read the Control File
function PushbuttonOpenControlFile_Callback(hObject, eventdata, handles)
[InputFilename, InputPathname, filterindex] = uigetfile( ...
    {'*.*',  'All Files (*.*)'}, ...
    'Pick control file');
if filterindex==1 % file selected
    % Show the Directory and File in the textbox
    set(handles.EditOpenControlFile, 'String', [InputPathname,InputFilename]);
    set(handles.EditSaveControlFile, 'String', '');
    handles.PathControlFile=[InputPathname,InputFilename];
    guidata(hObject, handles);
    
    % Update the data
    % Read the Control File
    fid=fopen(handles.PathControlFile,'r');
    num_file        = fscanf(fid,'%d',1);
    radius_filter   = fscanf(fid,'%d',1);
    destrip_method  = fscanf(fid,'%s',1);
    option_gia      = fscanf(fid,'%s',1);
    InputFileFormat = fscanf(fid,'%s',1);
    OutputFileFormat_1= fscanf(fid,'%s',1);
    if strcmp(OutputFileFormat_1,'SH_MAT')
        OutputFileFormat_2= fscanf(fid,'%f',1); % max degree
        OutputFileFormat_3= fscanf(fid,'%s',1); % output name
    end
    if  strcmp(OutputFileFormat_1,'GRID_MAT')
        OutputFileFormat_2= fscanf(fid,'%f',1); % max degree
        OutputFileFormat_3= fscanf(fid,'%s',1); % output name
        OutputFileFormat_4= fscanf(fid,'%f',1); % spatial resolution
    end    
    dir_c20         = fscanf(fid,'%s',1);
    dir_c21_s21     = fscanf(fid,'%s',1);
    dir_c22_s22     = fscanf(fid,'%s',1);
    dir_degree_1    = fscanf(fid,'%s',1);
    dir_in          = fscanf(fid,'%s',1);
    dir_out         = fscanf(fid,'%s',1);
    file_name=cell(num_file,1);% Name list of input GRACE files
    for i = 1:num_file % read the name list of files in control file
        file_name{i} = fscanf(fid,'%s',1);
    end
    fclose(fid);
    handles.InputPathname=dir_in;
    % Show input files
    set(handles.InputFileListbox,'String',file_name);
    % Show C20 C21 S21 C22 S22 Degree 1 pathname
    if exist(dir_c20,'file') && ~strcmp(dir_c20,'NAN')
        set(handles.CheckboxC20,'Value',1);
        set(handles.EditC20,'String',dir_c20);
    end
    if exist(dir_c21_s21,'file') && ~strcmp(dir_c21_s21,'NAN')
        set(handles.CheckboxC21S21,'Value',1);
        set(handles.EditC21S21,'String',dir_c21_s21);
    end
    if exist(dir_c22_s22,'file') && ~strcmp(dir_c22_s22,'NAN')
        set(handles.CheckboxC22S22,'Value',1);
        set(handles.EditC22S22,'String',dir_c22_s22);
    end
    if exist(dir_degree_1,'file') && ~strcmp(dir_degree_1,'NAN')
        set(handles.CheckboxDegree1,'Value',1);
        set(handles.EditDegree1,'String',dir_degree_1);
    end
    % Show destriping option
    if strcmp(destrip_method,'NONE')
        set(handles.RadiobuttonNonedestriping,'Value',1);
        set(handles.RadiobuttonSwenson,'Value',0);
        set(handles.RadiobuttonChambers2007,'Value',0);
        set(handles.RadiobuttonChambers2012,'Value',0);
        set(handles.RadiobuttonChenP3M6,'Value',0);
        set(handles.RadiobuttonChenP4M6,'Value',0);
        set(handles.RadiobuttonDuan,'Value',0);
    elseif strcmp(destrip_method,'SWENSON')
        set(handles.RadiobuttonNonedestriping,'Value',0);
        set(handles.RadiobuttonSwenson,'Value',1);
        set(handles.RadiobuttonChambers2007,'Value',0);
        set(handles.RadiobuttonChambers2012,'Value',0);
        set(handles.RadiobuttonChenP3M6,'Value',0);
        set(handles.RadiobuttonChenP4M6,'Value',0);
        set(handles.RadiobuttonDuan,'Value',0);    
    elseif strcmp(destrip_method,'CHAMBERS2007')
        set(handles.RadiobuttonNonedestriping,'Value',0);
        set(handles.RadiobuttonSwenson,'Value',0);
        set(handles.RadiobuttonChambers2007,'Value',1);
        set(handles.RadiobuttonChambers2012,'Value',0);
        set(handles.RadiobuttonChenP3M6,'Value',0);
        set(handles.RadiobuttonChenP4M6,'Value',0);
        set(handles.RadiobuttonDuan,'Value',0); 
    elseif strcmp(destrip_method,'CHAMBERS2012')
        set(handles.RadiobuttonNonedestriping,'Value',0);
        set(handles.RadiobuttonSwenson,'Value',0);
        set(handles.RadiobuttonChambers2007,'Value',0);
        set(handles.RadiobuttonChambers2012,'Value',1);
        set(handles.RadiobuttonChenP3M6,'Value',0);
        set(handles.RadiobuttonChenP4M6,'Value',0);
        set(handles.RadiobuttonDuan,'Value',0);     
    elseif strcmp(destrip_method,'CHENP3M6')
        set(handles.RadiobuttonNonedestriping,'Value',0);
        set(handles.RadiobuttonSwenson,'Value',0);
        set(handles.RadiobuttonChambers2007,'Value',0);
        set(handles.RadiobuttonChambers2012,'Value',0);
        set(handles.RadiobuttonChenP3M6,'Value',1);
        set(handles.RadiobuttonChenP4M6,'Value',0);
        set(handles.RadiobuttonDuan,'Value',0); 
    elseif strcmp(destrip_method,'CHENP4M6')
        set(handles.RadiobuttonNonedestriping,'Value',0);
        set(handles.RadiobuttonSwenson,'Value',0);
        set(handles.RadiobuttonChambers2007,'Value',0);
        set(handles.RadiobuttonChambers2012,'Value',0);
        set(handles.RadiobuttonChenP3M6,'Value',0);
        set(handles.RadiobuttonChenP4M6,'Value',1);
        set(handles.RadiobuttonDuan,'Value',0); 
    elseif strcmp(destrip_method,'DUAN')
        set(handles.RadiobuttonNonedestriping,'Value',0);
        set(handles.RadiobuttonSwenson,'Value',0);
        set(handles.RadiobuttonChambers2007,'Value',0);
        set(handles.RadiobuttonChambers2012,'Value',0);
        set(handles.RadiobuttonChenP3M6,'Value',0);
        set(handles.RadiobuttonChenP4M6,'Value',0);
        set(handles.RadiobuttonDuan,'Value',1); 
    end
    % Show GIA option
    if strcmp(option_gia,'GIA_notRemoved')
        set(handles.RadiobuttonGIAnotRemoved,'Value',1);
        set(handles.RadiobuttonGIARemovedGeru,'Value',0);
    elseif strcmp(option_gia,'GIA_Removed_Geru')
        set(handles.RadiobuttonGIAnotRemoved,'Value',0);
        set(handles.RadiobuttonGIARemovedGeru,'Value',1);
    end
    % Show filter radius
    set(handles.EditFilterRadius,'string',radius_filter);
    % Show output directory
    if exist(dir_out,'dir')
        set(handles.EditOutputDirectory,'String',dir_out);
    end
    
    
    % Show output format
    if strcmp(OutputFileFormat_1,'SH_MAT')
        set(handles.RadiobuttonOutputFormatSH,'value',1);
        set(handles.RadiobuttonOutputFormatGrid1degree,'value',0);
        set(handles.RadiobuttonOutputFormatGrid025degree,'value',0);
        set(handles.EditMaxDegreeOutput,'String',OutputFileFormat_2);
        set(handles.EditOutputSHFilename,'String',OutputFileFormat_3);
    elseif strcmp(OutputFileFormat_1,'GRID_MAT') && OutputFileFormat_4==1
        set(handles.RadiobuttonOutputFormatSH,'value',0);
        set(handles.RadiobuttonOutputFormatGrid1degree,'value',1);
        set(handles.RadiobuttonOutputFormatGrid025degree,'value',0);
        set(handles.EditMaxDegreeOutput,'String',OutputFileFormat_2);
        set(handles.EditOutputGridFilename1degree,'String',OutputFileFormat_3);
    elseif strcmp(OutputFileFormat_1,'GRID_MAT') && OutputFileFormat_4==0.25
        set(handles.RadiobuttonOutputFormatSH,'value',0);
        set(handles.RadiobuttonOutputFormatGrid1degree,'value',0);
        set(handles.RadiobuttonOutputFormatGrid025degree,'value',1);
        set(handles.EditMaxDegreeOutput,'String',OutputFileFormat_2);
        set(handles.EditOutputGridFilename025degree,'String',OutputFileFormat_3);
    end
    
    
    guidata(hObject, handles);

end




function EditOpenControlFile_Callback(hObject, eventdata, handles)
function EditOpenControlFile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
