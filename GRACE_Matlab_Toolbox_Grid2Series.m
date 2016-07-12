function varargout = GRACE_Matlab_Toolbox_Grid2Series(varargin)
% GRACE_MATLAB_TOOLBOX_GRID2SERIES MATLAB code for GRACE_Matlab_Toolbox_Grid2Series.fig
%      GRACE_MATLAB_TOOLBOX_GRID2SERIES, by itself, creates a new GRACE_MATLAB_TOOLBOX_GRID2SERIES or raises the existing
%      singleton*.
%
%      H = GRACE_MATLAB_TOOLBOX_GRID2SERIES returns the handle to a new GRACE_MATLAB_TOOLBOX_GRID2SERIES or the handle to
%      the existing singleton*.
%
%      GRACE_MATLAB_TOOLBOX_GRID2SERIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRACE_MATLAB_TOOLBOX_GRID2SERIES.M with the given input arguments.
%
%      GRACE_MATLAB_TOOLBOX_GRID2SERIES('Property','Value',...) creates a new GRACE_MATLAB_TOOLBOX_GRID2SERIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GRACE_Matlab_Toolbox_Grid2Series_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GRACE_Matlab_Toolbox_Grid2Series_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GRACE_Matlab_Toolbox_Grid2Series

% Last Modified by GUIDE v2.5 24-Mar-2015 18:05:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GRACE_Matlab_Toolbox_Grid2Series_OpeningFcn, ...
                   'gui_OutputFcn',  @GRACE_Matlab_Toolbox_Grid2Series_OutputFcn, ...
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


% --- Executes just before GRACE_Matlab_Toolbox_Grid2Series is made visible.
function GRACE_Matlab_Toolbox_Grid2Series_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GRACE_Matlab_Toolbox_Grid2Series (see VARARGIN)

% Choose default command line output for GRACE_Matlab_Toolbox_Grid2Series
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GRACE_Matlab_Toolbox_Grid2Series wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GRACE_Matlab_Toolbox_Grid2Series_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EditInputGrid_Callback(hObject, eventdata, handles)
% hObject    handle to EditInputGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditInputGrid as text
%        str2double(get(hObject,'String')) returns contents of EditInputGrid as a double


% --- Executes during object creation, after setting all properties.
function EditInputGrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditInputGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Open Grid file ---
function PushbuttonInputGrid_Callback(hObject, eventdata, handles)
[InputFilename, InputPathname, filterindex] = uigetfile( ...
    {'*.mat','Matlab-files (*.mat)'}, ...
    'Pick grid file');
if filterindex==1 % file selected
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.mat')
        load(strcat(InputPathname,InputFilename));
        if exist('grid_data_grace','var')% there is 'cs_grace' variable in SH mat files 
            handles.grid_data=grid_data_grace;
        elseif exist('grid_data','var') % no 'cs' variable in SH mat files
            handles.grid_data=grid_data;
        else    
            warndlg('There is no ''grid_data'' or ''grid_data_grace'' variable in .mat file!','Warning');
        end
        if exist('time','var')
            handles.time=time;
        end
        if exist('time_grace','var')
            handles.time=time_grace;
        end
        if exist('str_year','var')
            handles.str_year=str_year;
        end
        if exist('str_month','var')
            handles.str_month=str_month;
        end
    end
    guidata(hObject,handles);
    set(handles.EditInputGrid,'String',strcat(InputPathname,InputFilename));
end

% --- Open boundary file ---
function PushbuttonBoundaryFile_Callback(hObject, eventdata, handles)
[InputFilename, InputPathname, filterindex] = uigetfile( ...
    {'*.bln','Boundary-files (*.bln)'}, ...
    'Pick boundary file');
if filterindex==1 % file selected
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.bln')
        handles.boundaryfile=strcat(InputPathname,InputFilename);
    end
    guidata(hObject,handles);
    set(handles.EditBoundaryFile,'String',strcat(InputPathname,InputFilename));
end
    
% --- Open mask file ---
function PushbuttonMaskFile_Callback(hObject, eventdata, handles)
[InputFilename, InputPathname, filterindex] = uigetfile( ...
    {'*.*','All files (*.*)'}, ...
    'Pick mask file');
if filterindex==1 % file selected
%     [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
%     if strcmp(FILE_TYPE,'.bln')
        handles.boundaryfile=strcat(InputPathname,InputFilename);
%     end
    guidata(hObject,handles);
    set(handles.EditMaskFile,'String',strcat(InputPathname,InputFilename));
end

function RadiobuttonBoundaryfile_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonBoundaryfile,'value',1);
set(handles.RadiobuttonMaskfile,'value',0);
function RadiobuttonMaskfile_Callback(hObject, eventdata, handles)
set(handles.RadiobuttonBoundaryfile,'value',0);
set(handles.RadiobuttonMaskfile,'value',1);

function EditBoundaryFile_Callback(hObject, eventdata, handles)
function EditBoundaryFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditBoundaryFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EditMaskFile_Callback(hObject, eventdata, handles)
function EditMaskFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditMaskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Save series file ---
function PushbuttonOutputSeries_Callback(hObject, eventdata, handles)

[InputFilename,InputPathname,filterindex] = uiputfile('time_series.mat','Save time series file');

if filterindex==1 % file selected
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.mat')
        set(handles.EditOutputSeries,'String',strcat(InputPathname,InputFilename));
        hh=msgbox('Grid2Series is in processing.');
        pause(0.3);
        close(hh);
        %     end
        if get(handles.RadiobuttonBoundaryfile,'value')
            time_series=gmt_grid2series(handles.grid_data,handles.boundaryfile,'line');
        elseif get(handles.RadiobuttonMaskfile,'value')
            time_series=gmt_grid2series(handles.grid_data,handles.maskfile,'mask');
        end
        save(strcat(InputPathname, InputFilename),'time_series');
        % add other variables if they are in input SH file *.mat
        if isfield(handles,'str_year') % whether str_year exists in handles
            str_year= handles.str_year;
            save(strcat(InputPathname, InputFilename),'-append','str_year');
        end
        if isfield(handles,'str_month')
            str_month= handles.str_month;
            save(strcat(InputPathname, InputFilename),'-append','str_month');
        end
        if isfield(handles,'time')
            time= handles.time;
            save(strcat(InputPathname, InputFilename),'-append','time');
        end
        hh=msgbox('Grid2Series is done.');
        pause(0.5);
        close(hh);
    end
end

function EditOutputSeries_Callback(hObject, eventdata, handles)
function EditOutputSeries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditOutputSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PushbuttonOutputSeries.
% hObject    handle to PushbuttonOutputSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
