function varargout = GRACE_Matlab_Toolbox_SHGrid(varargin)
% GRACE_MATLAB_TOOLBOX_SHGRID MATLAB code for GRACE_Matlab_Toolbox_SHGrid.fig
%      GRACE_MATLAB_TOOLBOX_SHGRID, by itself, creates a new GRACE_MATLAB_TOOLBOX_SHGRID or raises the existing
%      singleton*.
%
%      H = GRACE_MATLAB_TOOLBOX_SHGRID returns the handle to a new GRACE_MATLAB_TOOLBOX_SHGRID or the handle to
%      the existing singleton*.
%
%      GRACE_MATLAB_TOOLBOX_SHGRID('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRACE_MATLAB_TOOLBOX_SHGRID.M with the given input arguments.
%
%      GRACE_MATLAB_TOOLBOX_SHGRID('Property','Value',...) creates a new GRACE_MATLAB_TOOLBOX_SHGRID or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GRACE_Matlab_Toolbox_SHGrid_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GRACE_Matlab_Toolbox_SHGrid_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GRACE_Matlab_Toolbox_SHGrid

% Last Modified by GUIDE v2.5 24-Mar-2015 16:22:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GRACE_Matlab_Toolbox_SHGrid_OpeningFcn, ...
                   'gui_OutputFcn',  @GRACE_Matlab_Toolbox_SHGrid_OutputFcn, ...
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
end
% --- Executes just before GRACE_Matlab_Toolbox_SHGrid is made visible.
function GRACE_Matlab_Toolbox_SHGrid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GRACE_Matlab_Toolbox_SHGrid (see VARARGIN)

% Choose default command line output for GRACE_Matlab_Toolbox_SHGrid
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GRACE_Matlab_Toolbox_SHGrid wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = GRACE_Matlab_Toolbox_SHGrid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end


% --- Open SH coefficients file ---
function PushbuttonOpenSH_Callback(hObject, eventdata, handles)
[InputFilename, InputPathname, filterindex] = uigetfile( ...
    {'*.mat','Matlab-files (*.mat)'; ...
    '*.gfc','ICGEM-files (*.gfc)'}, ...
    'Pick SH coefficients file');
if filterindex==1 % file selected
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.gfc') || strcmp(FILE_TYPE,'.GFC')
        [cs,~]=gmt_readgfc(strcat(InputPathname,InputFilename));
    elseif strcmp(FILE_TYPE,'.mat')
        load(strcat(InputPathname,InputFilename));
        if exist('cs_grace','var')% there is 'cs_grace' variable in SH mat files
            handles.cs=cs_grace;
        elseif exist('cs','var')
            handles.cs=cs;
        else
            warndlg('There is no SH coefficients variable ''cs'' or ''cs_grace'' in .mat file!','Warning');
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
    else
        warndlg('Input file format is wrong!','Warning');
    end
    guidata(hObject,handles);
    set(handles.EditOpenSH,'String',strcat(InputPathname,InputFilename));
end
end

% --- Save Grid file ---
function PushbuttonSaveGrid_Callback(hObject, eventdata, handles)
[InputFilename,InputPathname,FilterIndex] = uiputfile('grid_data.mat','Save grid file');
if get(handles.Radiobutton1Degree,'value')
    type=1;
else
    type=0.25;
end
if FilterIndex==1 && isfield(handles,'cs')% create output file if there's input file
    set(handles.EditSaveGrid,'string',strcat(InputPathname,InputFilename));
    hh=msgbox('Sphercial Harmonic Synthesis is in processing.');
    pause(0.3);
    close(hh);
    grid_data=gmt_cs2grid(handles.cs,0,type);
    save(strcat(InputPathname, InputFilename),'grid_data');
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
    hh=msgbox('Spherical Harmonic Synthesis is done.');
    pause(0.5);
    close(hh);
end
end

% --- Open Grid file ---
function PushbuttonOpenGrid_Callback(hObject, eventdata, handles)
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
    set(handles.EditOpenGrid,'String',strcat(InputPathname,InputFilename));
end
end

% --- Save SH coefficents file ---
function PushbuttonSaveSH_Callback(hObject, eventdata, handles)

degree_max=str2double(get(handles.EditMaxDegree,'string'));
if isnan(degree_max) || degree_max~=fix(degree_max) % not integer
    warndlg('Input maximum degree is wrong!','Warning');
    return;
end

[InputFilename,InputPathname,FilterIndex] = uiputfile('cs.mat','Save SH coefficients file');

if FilterIndex==1 && isfield(handles,'grid_data') % create output file
    set(handles.EditSaveSH,'string',strcat(InputPathname,InputFilename));
    hh=msgbox('Spherical Harmonic Analysis is in processing.');
    pause(0.3);
    close(hh);
    cs=gmt_grid2cs(handles.grid_data,degree_max);
    save(strcat(InputPathname, InputFilename),'cs');
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
    hh=msgbox('Spherical Harmonic Analysis is done.');
    pause(0.5);
    close(hh);
end
end


function Radiobutton1Degree_Callback(hObject, eventdata, handles)
    set(handles.Radiobutton1Degree,'value',1);
    set(handles.Radiobutton025Degree,'value',0);

end
function Radiobutton025Degree_Callback(hObject, eventdata, handles)
    set(handles.Radiobutton1Degree,'value',0);
    set(handles.Radiobutton025Degree,'value',1);
end


function EditOpenSH_Callback(hObject, eventdata, handles)
end
function EditOpenSH_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function EditSaveGrid_Callback(hObject, eventdata, handles)
end
function EditSaveGrid_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function EditOpenGrid_Callback(hObject, eventdata, handles)

end
function EditOpenGrid_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function EditSaveSH_Callback(hObject, eventdata, handles)
end
function EditSaveSH_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end




function EditMaxDegree_Callback(hObject, eventdata, handles)

end
function EditMaxDegree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditMaxDegree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
