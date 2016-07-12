function varargout = GRACE_Matlab_Toolbox_LeakageReductionSpatial(varargin)
% GRACE_MATLAB_TOOLBOX_LEAKAGEREDUCTIONSPATIAL MATLAB code for GRACE_Matlab_Toolbox_LeakageReductionSpatial.fig
%      GRACE_MATLAB_TOOLBOX_LEAKAGEREDUCTIONSPATIAL, by itself, creates a new GRACE_MATLAB_TOOLBOX_LEAKAGEREDUCTIONSPATIAL or raises the existing
%      singleton*.
%
%      H = GRACE_MATLAB_TOOLBOX_LEAKAGEREDUCTIONSPATIAL returns the handle to a new GRACE_MATLAB_TOOLBOX_LEAKAGEREDUCTIONSPATIAL or the handle to
%      the existing singleton*.
%
%      GRACE_MATLAB_TOOLBOX_LEAKAGEREDUCTIONSPATIAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRACE_MATLAB_TOOLBOX_LEAKAGEREDUCTIONSPATIAL.M with the given input arguments.
%
%      GRACE_MATLAB_TOOLBOX_LEAKAGEREDUCTIONSPATIAL('Property','Value',...) creates a new GRACE_MATLAB_TOOLBOX_LEAKAGEREDUCTIONSPATIAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GRACE_Matlab_Toolbox_LeakageReductionSpatial_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GRACE_Matlab_Toolbox_LeakageReductionSpatial_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GRACE_Matlab_Toolbox_LeakageReductionSpatial

% Last Modified by GUIDE v2.5 05-Sep-2015 17:37:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GRACE_Matlab_Toolbox_LeakageReductionSpatial_OpeningFcn, ...
                   'gui_OutputFcn',  @GRACE_Matlab_Toolbox_LeakageReductionSpatial_OutputFcn, ...
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


% --- Executes just before GRACE_Matlab_Toolbox_LeakageReductionSpatial is made visible.
function GRACE_Matlab_Toolbox_LeakageReductionSpatial_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GRACE_Matlab_Toolbox_LeakageReductionSpatial (see VARARGIN)

% Choose default command line output for GRACE_Matlab_Toolbox_LeakageReductionSpatial
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GRACE_Matlab_Toolbox_LeakageReductionSpatial wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GRACE_Matlab_Toolbox_LeakageReductionSpatial_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EditOpenSH_Callback(hObject, eventdata, handles)
% hObject    handle to EditOpenSH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditOpenSH as text
%        str2double(get(hObject,'String')) returns contents of EditOpenSH as a double


% --- Executes during object creation, after setting all properties.
function EditOpenSH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditOpenSH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PushbuttonOpenSH.
function PushbuttonOpenSH_Callback(hObject, eventdata, handles)
% hObject    handle to PushbuttonOpenSH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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


function EditSaveSH_Callback(hObject, eventdata, handles)
% hObject    handle to EditSaveSH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSaveSH as text
%        str2double(get(hObject,'String')) returns contents of EditSaveSH as a double


% --- Executes during object creation, after setting all properties.
function EditSaveSH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSaveSH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PushbuttonSaveSH.
function PushbuttonSaveSH_Callback(hObject, eventdata, handles)

[InputFilename,InputPathname,FilterIndex] = uiputfile('cs_leakage_removed.mat','Save SH coefficients file');

if FilterIndex==1 && isfield(handles,'cs') % create output file
    set(handles.EditSaveSH,'string',strcat(InputPathname,InputFilename));
    hh=msgbox('Leakage reduction in spatial domain is in processing.');
    pause(0.3);
    close(hh);
    % Processing leakage reduction
    %     cs=gmt_grid2cs(handles.grid_data,degree_max);
    
    % Get leakage reduction method
    if get(handles.RadiobuttonOceanLeakageReduction,'value')
        type='ocean';
    else
        type='land';
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
    radius_filter=get(handles.EditFilterRadius,'String');
    cs=gmt_cs2leakagefreecs(handles.cs,type,option_destriping,radius_filter);

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
    hh=msgbox('Leakage reduction in spatial domain is done.');
    pause(0.5);
    close(hh);
end



% --- Executes on button press in RadiobuttonOceanLeakageReduction.
function RadiobuttonOceanLeakageReduction_Callback(hObject, eventdata, handles)
    set(handles.RadiobuttonOceanLeakageReduction,'value',1);
    set(handles.RadiobuttonLandLeakageReduction,'value',0);


% --- Executes on button press in RadiobuttonLandLeakageReduction.
function RadiobuttonLandLeakageReduction_Callback(hObject, eventdata, handles)
    set(handles.RadiobuttonOceanLeakageReduction,'value',0);
    set(handles.RadiobuttonLandLeakageReduction,'value',1);


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



function EditFilterRadius_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

function EditFilterRadius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
