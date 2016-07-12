function varargout = GRACE_Matlab_Toolbox_HarmonicAnalysis(varargin)
% GRACE_MATLAB_TOOLBOX_HARMONICANALYSIS MATLAB code for GRACE_Matlab_Toolbox_HarmonicAnalysis.fig
%      GRACE_MATLAB_TOOLBOX_HARMONICANALYSIS, by itself, creates a new GRACE_MATLAB_TOOLBOX_HARMONICANALYSIS or raises the existing
%      singleton*.
%
%      H = GRACE_MATLAB_TOOLBOX_HARMONICANALYSIS returns the handle to a new GRACE_MATLAB_TOOLBOX_HARMONICANALYSIS or the handle to
%      the existing singleton*.
%
%      GRACE_MATLAB_TOOLBOX_HARMONICANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRACE_MATLAB_TOOLBOX_HARMONICANALYSIS.M with the given input arguments.
%
%      GRACE_MATLAB_TOOLBOX_HARMONICANALYSIS('Property','Value',...) creates a new GRACE_MATLAB_TOOLBOX_HARMONICANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GRACE_Matlab_Toolbox_HarmonicAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GRACE_Matlab_Toolbox_HarmonicAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GRACE_Matlab_Toolbox_HarmonicAnalysis

% Last Modified by GUIDE v2.5 25-Mar-2015 21:59:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GRACE_Matlab_Toolbox_HarmonicAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @GRACE_Matlab_Toolbox_HarmonicAnalysis_OutputFcn, ...
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


% --- Executes just before GRACE_Matlab_Toolbox_HarmonicAnalysis is made visible.
function GRACE_Matlab_Toolbox_HarmonicAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GRACE_Matlab_Toolbox_HarmonicAnalysis (see VARARGIN)

% Choose default command line output for GRACE_Matlab_Toolbox_HarmonicAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GRACE_Matlab_Toolbox_HarmonicAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GRACE_Matlab_Toolbox_HarmonicAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EditOpenSeries_Callback(hObject, eventdata, handles)
% hObject    handle to EditOpenSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditOpenSeries as text
%        str2double(get(hObject,'String')) returns contents of EditOpenSeries as a double


% --- Executes during object creation, after setting all properties.
function EditOpenSeries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditOpenSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PushbuttonOpenSeries.
function PushbuttonOpenSeries_Callback(hObject, eventdata, handles)
[InputFilename, InputPathname, filterindex] = uigetfile( ...
    {'*.mat','Matlab-files (*.mat)'}, ...
    'Pick series file');
if filterindex==1 % file selected
    set(handles.EditOpenSeries,'string',strcat(InputPathname,InputFilename));
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.mat')
        load(strcat(InputPathname,InputFilename));
        if exist('time_series_grace','var')
            time_series=time_series_grace;
        end
        if ~exist('time_series','var')
           warndlg('There is no ''time_series'' or ''time_series_grace'' variable in .mat file!','Warning');
           return;
        end

        if exist('time_grace','var')
            time=time_grace;
        end
        if ~exist('time','var')
           warndlg('There is no ''time'' or ''time_grace'' variable in .mat file!','Warning');
           return;
        end
        
        grid_data(1,1,:)=time_series;
        
        [ Amplitude1, Amplitude1_std, Phase1,Phase1_std, Amplitude2,...
            Amplitude2_std, Phase2, Phase2_std, Trend, Trend_std,...
            ~, ~, ~] = gmt_harmonic(time,[],grid_data);

        set(handles.EditAnnualAmplitude,'string',strcat(num2str(Amplitude1,'%3.2e'),'+/-',num2str(Amplitude1_std,'%3.2e')));
        set(handles.EditAnnualPhase,'string',strcat(num2str(Phase1,'%4.1f'),'+/-',num2str(Phase1_std,'%4.1f')));
        set(handles.EditSemiannualAmplitude,'string',strcat(num2str(Amplitude2,'%3.2e'),'+/-',num2str(Amplitude2_std,'%3.2e')));
        set(handles.EditSemiannualPhase,'string',strcat(num2str(Phase2,'%4.1f'),'+/-',num2str(Phase2_std,'%4.1f')));
        set(handles.EditTrend,'string',strcat(num2str(Trend,'%3.2e'),'+/-',num2str(Trend_std,'%3.2e')));        
    end
end
        

function EditAnnualAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to EditAnnualAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditAnnualAmplitude as text
%        str2double(get(hObject,'String')) returns contents of EditAnnualAmplitude as a double


% --- Executes during object creation, after setting all properties.
function EditAnnualAmplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditAnnualAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSemiannualAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to EditSemiannualAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSemiannualAmplitude as text
%        str2double(get(hObject,'String')) returns contents of EditSemiannualAmplitude as a double


% --- Executes during object creation, after setting all properties.
function EditSemiannualAmplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSemiannualAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditTrend_Callback(hObject, eventdata, handles)
% hObject    handle to EditTrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditTrend as text
%        str2double(get(hObject,'String')) returns contents of EditTrend as a double


% --- Executes during object creation, after setting all properties.
function EditTrend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditTrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditAnnualPhase_Callback(hObject, eventdata, handles)
% hObject    handle to EditAnnualPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditAnnualPhase as text
%        str2double(get(hObject,'String')) returns contents of EditAnnualPhase as a double


% --- Executes during object creation, after setting all properties.
function EditAnnualPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditAnnualPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSemiannualPhase_Callback(hObject, eventdata, handles)
% hObject    handle to EditSemiannualPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSemiannualPhase as text
%        str2double(get(hObject,'String')) returns contents of EditSemiannualPhase as a double


% --- Executes during object creation, after setting all properties.
function EditSemiannualPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSemiannualPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSaveSemiannualphase_Callback(hObject, eventdata, handles)
% hObject    handle to EditSaveSemiannualphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSaveSemiannualphase as text
%        str2double(get(hObject,'String')) returns contents of EditSaveSemiannualphase as a double


% --- Executes during object creation, after setting all properties.
function EditSaveSemiannualphase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSaveSemiannualphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSaveAnnualphase_Callback(hObject, eventdata, handles)
% hObject    handle to EditSaveAnnualphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSaveAnnualphase as text
%        str2double(get(hObject,'String')) returns contents of EditSaveAnnualphase as a double


% --- Executes during object creation, after setting all properties.
function EditSaveAnnualphase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSaveAnnualphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSaveTrend_Callback(hObject, eventdata, handles)
% hObject    handle to EditSaveTrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSaveTrend as text
%        str2double(get(hObject,'String')) returns contents of EditSaveTrend as a double


% --- Executes during object creation, after setting all properties.
function EditSaveTrend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSaveTrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSaveSemiannualamplitude_Callback(hObject, eventdata, handles)
% hObject    handle to EditSaveSemiannualamplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSaveSemiannualamplitude as text
%        str2double(get(hObject,'String')) returns contents of EditSaveSemiannualamplitude as a double


% --- Executes during object creation, after setting all properties.
function EditSaveSemiannualamplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSaveSemiannualamplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditSaveAnnualamplitude_Callback(hObject, eventdata, handles)
% hObject    handle to EditSaveAnnualamplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditSaveAnnualamplitude as text
%        str2double(get(hObject,'String')) returns contents of EditSaveAnnualamplitude as a double


% --- Executes during object creation, after setting all properties.
function EditSaveAnnualamplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSaveAnnualamplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PushbuttonOpenGrid.
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
            return;
        end
        if exist('time','var')
            handles.time=time;
        elseif exist('time_grace','var')
            handles.time=time_grace;
        else
            warndlg('There is no ''time'' or ''time_grace'' variable in .mat file!','Warning');
            return;
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

function EditOpenGrid_Callback(hObject, eventdata, handles)
% hObject    handle to EditOpenGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditOpenGrid as text
%        str2double(get(hObject,'String')) returns contents of EditOpenGrid as a double


% --- Executes during object creation, after setting all properties.
function EditOpenGrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditOpenGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Save Annual amplitude ---
function PushbuttonSaveAnnualamplitude_Callback(hObject, eventdata, handles)
[InputFilename,InputPathname,filterindex] = uiputfile('annual_amplitude.mat','Save annual amplitude file');
if filterindex==1 % file selected
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.mat')
        set(handles.EditSaveAnnualamplitude,'String',strcat(InputPathname,InputFilename));
    end
end
% --- Save Annual phase ---
function PushbuttonSaveAnnualphase_Callback(hObject, eventdata, handles)
[InputFilename,InputPathname,filterindex] = uiputfile('annual_phase.mat','Save annual phase file');
if filterindex==1 % file selected
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.mat')
        set(handles.EditSaveAnnualphase,'String',strcat(InputPathname,InputFilename));
    end
end

% --- Save Semiannual amplitude ---
function PushbuttonSaveSemiannualamplitude_Callback(hObject, eventdata, handles)
[InputFilename,InputPathname,filterindex] = uiputfile('semiannual_amplitude.mat','Save semi-annual amplitude file');
if filterindex==1 % file selected
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.mat')
        set(handles.EditSaveSemiannualamplitude,'String',strcat(InputPathname,InputFilename));
    end
end

% --- Save Semiannual phase ---
function PushbuttonSaveSemiannualphase_Callback(hObject, eventdata, handles)
[InputFilename,InputPathname,filterindex] = uiputfile('semiannual_phase.mat','Save semi-annual phase file');
if filterindex==1 % file selected
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.mat')
        set(handles.EditSaveSemiannualphase,'String',strcat(InputPathname,InputFilename));
    end
end

% --- Save Trend ---
function PushbuttonSaveTrend_Callback(hObject, eventdata, handles)
[InputFilename,InputPathname,filterindex] = uiputfile('trend.mat','Save trend file');
if filterindex==1 % file selected
    [~,~,FILE_TYPE]=fileparts(strcat(InputPathname,InputFilename));
    if strcmp(FILE_TYPE,'.mat')
        set(handles.EditSaveTrend,'String',strcat(InputPathname,InputFilename));
    end
end

% --- Executes on button press in PushbuttonCalculate.
function PushbuttonCalculate_Callback(hObject, eventdata, handles)

    hh=msgbox('Spatial Harmonic Analysis is in processing.');
    pause(0.3);
    close(hh);
[ Amplitude1, Amplitude1_std, Phase1,Phase1_std, Amplitude2,...
    Amplitude2_std, Phase2, Phase2_std, Trend, Trend_std,...
    ~, ~, ~] = gmt_harmonic(handles.time,[],handles.grid_data);

if ~isempty(get(handles.EditSaveAnnualamplitude,'string'))
    annual_amplitude=Amplitude1;
    save(get(handles.EditSaveAnnualamplitude,'string'),'annual_amplitude');
end

if ~isempty(get(handles.EditSaveAnnualphase,'string'))
    annual_phase=Phase1;
    save(get(handles.EditSaveAnnualphase,'string'),'annual_phase');
end

if ~isempty(get(handles.EditSaveSemiannualamplitude,'string'))
    semiannual_amplitude=Amplitude2;
    save(get(handles.EditSaveSemiannualamplitude,'string'),'semiannual_amplitude');
end

if ~isempty(get(handles.EditSaveSemiannualphase,'string'))
    semiannual_phase=Phase2;
    save(get(handles.EditSaveSemiannualphase,'string'),'semiannual_phase');
end

if ~isempty(get(handles.EditSaveTrend,'string'))
    trend=Trend;
    save(get(handles.EditSaveTrend,'string'),'trend');
end

hh=msgbox('Spatial harmonic analysis is done.');
pause(1);
close(hh);
