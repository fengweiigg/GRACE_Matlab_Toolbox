function varargout = GRACE_Matlab_Toolbox(varargin)
% GRACE_MATLAB_TOOLBOX MATLAB code for GRACE_Matlab_Toolbox.fig
%      GRACE_MATLAB_TOOLBOX, by itself, creates a new GRACE_MATLAB_TOOLBOX or raises the existing
%      singleton*.
%
%      H = GRACE_MATLAB_TOOLBOX returns the handle to a new GRACE_MATLAB_TOOLBOX or the handle to
%      the existing singleton*.
%
%      GRACE_MATLAB_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRACE_MATLAB_TOOLBOX.M with the given input arguments.
%
%      GRACE_MATLAB_TOOLBOX('Property','Value',...) creates a new GRACE_MATLAB_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GRACE_Matlab_Toolbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GRACE_Matlab_Toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GRACE_Matlab_Toolbox

% Last Modified by GUIDE v2.5 04-Sep-2015 23:20:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GRACE_Matlab_Toolbox_OpeningFcn, ...
                   'gui_OutputFcn',  @GRACE_Matlab_Toolbox_OutputFcn, ...
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


% --- Executes just before GRACE_Matlab_Toolbox is made visible.
function GRACE_Matlab_Toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GRACE_Matlab_Toolbox (see VARARGIN)

% Choose default command line output for GRACE_Matlab_Toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GRACE_Matlab_Toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GRACE_Matlab_Toolbox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% --------------------------------------------------------------------
function GSMGADProcessing_Callback(hObject, eventdata, handles)
% open('GRACE_Matlab_Toolbox_preprocessing.fig');  %result??????????
run GRACE_Matlab_Toolbox_preprocessing
% hObject    handle to GSMGADProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function GSMGADProcessing_ClickedCallback(hObject, eventdata, handles)
% open('GRACE_Matlab_Toolbox_preprocessing.fig');  %result??????????
run GRACE_Matlab_Toolbox_preprocessing


% --------------------------------------------------------------------
function LeakageReductionSpatial_Callback(hObject, eventdata, handles)
% hObject    handle to LeakageReductionSpatial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run GRACE_Matlab_Toolbox_LeakageReductionSpatial

% --------------------------------------------------------------------
function SHGrid_Callback(hObject, eventdata, handles)
% hObject    handle to SHGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% open('GRACE_Matlab_Toolbox_SHGrid.fig');  %result??????????
run GRACE_Matlab_Toolbox_SHGrid

% --------------------------------------------------------------------
function Grid2Series_Callback(hObject, eventdata, handles)
% hObject    handle to Grid2Series (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run GRACE_Matlab_Toolbox_Grid2Series
% --------------------------------------------------------------------
function HarmonicAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to HarmonicAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% open('GRACE_Matlab_Toolbox_HarmonicAnalysis.fig');  %result??????????

run GRACE_Matlab_Toolbox_HarmonicAnalysis


% --------------------------------------------------------------------
function GRACEProcessing_Callback(hObject, eventdata, handles)
% hObject    handle to GRACEProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function DataAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to DataAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
