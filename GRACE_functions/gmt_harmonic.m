function [ Amplitude1, Amplitude1_std, Phase1,Phase1_std, Amplitude2, Amplitude2_std, Phase2, Phase2_std, Trend, Trend_std, Trend_line, Resid, Interp] = gmt_harmonic(t,t1,grid_data,grid_data_std)

% Estimate annual, semi-annual and linear trend of time series
% 
% INPUT:
%   t              time epoch
%   t1             the time epoch for the predicted (missing) months
%   grid_data      the input field grid_data(lat,lon,epochs), 2-D or 3-D matrix
%   grid_data_std  the standard deviations of input field
%
% OUTPUT:
%   Amplitude1     amplitude of estimated annual cycle
%   Amplitude1_std 2*standard deviation of estimated annual amplitude
%   Phase1         phase of estimated cycle, angle of cosine
%   Phase1_std     2*standard deviation of estimated phase
%   Amplitude2     amplitude of estimated semi-annual cycle
%   Amplitude2_std 2*standard deviation of estimated semi-annual amplitude
%   Phase2         phase of estimated cycle, angle of cosine
%   Phase2_std     2*standard deviation of estimated phase
%   Trend          trend of time series
%   Trend_std      2*standard deviation of estimated trend
%   Trend_line     predicted trend line
%   Interp         predicted trend + seasonal cycles
%   Resid          residuals after removing trend and seasonal cycles
%
% If standard deviations of input field are not provided, the uncertainties 
% (Amplitude1_std,Phase1_std,Amplitude2_std,Phase2_std, Trend_std) have been estimated 
% as two standard deviations after propagation of monthly value errors in the least squares fit procedure, 
% which represent the 95% confidence interval.
% If standard deviations of input field are not provided, the uncertainties
% here represent formal errors in the least square fit procedure.
%
% FENG Wei 25/03/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

if ndims(grid_data)==2  % grid is 2-D matrix
    rows=1;
    cols=1;
    epochs= max(size(grid_data));
    grid_data=reshape(grid_data,[rows,cols,epochs]);
    if nargin<3 || nargin>4
        error('Number of input parameters in gmt_harmonic function is wrong!');
    end
    if nargin==3
        grid_data_std = [];
    else
        grid_data_std=reshape(grid_data_std,[rows,cols,epochs]);
    end
end

[rows,cols,epochs] = size(grid_data);

if nargin<3 || nargin>4
    error('Number of input parameters in gmt_harmonic function is wrong!');
end

if nargin==3
    grid_data_std = [];
end

omega = 2*pi;

% ------------------
% initial
% ------------------
time_series     = zeros(1,epochs);
Trend           = zeros(rows, cols);
Trend_std       = zeros(rows, cols);
Amplitude1      = zeros(rows, cols);
Amplitude1_std  = zeros(rows, cols);
Phase1          = zeros(rows, cols);
Phase1_std      = zeros(rows, cols);
Amplitude2      = zeros(rows, cols);
Amplitude2_std  = zeros(rows, cols);
Phase2          = zeros(rows, cols);
Phase2_std      = zeros(rows, cols);
epochs2=max(size(t1));
Trend_line      = zeros(rows, cols, epochs2);
Resid           = zeros(rows, cols, epochs);
Interp          = zeros(rows, cols, epochs2);
% ------------------
% end
% ------------------


% create normal matrix
A = [];
% ones function create unit matrix with epochs X 1
A(:,1) = (ones(epochs,1));
A(:,2) = t';
A(:,3) = cos(omega*t'); % annual
A(:,4) = sin(omega*t');
A(:,5) = cos(omega*2*t'); % semi-annual
A(:,6) = sin(omega*2*t');
NM = (A'*A)\A';

% No standard deviations of input field
if nargin==3
    for i = 1:rows
        for j = 1:cols
            for k = 1:epochs
                time_series(k) = grid_data(i,j,k);
            end
            if sum(isnan(time_series))==0 % there is no NaN in time series
                % normal equation, x contains all estimated parameters
                x    = NM*time_series';
                % cofactor matrix
                Exx = (A'*A)\eye(size(A'*A));
                %Exx = inv((A'*A));
                
                %
                % annual ampltidue & phase
                %
                Ampl = sqrt(x(3)^2 + x(4)^2);
                if x(4)>0  %  sin >0 [0, pi/2], [pi/2, pi]
                    Pha = acos(x(3)/Ampl);
                elseif x(4)<0 % sin<0 [pi, 3*pi/2], [3*pi/2, 2*pi]
                    Pha = 2*pi - acos(x(3)/Ampl);
                elseif x(4)==0 && x(3)==1
                    Pha = 0;
                elseif x(4)==0 && x(3)==-1
                    Pha = pi;
                elseif x(4)==0 && x(3)==0
                    Pha = 0;
                end
                % from radian measure to angle measure
                Pha  = Pha*180/pi;
                % from angle measure to day in the year
                % Pha  = Pha*365/360;
                Amplitude1(i,j) = Ampl; % annual amplitude
                Phase1(i,j)     = Pha; % annual phase
                
                %
                % semi-annual ampltidue & phase
                %
                Ampl = sqrt(x(5)^2 + x(6)^2);
                if x(6)>0  %  sin >0 [0, pi/2], [pi/2, pi]
                    Pha = acos(x(5)/Ampl);
                elseif x(6)<0 % sin<0 [pi, 3*pi/2], [3*pi/2, 2*pi]
                    Pha = 2*pi - acos(x(5)/Ampl);
                elseif x(6)==0 && x(5)==1
                    Pha = 0;
                elseif x(6)==0 && x(5)==-1
                    Pha = pi;
                elseif x(6)==0 && x(5)==0
                    Pha = 0;
                end
                % from radian measure to angle measure
                Pha  = Pha*180/pi;
                % from angle measure to day in the year
                % Pha  = Pha*365/360/2;
                Amplitude2(i,j) = Ampl; % semi-annual amplitude
                Phase2(i,j)      = Pha; % semi-annual phase
                
                % predict time series at epochs t1 using bias, trend, seasonal
                % cycles
                if ~isempty(t1)
                    Trend_line(i,j,:) = x(1) + x(2)*t1';
                    Oscilation_int=  x(1) + x(2)*t1' + x(3)*cos(omega*t1') + x(4)*sin(omega*t1') + x(5)*cos(2*omega*t1') + x(6)*sin(2*omega*t1');
                    Interp(i,j,:)=Oscilation_int;
                end
                
                % residual
                res = time_series' - A*x;
                Resid(i,j,:) = res;
                
                % aposteriori variance estimate: Error^2/(number of obs - number of param)
                % aposteriori unit weight mean error 验后单位权中误差
                var_est = res'*res/(length(time_series) - 6);
                % aposteriori covariance matrix for estimated parameters from aposteriori
                % unit weight mean error and co-factor matrix
                % 协因数阵乘以单位权中误差得到验后的观测值权中误差矩阵
                Exx_scal = var_est*Exx;
                
                Trend(i,j)     = x(2);
                Trend_std(i,j) = 2*sqrt(Exx_scal(2,2));
                
                % 1 sigma: 68.3% ; 2 sigma: 95.4%; 3 sigma: 99.7%
                Amplitude1_std(i,j) = 2*sqrt( ( x(3)*x(3) * Exx_scal(3,3) + x(4)*x(4) * Exx_scal(4,4) )/(Amplitude1(i,j)*Amplitude1(i,j)));
                if abs(x(4)/x(3))<=1
                    x_43 = x(4)/x(3);
                    Phase1_std(i,j) = (1-x_43^2+x_43^4-x_43^6+x_43^8-x_43^10)^2 * ( (x(4)/x(3)^2)^2 * Exx_scal(3,3) + (1/x(3))^2 * Exx_scal(4,4)) ;
                    Phase1_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                else
                    x_34 = x(3)/x(4);
                    Phase1_std(i,j) =( x_34^2-x_34^4+x_34^6-x_34^8+x_34^10)^2 * ( (x(4)/x(3)^2)^2 * Exx_scal(3,3) + (1/x(3))^2 * Exx_scal(4,4)) ;
                    Phase1_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                end
                
                Amplitude2_std(i,j) = 2*sqrt( ( x(5)*x(5) * Exx_scal(5,5) + x(6) * x(6) * Exx_scal(6,6) )/(Amplitude2(i,j)*Amplitude2(i,j)));
                if abs(x(6)/x(5))<=1
                    x_65 = x(6)/x(5);
                    Phase2_std(i,j) = (1-x_65^2+x_65^4-x_65^6+x_65^8-x_65^10)^2 * ( (x(6)/x(5)^2)^2 * Exx_scal(5,5) + (1/x(5))^2 * Exx_scal(6,6)) ;
                    Phase2_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                else
                    x_56 = x(5)/x(6);
                    Phase2_std(i,j) =( x_56^2-x_56^4+x_56^6-x_56^8+x_56^10)^2 * ( (x(6)/x(5)^2)^2 * Exx_scal(5,5) + (1/x(5))^2 * Exx_scal(6,6)) ;
                    Phase2_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                end
            else % there is NaN in time series
                Amplitude1(i,j) = NaN; % annual amplitude
                Amplitude1_std(i,j)=NaN;
                Phase1(i,j)     = NaN; % annual phase
                Phase1_std(i,j) = NaN;
                Amplitude2(i,j) = NaN; % semi-annual amplitude
                Phase2(i,j)      = NaN; % semi-annual phase
                Phase2_std(i,j) = NaN;
                Trend(i,j)     = NaN;
                Trend_std(i,j) = NaN;
                Trend_line(i,j,1:epochs)=NaN;
                Resid(i,j,1:epochs) =NaN;
                Interp(i,j,1:epochs)=NaN;
            end
        end
    end
    
    % There are standard deviations of input field
elseif nargin==4
    for i = 1:rows
        for j = 1:cols
            for k = 1:epochs
                time_series(k) = grid_data(i,j,k);
                time_series_var(k) = grid_data_std(i,j,k)^2;
            end
            if sum(isnan(time_series))==0 % there is no NaN in time series
                % set a priori unit weight mean error is ONE
                % weights
                P = diag(1./time_series_var);
                
                % normal equation
                % x    = (A'*P*A)\eye(size(A'*P*A))*A'*P*time_series';
                x    = (A'*P*A)\A'*P*time_series';
                
                %  cofactor matrix
                %Exx = (A'*P*A)\eye(size(A'*P*A));
                Exx = inv((A'*P*A));
                
                %
                % annual ampltidue & phase
                %
                Ampl = sqrt(x(3)^2 + x(4)^2);
                if x(4)>0  %  sin >0 [0, pi/2], [pi/2, pi]
                    Pha = acos(x(3)/Ampl);
                elseif x(4)<0 % sin<0 [pi, 3*pi/2], [3*pi/2, 2*pi]
                    Pha = 2*pi - acos(x(3)/Ampl);
                elseif x(4)==0 && x(3)==1
                    Pha = 0;
                elseif x(4)==0 && x(3)==-1
                    Pha = pi;
                elseif x(4)==0 && x(3)==0
                    Pha = 0;
                end
                % from radian measure to angle measure
                Pha  = Pha*180/pi;
                % from angle measure to day in the year
                % Pha  = Pha*365/360;
                Amplitude1(i,j) = Ampl;
                Phase1(i,j)      = Pha;
                
                %
                % semi-annual  ampltidue & phase
                %
                Ampl = sqrt(x(5)^2 + x(6)^2);
                if x(6)>0  %  sin >0 [0, pi/2], [pi/2, pi]
                    Pha = acos(x(5)/Ampl);
                elseif x(6)<0 % sin<0 [pi, 3*pi/2], [3*pi/2, 2*pi]
                    Pha = 2*pi - acos(x(5)/Ampl);
                elseif x(6)==0 && x(5)==1
                    Pha = 0;
                elseif x(6)==0 && x(5)==-1
                    Pha = pi;
                elseif x(6)==0 && x(5)==0
                    Pha = 0;
                end
                % from radian measure to angle measure
                Pha  = Pha*180/pi;
                % from angle measure to day in the year
                % Pha  = Pha*365/360/2;
                Amplitude2(i,j) = Ampl;
                Phase2(i,j)      = Pha;
                
                
                % predict time series at epochs t1 using bias, trend, seasonal
                % cycles
                if ~isempty(t1)
                    Trend_line(i,j,:) = x(1) + x(2)*t1';
                    Oscilation_int=  x(1) + x(2)*t1' + x(3)*cos(omega*t1') + x(4)*sin(omega*t1') + x(5)*cos(2*omega*t1') + x(6)*sin(2*omega*t1');
                    Interp(i,j,:)=Oscilation_int;
                end
                
                % residual
                res = time_series' - A*x;
                Resid(i,j,:) = res;
                
                % covariance matrix for estimated parameters
                Exx_scal = Exx; % 用验前的单位权中误差
                
                Trend(i,j)     = x(2);
                Trend_std(i,j) = 2*sqrt(Exx_scal(2,2));
                
                % 1 sigma: 68.3% ; 2 sigma: 95.4%; 3 sigma: 99.7%
                Amplitude1_std(i,j) = 2*sqrt( ( x(3)*x(3) * Exx_scal(3,3) + x(4)*x(4) * Exx_scal(4,4) )/(Amplitude1(i,j)*Amplitude1(i,j)));
                if abs(x(4)/x(3))<=1
                    x_43 = x(4)/x(3);
                    Phase1_std(i,j) = (1-x_43^2+x_43^4-x_43^6+x_43^8-x_43^10)^2 * ( (x(4)/x(3)^2)^2 * Exx_scal(3,3) + (1/x(3))^2 * Exx_scal(4,4)) ;
                    Phase1_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                else
                    x_34 = x(3)/x(4);
                    Phase1_std(i,j) =( x_34^2-x_34^4+x_34^6-x_34^8+x_34^10)^2 * ( (x(4)/x(3)^2)^2 * Exx_scal(3,3) + (1/x(3))^2 * Exx_scal(4,4)) ;
                    Phase1_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                end
                
                Amplitude2_std(i,j) = 2*sqrt( ( x(5)*x(5) * Exx_scal(5,5) + x(6) * x(6) * Exx_scal(6,6) )/(Amplitude2(i,j)*Amplitude2(i,j)));
                if abs(x(6)/x(5))<=1
                    x_65 = x(6)/x(5);
                    Phase2_std(i,j) = (1-x_65^2+x_65^4-x_65^6+x_65^8-x_65^10)^2 * ( (x(6)/x(5)^2)^2 * Exx_scal(5,5) + (1/x(5))^2 * Exx_scal(6,6)) ;
                    Phase2_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                else
                    x_56 = x(5)/x(6);
                    Phase2_std(i,j) =( x_56^2-x_56^4+x_56^6-x_56^8+x_56^10)^2 * ( (x(6)/x(5)^2)^2 * Exx_scal(5,5) + (1/x(5))^2 * Exx_scal(6,6)) ;
                    Phase2_std(i,j) = 2*sqrt(Phase1_std(i,j))*180/pi;
                end
            else % there is NaN in time series
                Amplitude1(i,j) = NaN; % annual amplitude
                Amplitude1_std(i,j)=NaN;
                Phase1(i,j)     = NaN; % annual phase
                Phase1_std(i,j) = NaN;
                Amplitude2(i,j) = NaN; % semi-annual amplitude
                Phase2(i,j)      = NaN; % semi-annual phase
                Phase2_std(i,j) = NaN;
                Trend(i,j)     = NaN;
                Trend_std(i,j) = NaN;
                Trend_line(i,j,1:epochs)=NaN;
                Resid(i,j,1:epochs) =NaN;
                Interp(i,j,1:epochs)=NaN;
            end
        end
    end
end
    
