function W = gmt_gaussian(max_degree,radius_filter)

% Compute weights for every degree in Gaussian smoothing
% 
% INPUT:
%   max_degree      [1 x 1]  maximum degree
%   radius_filter   [1 x 1]  Radius of Gaussian smoothing, unit: km 
%                            (i.e., half width of Gaussian weighting function)
%
% OUTPUT:
%   W    [max_degree+1 x 1]  weighting coefficients
%
% The gaussian weighting function is calculated based on equation (34) of Chambers et al.
% JGR, 2006, Observing seasonal steric sea level variations with GRACE and
% satellite altimetry
%
% FENG Wei 18/12/2015
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

% Input checking
% -------------------------
if length(max_degree)   ~= 1, error('Degree must be scalar.'); end
if max_degree < 2, error('Maximum degree must be higher than 2.'); end
if length(radius_filter) ~= 1, error('Cap size must be scalar.'); end

W = zeros(max_degree+1,1);

for i=1:max_degree+1
    W(i)=exp(-((i-1)*radius_filter/6371)^2/(4*log(2)));
end

% % Method2: Recursive algorithm according to Wahr et. al. (1998) equation (34) is
% % not recommended since its instabilities beyond degree 50.
%
% b = log(2)/(1-cos(radius_filter/6371));
% % recursive calculation of the weighting coefficients
% W(1) = 1;
% W(2) =  (1+exp(-2*b))/(1-exp(-2*b)) - 1/b ;
% for i = 2:max_degree
%     W(i+1) = W(i-1) - (2*(i-1)+1)/b * W(i);
% %     W(i+1) = W(i-1) - (2*i+1)/b * W(i); % Wahr et. al. (1998) equation (34) is wrong 
%     if W(i+1) > W(i) || W(i+1) < 0, W(i+1) = 0; end
% end

    

