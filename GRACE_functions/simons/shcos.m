function [C,b,e]=shcos(lmcosi,l)
% [C,b,e]=SHCOS(lmcosi,l)
%
% Finds COSINE coefficients of degree l from a matrix (or vector) listing
% (l m) Ccos (Csin)  (sorted but can start anywhere)
% Ccos               (sorted and starting at zero)
%
% OUTPUT
% 
% C     The relevant COSINE expansion coefficients
% b     Row number of the first one in lmcosi
% e     Row number of the last one in lmcosi
%
% See also SHSIN
%
% Last modified by fjsimons-at-alum.mit.edu, Feb 16th 2004

% Vector or matrix?
if size(lmcosi,2)==1
  lmin=0;
else
  lmin=min(lmcosi(:,1));
end

% Get the indices
b=addmup(l-1)+1-addmup(lmin-1);
e=addmup(l)-addmup(lmin-1);

% Get the coefficients
if size(lmcosi,2)==1
  C=lmcosi(b:e,1);
else
  C=lmcosi(b:e,3);
end

