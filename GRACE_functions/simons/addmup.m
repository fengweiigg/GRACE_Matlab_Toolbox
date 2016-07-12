function nrmorl=addmup(nr,drk)
% nrmorl=ADDMUP(nr,drk)
%
% Calculates the number of real spherical harmonic orders that belong to
% an expansion from degree l=0 to L, OR vice versa.
% For arrays where m=0:l.
%
% INPUT:
%
% nr          The number of real spherical harmonic orders, OR
%             The degree of the expansion
% drk         'a' the input is the degree, calculate the number [default]
%             'r' the input is the number, calculate the degree
%
% OUTPUT:
%
% nrmorl      The degree of the expansion, OR
%             The number of real spherical harmonic orders
%
% NOTE:
%
% addmup(-1) equals 0, which is just as well
%
% See also ADDMOFF, ADDMOUT
%
% Last modified by fjsimons-at-alum.mit.edu, 05/17/2011

% Calculates the number of real spherical harmonic orders
% that belong to an expansion from degree l=0 to L - or vice versa. 
% Since m=0:l for real spherical harmonics, the number of orders per
% degree is l+1. This program thus calculates
% $\sum\limits_{l=0}^{L}(l+1)$ and its inverse, for a given L, using
% analytical expressions. Works for matrices, vectors as well as
% scalars. Note: degrees 0 and 1 are of course included.

% Check the behavior for the negatives should you ever need those

defval('nr',3)
defval('drk','a')

switch  drk
  case 'a'
    nrmorl=1/2*(nr+1).^2+1/2*nr+1/2;
    % nrmorl=1/2*nr.^2+3/2*nr+1;
  case 'r'
    nrmorl=-3/2+1/2*sqrt(1+8*nr);
    if prod((round(nrmorl)==nrmorl)+0)==0
      warning('Invalid entry - noninteger result')
    end
 otherwise
  error('Specify a valid option')
end
