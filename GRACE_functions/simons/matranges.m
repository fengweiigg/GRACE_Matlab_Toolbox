function indexmat=matranges(ranges)
% indexmat=MATRANGES(ranges)
%
% Makes an index vector with monotonically increasing indices between
% pairs of numbers supplied as input.
%
% See also ADDMOUT, ADDMON. This is a very fast routine.
%
% EXAMPLE:
%
% matranges([1 4 1 2 -1 2])
%
% Last modified by fjsimons-at-alum.mit.edu, 16.09.2005

ranges=ranges(:)';
if mod(length(ranges),2)
  error('Ranges must form pairs')
end

lower=ranges(1:2:end);
upper=ranges(2:2:end);
hulp1=ones(1,sum(upper-lower+1));
hulp2=cumsum(upper(1:end-1)-lower(1:end-1)+1)+1;
hulp1(hulp2)=diff([upper(1:end-1) ; lower(2:end)],1);
hulp1(1)=ranges(1);

indexmat=cumsum(hulp1);

