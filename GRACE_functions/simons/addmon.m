function [dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm,demin]=addmon(L,m)
% [dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm,demin]=ADDMON(L,m)
% Returns spherical harmonic degree and positive-order indexing arrays.
%
% INPUT:
%
% L         Maximal degree of the expansion
% m         An order requested [default: NaN, nothing requested]
%
% OUTPUT:
%
% dems      Vector of orders  m involved in the expansion, [0 01 012] 
% dels      Vector of degrees l involved in the expansion, [0 11 222]
% mz        Index to the m=0 locations for the zonal coefficients
% lmcosi    Matrix with zeroes where cos/sine coefficients will go
% mzi       Insertion location of m=0 sin terms in vector not listing them
% mzo       Running index into lmcosi(3:4) returning vector without
%           the sin coefficient at m=0, effectively unwrapping the
%           cosine and sine elements and reverting back to +/-m in the
%           ordering of orders [0 0-11 0-11-22] as in, e.g. KERNELC
% bigm      Single vector with ordered orders, [0 0-11 0-11-22] as in KERNELC
% bigl      Single vector with ordered degrees, [0 111 22222 3333333] etc
% rinm      Vector reindexing bigm back to [0 -101 -2-1012 -3-2-10123]
%           etc, as in GLMALPHA|PTO
% ronm      Running index into lmcosi(3:4) unwrapping the orders as in ADDMOUT
% demin     The index into vector dems listing all entries that belong to
%           the order m given at the input if a specific request is made
%
% SEE ALSO:
%
% XYZ2PLM, PLM2XYZ, MATRANGES, ADDMOUT, ADDMIN, ADDMUP
%
% EXAMPLE:
%
% [M,L,mz,BM]=addmout(10); [a,b,c,d,e,f,M2,L2,rinm]=addmon(10);
% difer(M-M2(rinm))
% difer(L-L2(rinm)) 
% difer(M(BM)-M2(rinm(BM)))
% difer(L(BM)-L2(rinm(BM)))
%
% EXAMPLE: 
%
% Find the positions of the order m in an lmcosi matrix up to degree L
%
% L=15; m=round(rand*L); 
% [dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm,demin]=addmon(L,m);
% dems(demin)
%
% Last modified by fjsimons-at-alum.mit.edu, 1/20/2011

defval('m',NaN)

matr=(repmat(0:L,2,1)'*diag([0 1]))';
dems=matranges(matr(:)')';
els=0:L;
dels=gamini(els,els+1)';

if nargout>=3
  mz=[1 cumsum(els(1:end-1)+1)+1]';
  if nargout>=4
    lmcosi=[dels dems zeros(length(dels),2)];
    if nargout>=5
      mzi=(mz*2+1)-[1:length(mz)]';
      if nargout>=6
	% Just number the [dels dems] vector sequentially down the rows
	indo=reshape(1:length(dems)*2,length(dems),2);
	% And stick a zero in where the sine coefficient goes
	indo(length(dels)+mz)=0; indo=indo'; 
	% Simply pick the nonzero elements on this index array
	mzo=indo(~~indo);
	if L>1
	  % Note that this never skips more than two positions
	  difer(unique(diff(sort(mzo)))-[1 2 3]',[],[],NaN)
	end
	if nargout>=7
	  fm=[-dems dems]; 
	  bigm=fm(mzo);
	  if nargout>=8
	    fl=[dels dels];
	    bigl=fl(mzo);
	    if nargout>=9
	      if  L>=2
		% Must do the first couple of degrees myself since gamini(x,0) is bad
		lopos=cumsum(2*els(1:end-1)+1)+1; 
		ka=[repmat(1,1,L-1) ; els(3:end)-1 ; repmat(1,1,L-1); els(3:end)];
		ko=[[1 0 -1 2] gamini(repmat([0 -2 -1 2],1,L-1),ka(:)')]';
		ko(lopos)=[2:2:2*L];
		rinm=cumsum(ko);  
	      elseif L==0
		rinm=1;
	      elseif L==1
		rinm=[1 3 2 4]';
	      end
	      if nargout>=10
		ronm=mzo(rinm);
		if nargout>=11
		  demin=addmup(m-1:L-1)+m+1;
		end
	      end
	    end
	  end
	end
      end
    end
  end
end
