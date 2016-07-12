function varargout=plm2xyz(lmcosi,degres,c11cmn,lmax,latmax,Plm)
% [r,lon,lat,Plm]=PLM2XYZ(lmcosi,degres,c11cmn,lmax,latmax,Plm)
% [r,lon,lat,Plm]=PLM2XYZ(lmcosi,lat,lon,lmax,latmax,Plm)
%
% Inverse (4*pi-normalized real) spherical harmonic transform.
%
% Compute a spatial field from spherical harmonic coefficients given as 
% [l m Ccos Csin] (not necessarily starting from zero, but sorted), with
% degree resolution 'degres' [default: approximate Nyquist degree].
%
% INPUT:
%
% lmcosi     Matrix listing l,m,cosine and sine expansion coefficients
%            e.g. those coming out of ADDMON
% degres     Longitude/ latitude spacing, in degrees [default: Nyquist] OR
%            "lat": a column vector with latitudes [degrees]
% c11cmn     Corner nodes of lon/lat grid [default: 0 90 360 -90] OR
%            "lon": a column vector with longitudes [degrees]
% lmax       Maximum bandwidth expanded at a time [default: 720]
% latmax     Maximum linear size of the latitude grid allowed [default: Inf]
% Plm        The appropriate Legendre polynomials should you already have them
% 
% OUTPUT:
%
% r          The field (matrix for a grid, vector for scattered points)
% lon,lat    The grid (matrix) or evaluation points (vector), in degrees
% Plm        The set of appropriate Legendre polynomials should you want them
% 
% EXAMPLE:
%
% plm2xyz('demo1') % Illustrates forward and inverse transform
% plm2xyz('demo2',fra) % with 'fra' a data fraction
% plm2xyz('demo3') % Plots EGM96 versus EGM2008 free-air anomaly
% plm2xyz('demo4') % Plots GTM3AR versus EGM2008 topography
% plm2xyz('demo5') % Plots EGM96 geopotential
% plm2xyz('demo6') % EGM2008 topography somewhere
% plm2xyz('demo7') % Tests the expansion to scattered points
% plm2xyz('demo8') % EGM2008 gravity globally
%
% NOTE:
%
% The degree range is now split intelligently in blocks of degrees whose
% memory requirements never exceed the initial expansion from 0 to lmax.
% For very high bandwidth models, specify a small region, and this
% together with 'degres' will determine memory requirements. Built-in
% maxima using 'latmax' and 'lmax'.
%
% NOTE:
%
% Compare to YLM which has a different (-1)^m phase factor.
%
% See also XYZ2PLM, PLM2SPEC, TH2PL, PL2TH, YLM
%
% Special thanks to kwlewis-at-princeton.edu for spotting a bug.
% Last modified by fjsimons-at-alum.mit.edu, 07/13/2012
if ~isstr(lmcosi)
  t0=clock;
  % Lowest degree of the expansion
  lmin=lmcosi(1);
  % Highest degree (bandwidth of the expansion)
  L=lmcosi(end,1);
  % Never use Libbrecht algorithm... found out it wasn't that good
  defval('libb',0)
  % Default resolution is the Nyquist degree; return equal sampling in
  % longitude and latitude; sqrt(L*(L+1)) is equivalent wavelength
  degN=180/sqrt(L*(L+1));
  defval('degres',degN);
  
  % When do you get a task bar?
  taskmax=100;
  taskinf=0;
  
  % Default grid is all of the planet
  defval('c11cmn',[0 90 360 -90]);
  
  % Build in maxima to save memory
  defval('latmax',Inf); 
  defval('lmax',720);

  % But is it a grid or are they merely scattered points?
  if length(degres)==length(c11cmn)
    % It's a bunch of points!
    nlat=length(degres);
    nlon=length(c11cmn);
    % Colatitude vector in radians
    theta=[90-degres(:)']*pi/180;
    % Longitude vector in radians
    phi=c11cmn(:)'*pi/180;

    % Initialize output vector
    r=repmat(0,nlat,1);
    
    % Now if this is too large reduce lmax, our only recourse to hardcode
    ntb=256;
    if round(sqrt(nlat)) >= ntb || round(sqrt(nlon)) >= ntb
      lmax=round(ntb);
    end
  elseif length(degres)==1 && length(c11cmn)==4
    % It's a grid
    if degres>degN
      disp('PLM2XYZ: You can do better! Ask for more spatial resolution')
      disp(sprintf('Spatial sampling ALLOWED: %8.3f ; REQUESTED: %6.3f',...
		   degN,degres))
    end
    % The number of longitude and latitude grid points that will be computed
    nlon=min(ceil([c11cmn(3)-c11cmn(1)]/degres+1),latmax);
    nlat=min(ceil([c11cmn(2)-c11cmn(4)]/degres+1),2*latmax+1);

    % Initialize output grid
    r=repmat(0,nlat,nlon);
    
    % Longitude grid vector in radians
    phi=linspace(c11cmn(1)*pi/180,c11cmn(3)*pi/180,nlon);
    % Colatitude grid vector in radians
    theta=linspace([90-c11cmn(2)]*pi/180,[90-c11cmn(4)]*pi/180,nlat);
    
%    disp(sprintf('Creating %i by %i grid with resolution %8.3f',nlat,nlon,degres))
  else
    error('Make up your mind - is it a grid or a list of points?')
  end

  % Here we were going to build an option for a polar grid
  % But abandon this for now
  % [thetap,phip]=rottp(theta,phi,pi/2,pi/2,0);

  % Piecemeal degree ranges
  % Divide the degree range increments spaced such that the additional
  % number of degrees does not exceed addmup(lmax)
  % If this takes a long time, abort it
  els=0; ind=0;
  while els<L
     ind=ind+1;
    % Take positive root
    els(ind+1)=min(floor(max(roots(...
	[1 3 -els(ind)^2-3*els(ind)-2*addmup(lmax)]))),L);
    if any(diff(els)==0)
      error('Increase lmax as you are not making progress')
    end
  end
  % Now els contains the breakpoints of the degrees
  if ~all(diff(addmup(els))<=addmup(lmax))
    error('The subdivision of the degree scale went awry')
  end

  % Here's the lspacings
  if length(els)>2
    els=pauli(els,2)+...
	[0 0 ; ones(length(els)-2,1) zeros(length(els)-2,1)];
  end

  for ldeg=1:size(els,1)
    ldown=els(ldeg,1);
    lup=els(ldeg,2);
    % Construct the filename
    fnpl=sprintf('%s/LSSM-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'LEGENDRE'),...
		 ldown,lup,nlat);
    % ONLY COMPLETE LINEARLY SPACED SAMPLED VECTORS ARE TO BE SAVED!
    if [exist(fnpl,'file')==2 && length(c11cmn)==4 && all(c11cmn==[0 90 360 -90])]...
	  & ~[size(els,1)==1 &&  exist('Plm','var')==1]
      % Get Legendre function values at linearly spaced intervals
      disp(sprintf('Using preloaded %s',fnpl))
      load(fnpl)
      % AND TYPICALLY ANYTHING ELSE WOULD BE PRECOMPUTED, BUT THE GLOBAL
      % ONES CAN TOO! The Matlabpool check doesn't seem to work inside 
    elseif size(els,1)==1 &&  exist('Plm','var')==1 && matlabpool('size')==0
      % disp(sprintf('Using precomputed workspace Legendre functions'))
    else
      % Evaluate Legendre polynomials at selected points
      try
        Plm=nan(length(theta),addmup(lup)-addmup(ldown-1));
      catch ME
        error(sprintf('\n %s \n\n Decrease lmax in PLM2XYZ \n',ME.message))
      end
      if [lup-ldown]>taskmax && length(degres) > taskinf
	h=waitbar(0,sprintf(...
	    'Evaluating Legendre polynomials between %i and %i',...
	    ldown,lup));
      end
      in1=0;
      in2=ldown+1;
      % Always start from the beginning in this array, regardless of lmin
      for l=ldown:lup
	if libb==0
	  Plm(:,in1+1:in2)=(legendre(l,cos(theta(:)'),'sch')*sqrt(2*l+1))';
	else
	  Plm(:,in1+1:in2)=(libbrecht(l,cos(theta(:)'),'sch')*sqrt(2*l+1))';
	end
	in1=in2;
	in2=in1+l+2;
	if [lup-ldown]>taskmax && length(degres) > taskinf
	  waitbar((l-ldown+1)/(lup-ldown+1),h)
	end
      end
      if [lup-ldown]>taskmax && length(degres) > taskinf
	delete(h)
      end
      if length(c11cmn)==4 && all(c11cmn==[0 90 360 -90])
	save(fnpl,'Plm','-v7.3')
	disp(sprintf('Saved %s',fnpl))
      end
    end

    % Loop over the degrees
    more off
%    disp(sprintf('PLM2XYZ Expansion from %i to %i',max(lmin,ldown),lup))
    for l=max(lmin,ldown):lup
      % Compute Schmidt-normalized Legendre functions at 
      % the cosine of the colatitude (=sin(lat)) and 
      % renormalize them to the area of the unit sphere
      
      % Remember the Plm vector always starts from ldown
      b=addmup(l-1)+1-addmup(ldown-1);
      e=addmup(l)-addmup(ldown-1);
      
      plm=Plm(:,b:e)';
      
      m=0:l;
      mphi=m(:)*phi(:)';
      
      % Normalization of the harmonics is to 4\ pi, the area of the unit
      % sphere: $\int_{0}^{\pi}\int_{0}^{2\pi}
      % |P_l^m(\cos\theta)\cos(m\phi)|^2\sin\theta\,d\theta d\phi=4\pi$.
      % Note the |cos(m\phi)|^2 d\phi contributes exactly \pi for m~=0
      % and 2\pi for m==0 which explains the absence of sqrt(2) there;
      % that fully normalized Legendre polynomials integrate to 1/2/pi
      % that regular Legendre polynomials integrate to 2/(2l+1),
      % Schmidt polynomials to 4/(2l+1) for m>0 and 2/(2l+1) for m==0,
      % and Schmidt*sqrt(2l+1) to 4 or 2. Note that for the integration of
      % the harmonics you get either 2pi for m==0 or pi+pi for cosine and
      % sine squared (cross terms drop out). This makes the fully
      % normalized spherical harmonics the only ones that consistently give
      % 1 for the normalization of the spherical harmonics.
      % Note this is using the cosines only; the "spherical harmonics" are
      % actually only semi-normalized.
      % Test normalization as follows (using inaccurate Simpson's rule):
      defval('tst',0)
      if tst
	f=(plm.*plm)'.*repmat(sin(theta(:)),1,l+1);
	c=(cos(mphi).*cos(mphi))';
	ntest=simpson(theta,f).*simpson(phi,c)/4/pi;
	disp(sprintf('Mean Normalization Error l= %3.3i: %8.3e',...
		     l,(sum(abs(1-ntest)))/(l+1)))
	% For a decent test you would use "legendreprodint"
      end

      % Find the cosine and sine coefficients for this degree
      clm=shcos(lmcosi,l);
      slm=shsin(lmcosi,l);
      
      fac1=repmat(clm,1,nlon).*cos(mphi)+...
	   repmat(slm,1,nlon).*sin(mphi);
      
      % Sum over all orders and (through loop) over all degrees
      if length(degres)==length(c11cmn)
	expa=sum(plm.*fac1,1)'; 
	% Or diag(plm'*fac1) if you will
      elseif length(degres)==1 & length(c11cmn)==4
	expa=plm'*fac1;
      end
      r=r+expa;
    end
  end
    
  lon=phi*180/pi;
  lat=90-theta*180/pi;

  % Prepare output
  vars={r,lon,lat,Plm};
  varargout=vars(1:nargout);
  
%  disp(sprintf('PLM2XYZ (Synthesis) took %8.4f s',etime(clock,t0)))
  
elseif strcmp(lmcosi,'demo1')
  lmax=30; [m,l,mzero]=addmon(lmax);
  c=randn(addmup(lmax),2).*([l l].^(-1)); 
  c(1)=3; c(mzero,2)=0; lmcosi=[l m c];
  L=30;
  [r,lon,lat]=plm2xyz(lmcosi,180/sqrt(L*(L+1)));
  C1=xyz2plm(r,L,'simpson');
  C2=xyz2plm(r,L,'gl');
  C3=xyz2plm(r,L,'im');
  lat=linspace(90,-90,size(r,1));
  C4=xyz2plm(r,L,'im',lat);
  C5=xyz2plm(r,L,'irr');
  clf
  ah(1)=subplot(211);
  p1(1)=plot(abs(lmcosi(:,3)-C1(1:addmup(lmax),3)),'b+-'); hold on
  p1(2)=plot(abs(lmcosi(:,3)-C2(1:addmup(lmax),3)),'rv-');
  p1(3)=plot(abs(lmcosi(:,3)-C3(1:addmup(lmax),3)),'kx-');
  p1(4)=plot(abs(lmcosi(:,3)-C4(1:addmup(lmax),3)),'go-'); 
  p1(5)=plot(abs(lmcosi(:,3)-C5(1:addmup(lmax),3)),'mo-'); hold off 
  ylim([-0.1 1]*1e-14)
  yl(1)=ylabel('Absolute error');
  xl(1)=xlabel('Cumulative degree and order');
  legend('simpson','gl','im','imlat','irr')

  ah(2)=subplot(212);
  p2(1)=plot(abs(lmcosi(:,3)-C1(1:addmup(lmax),3)),'b+-'); hold on
  p2(2)=plot(abs(lmcosi(:,3)-C2(1:addmup(lmax),3)),'rv-');
  p2(3)=plot(abs(lmcosi(:,3)-C3(1:addmup(lmax),3)),'kx-');
  p2(4)=plot(abs(lmcosi(:,3)-C4(1:addmup(lmax),3)),'go-'); 
  p2(5)=plot(abs(lmcosi(:,3)-C5(1:addmup(lmax),3)),'mo-'); hold off 
  axis tight
  yl(2)=ylabel('Absolute error');
  xl(2)=xlabel('Cumulative degree and order');

  longticks(ah,2)
  set([p1 p2],'MarkerS',4)
  set([xl yl],'FontS',15)
  fig2print(gcf,'landscape')
  figdisp
elseif strcmp(lmcosi,'demo2')
  lmax=10; L=10;
  [m,l,mzero]=addmon(lmax);
  c=randn(addmup(lmax),2).*([l l].^(-1)); 
  c(1)=3; c(mzero,2)=0; lmcosi=[l m c];
  [r,lon,lat]=plm2xyz(lmcosi,180/sqrt(L*(L+1)));
  tol=length(lon)*length(lat);
  defval('degres',0.4)
  fra=degres;
  unform=2;
  [LON,LAT]=meshgrid(lon,lat);
  if unform==1
    % Not really uniform
    indo=indeks(shuffle([1:tol]'),1:ceil(fra*tol));
    lonr=LON(indo);
    latr=LAT(indo);
  else 
    % Really uniform on the sphere
    [lonr,latr]=randsphere(ceil(fra*tol));
    indo=sub2ind(size(r),ceil(scale(latr,[1 length(lat)])),...
		 ceil(scale(lonr,[1 length(lon)])));
    lonr=LON(indo);
    latr=LAT(indo);
  end
  C6=xyz2plm(r(indo),L,'irr',latr,lonr);
  clf
  ah(1)=subplot(221);
  plotplm(r,lon*pi/180,lat*pi/180,1)
  xl(1)=title('Input');
  hold on
  [x,y]=mollweide(lonr*pi/180,latr*pi/180);
  p1=plot(x,y,'o');
  ck=caxis;
  cb(1)=colorbar('hor');
  ah(2)=subplot(222);
  rec=plm2xyz(C6);
  plotplm(rec,lon*pi/180,lat*pi/180,1)
  hold on
  p4=plot(x,y,'o');
  xl(2)=title('Reconstruction');
  caxis(ck)
  cb(2)=colorbar('hor');
  ah(3)=subplot(223);
  ps(1)=plot([0:L],plm2spec(lmcosi),'bo-');
  hold on
  ps(2)=plot([0:L],plm2spec(C6),'rv-');
  set(ps(1),'MarkerE','b','MarkerF','r','LineW',2,'MarkerS',4)
  set(ps(2),'MarkerE','r','MarkerF','r','LineW',2,'MarkerS',4)
  set(gca,'Yscale','log'); grid
  xl(3)=xlabel('Degree');
  yl(1)=ylabel('Power');
  xl(5)=title(sprintf('Spectral comparison, %s=%3.1f','\alpha',fra));
  lg=legend('Input','Recovered');
  ah(4)=subplot(224);
  difo=abs(r-rec);
  difo(difo<1e-12)=NaN;
  plotplm(difo,lon*pi/180,lat*pi/180,1)
  hold on
  p3=plot(x,y,'o');
  set([p3 p1 p4],'MarkerE','k','MarkerF','k','MarkerS',2)
  xl(4)=title('Difference > 10^{-12}');
  caxis(ck)
  cb(3)=colorbar('hor');
  set([xl yl],'FontS',15)
  longticks(ah,2)
  movev(cb,-.03)
  kelicol
  fig2print(gcf,'landscape')
  figdisp('plm2xyz_demo')
elseif strcmp(lmcosi,'demo3')
  % Load some constants - they happen to be the same for both models
  a96=fralmanac('a_EGM96','Earth');
  GM96=fralmanac('GM_EGM96','Earth');
  a2008=fralmanac('a_EGM2008','Earth');
  GM2008=fralmanac('GM_EGM2008','Earth');

  % Bypass FRALMANAC for efficiency
  load(fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS','SHM'))

  % Load the geopotential coefficients for EGM96
  v96=SHM.EGM96;
  % Load the geopotential coefficients for EGM2008
  v2008=SHM.EGM2008_ZeroTide;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOP ROW
  
  % Take out the useless error terms and go up to 360 only
  v96=v96(1:addmup(360)-addmup(1),1:4); 
  v2008=v2008(1:addmup(360)-addmup(1),1:4);

  % Convert to free-air gravity anomaly wrt to the C20 term, bypassing PLM2POT
  v96(1,3)=0; 
  v96(:,3:4)=v96(:,3:4)*GM96/a96.*repmat([(v96(:,1)-1)/a96],1,2);
  v2008(1,3)=0;
  v2008(:,3:4)=v2008(:,3:4)*GM2008/a2008.*repmat([(v2008(:,1)-1)/a2008],1,2);
  
  % Expand, convert to mGal
  v96xyz=plm2xyz(v96)*1e5;
  disp(' ')
  v2008xyz=plm2xyz(v2008)*1e5;
  disp(' ')
  
  figure(gcf); clf
  [ah,ha]=krijetem(subnum(2,2));
  axes(ah(1)); imagef([],[],v96xyz); caxis([-100 100]); axis image
  cb(1)=colorbar('hor'); axes(cb(1)); xlabel('EGM96 free-air anomaly (mgal)')
  axes(ah(2)); imagef([],[],v2008xyz); caxis([-100 100]); axis image
  cb(2)=colorbar('hor'); axes(cb(2)); xlabel('EGM2008 free-air anomaly (mgal)')
  shrink(cb(1:2),2,1.5); movev(cb(1:2),-.075)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOTTOM ROW
  % Load the geopotential coefficients for EGM96 - again
  v96=SHM.EGM96;
  % Load the geopotential coefficients for EGM2008
  v2008=SHM.EGM2008_ZeroTide;
  % Take out the useless error terms but increase resolution
  v96=v96(:,1:4); 
  v2008=v2008(1:addmup(1024)-addmup(1),1:4);

  % Convert to free-air gravity anomaly wrt to the C20 term, bypassing PLM2POT
  v96(1,3)=0; 
  v96(:,3:4)=v96(:,3:4)*GM96/a96.*repmat([(v96(:,1)-1)/a96],1,2);
  v2008(1,3)=0;
  v2008(:,3:4)=v2008(:,3:4)*GM2008/a2008.*repmat([(v2008(:,1)-1)/a2008],1,2);
  
  % Now specify a certain longitude-latitude grid for partial expansion
  c11cmn=[76 14 160 -32];

  % Expand, convert to mGal 
  v96xyz=plm2xyz(v96,[],c11cmn)*1e5;
  disp(' ')
  v2008xyz=plm2xyz(v2008,[],c11cmn)*1e5;
  disp(' ')
    
  axes(ah(3)); imagef([],[],v96xyz); caxis([-100 100]); axis image
  cb(3)=colorbar('hor'); axes(cb(3)); xlabel('EGM96 free-air anomaly (mgal)')
  axes(ah(4)); imagef([],[],v2008xyz); caxis([-100 100]); axis image
  cb(4)=colorbar('hor'); axes(cb(4)); xlabel('EGM2008 free-air anomaly (mgal)')
  shrink(cb(3:4),2,1.5); movev(cb(3:4),-.075)
  
  fig2print(gcf,'landscape')
  
elseif strcmp(lmcosi,'demo4')
  % Bypass FRALMANAC for efficiency
  load(fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS','SHM'))

  % Load the topography coefficients for GTM3AR
  vgtm=SHM.GTM3AR;
  % Load the topography coefficients for EGM2008
  v2008=SHM.EGM2008_Topography;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOP ROW
  
  % Take out the useless error terms and go up to 720 only
  vgtm=vgtm(1:addmup(720),1:4); 
  v2008=v2008(1:addmup(720),1:4);

  % Take out the C00-C20 terms to ignore mean and first-order flattening
  vgtm(1:addmup(2),3:4)=0; 
  v2008(1:addmup(2),3:4)=0;

  % Expand, convert to km
  vgtmxyz=plm2xyz(vgtm)*1e-3;
  v2008xyz=plm2xyz(v2008)*1e-3;
  
  figure(gcf); clf
  [ah,ha]=krijetem(subnum(2,2));
  axes(ah(1)); imagef([],[],vgtmxyz); caxis([-5 5]); axis image
  cb(1)=colorbar('hor'); axes(cb(1)); xlabel('GTM3AR topography (km)')
  axes(ah(2)); imagef([],[],v2008xyz); caxis([-5 5]); axis image
  cb(2)=colorbar('hor'); axes(cb(2)); xlabel('EGM2008 topography (km)')
  shrink(cb(1:2),2,1.5); movev(cb(1:2),-.075)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOTTOM ROW
  % Load the topography coefficients for GTM3AR - again
  vgtm=SHM.GTM3AR;
  % Load the topography coefficients for EGM2008
  v2008=SHM.EGM2008_Topography;
  % Take out the useless error terms but increase resolution
  vgtm=vgtm(:,1:4); 
  v2008=v2008(1:addmup(1024),1:4);

  % Convert to free-air gravity anomaly wrt to the C20 term, bypassing PLM2POT
  vgtm(1:addmup(2),3:4)=0; 
  v2008(1:addmup(3:4),3:4)=0;
  
  % Now specify a certain longitude-latitude grid for partial expansion
  c11cmn=[76 14 160 -32];

  % Expand, convert to km
  vgtmxyz=plm2xyz(vgtm,[],c11cmn)*1e-3;
  v2008xyz=plm2xyz(v2008,[],c11cmn)*1e-3;
  
  axes(ah(3)); imagef([],[],vgtmxyz); caxis([-5 5]); axis image
  cb(3)=colorbar('hor'); axes(cb(3)); xlabel('GTM3AR topography (km)')
  axes(ah(4)); imagef([],[],v2008xyz); caxis([-5 5]); axis image
  cb(4)=colorbar('hor'); axes(cb(4)); xlabel('EGM2008 topography (km)')
  shrink(cb(3:4),2,1.5); movev(cb(3:4),-.075)
  
  fig2print(gcf,'landscape')
  % I suppose this shows that the two C20 coefficients are very
  % different! Look at Australia.
elseif strcmp(lmcosi,'demo5')
  % Checks that the latest modifications have been put in right
  % Load some constants - they happen to be the same for both models
  a96=fralmanac('a_EGM96','Earth');
  GM96=fralmanac('GM_EGM96','Earth');

  % Bypass FRALMANAC for efficiency
  load(fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS','SHM'))

  % Load the geopotential coefficients for EGM96
  v96=SHM.EGM96;

  % Take out the useless error terms but increase resolution
  v96=v96(:,1:4); 

  % Convert to free-air gravity anomaly wrt to the C20 term, bypassing PLM2POT
  v96(1,3)=0; 
  v96(:,3:4)=v96(:,3:4)*GM96/a96.*repmat([(v96(:,1)-1)/a96],1,2);

  % In chunks of 72
  r1=plm2xyz(v96,[],[],90);
  % In chunks of 144
  r2=plm2xyz(v96,[],[],180);
  % All at once
  r3=plm2xyz(v96,[],[],360);
  
  % Check that the results are identical
  difer(r1-r2); difer(r2-r3); difer(r3-r1)
elseif strcmp(lmcosi,'demo6')
  v=fralmanac('EGM2008_Topography','SHM');
  v=v(1:addmup(1000)-addmup(v(1)-1),:);
  [r,lon,lat,Plm]=plm2xyz(v,[],[8 45 20 37]);
  if nargout==0
    plotplm(r,lon*pi/180,lat*pi/180,4);
  end
  % Now with Dongs' which looks good, passed all tests 4/21/2010
  % [rdw,londw,latdw]=plm2xyzdw(v,[],[45 37 8 20],0,0);
  
  % Prepare output
  vars={r,lon,lat,Plm};
  varargout=vars(1:nargout);
elseif strcmp(lmcosi,'demo7')
  % Load the model
  load(fullfile(getenv('IFILES'),...
		       'EARTHMODELS','POMME-4','pomme-4.2s-nosecular.mat'));
  % Bandwidth
  L=72;
  
  % Restrict the model to degree L
  lmcosi=lmcosi(1:addmup(L)-addmup(lmcosi(1)-1),:);
  
  % Convert to radial-component magnetic field on the reference surface
  lmcosip=plm2mag(lmcosi);
  
  % Generate a random set of locations inside a spherical patch
  Nd=1000;
  TH=15;
  phi0=15;
  theta0=70;
  [lon,lat]=randpatch(Nd,TH,phi0,theta0);

  % Perform the expansion to the complete grid
  [r1,long,latg]=plm2xyz(lmcosip,5);

  [LONG,LATG]=meshgrid(long,latg);

  % Now need to convert lmcosi into a format ready for the degree and
  % order ordering implicit in ylm - using ADDMON
  [~,~,~,~,~,~,~,~,~,ronm]=addmon(L);

  % Watch for the degree 0 term which isn't in there
  lmcosi=[0 0 0 0 ; lmcosi];
  lmcosip=[0 0 0 0 ; lmcosip];

  % And compare with a full-on gridded version using YLM
  [Y,theta,phi,dems,dels]=ylm([0 L],[],(90-latg)*pi/180,long*pi/180);
  % Don't forget to fix the phase and the normalization
  Y=[Y.*repmat((-1).^dems,1,size(Y,2))]*sqrt(4*pi);
  % Expand and reshape
  r2=reshape(Y'*...
	     lmcosip(2*size(lmcosip,1)+ronm),length(latg),length(long));

  % Check again using the irregular version of PLM2XYZ
  r3=reshape(plm2xyz(lmcosip,LATG(:),LONG(:)),length(latg),length(long));

  % These guys are not identical, but many of them are... and the rest
  % are very close
  difer(r1-r2,6)
  difer(r1-r3,6)
  difer(r2-r3,6)

  tic
  % Now check the expansion on the irregular set
  r4=plm2xyz(lmcosip,lat,lon);
  toc
  
  tic
  % And compare with a full-on irregular version using YLM
  [Y,theta,phi,dems,dels]=...
      ylm([0 L],[],(90-lat)*pi/180,lon*pi/180,[],[],[],1);
  % Don't forget to fix the phase and the normalization
  Y=[Y.*repmat((-1).^dems,1,size(Y,2))]*sqrt(4*pi);
  
  % Then perform the expansion watching the extra phase factor
  r5=Y'*lmcosip(2*size(lmcosip,1)+ronm);
  toc
  
  % Check the result
  difer(r4-r5,6)
  
  % Check it out using SCATTER.
  
  % Plot CRUSTAL model using PLOTPLM, compare with SAARIMAKI1
  lmcosip(1:addmup(15),3:4)=0;
  [r1,long,latg]=plm2xyz(lmcosip,0.5);
  % Obviously need to check the units here
  imagefnan([0 90],[360 -90],r1,[],[-1000 1000]); axis tight
  plotcont
  minmax(r1/1000)
elseif strcmp(lmcosi,'demo8')
  v=fralmanac('EGM2008_ZeroTide','SHM');
  % Note that gravity does not start at zero
  % Geoid = 3 % Free-air gravity anomaly = 2
  v=plm2pot(v(1:addmup(720)-addmup(v(1)-1),:),[],[],[],3);
  [r,lon,lat,Plm]=plm2xyz(v);
  plotplm(r)
  % Prepare output
  vars={r,lon,lat,Plm};
  varargout=vars(1:nargout);
end
