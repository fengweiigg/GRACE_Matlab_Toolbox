
function dataDDK=gmt_destriping_ddk(number,data)

% DDK filtering
% 
% INPUT:
%   number            the type of DDK filter
%   data              spherical harmonic coefficients before filtering
%
% OUTPUT:
%   grid_filter       equi-angular grid N*2N, N=180 or 720
% 
% DDK1d11: filtered with inverse signal degree power law 1e11*deg^4 (DDK5) weakest smoothing
% DDK5d11: filtered with inverse signal degree power law 5e11*deg^4 (DDK4)         |
% DDK1d12: filtered with inverse signal degree power law 1e12*deg^4 (DDK3)         |
% DDK1d13: filtered with inverse signal degree power law 1e13*deg^4 (DDK2)         |
% DDK1d14: filtered with inverse signal degree power law 1e14*deg^4 (DDK1) strongest smoothing

% FENG Wei 18/10/2016
% 
% State Key Laboratory of Geodesy and Earth's Dynamics
% Institute of Geodesy and Geophysics, Chinese Academy of Sciences
% fengwei@whigg.ac.cn

% This function is created based on the code from Roelof Rietbroek.
% Copyright Roelof Rietbroek 2016
% This software is licensed under the MIT license see https://github.com/strawpants/GRACE-filter/blob/master/LICENSE
% URL: https://github.com/strawpants/GRACE-filter


switch number
    case 1 %strongest DDK1
        file='Wbd_2-120.a_1d14p_4';
    case 2 % DDK2
        file='Wbd_2-120.a_1d13p_4';
    case 3 % DDK3
        file='Wbd_2-120.a_1d12p_4';
    case 4 % DDK4
        file='Wbd_2-120.a_5d11p_4';
    case 5 % DDK5
        file='Wbd_2-120.a_1d11p_4';
	case 6
		file='Wbd_2-120.a_5d10p_4';
	case 7
		file='Wbd_2-120.a_1d10p_4';
	case 8
		file='Wbd_2-120.a_5d9p_4';		
end

% rep=pwd;
% cd('C:\_PROGZ\GRACE\DDK\filtercoef\')
dat=read_BIN([file]);
% cd(rep)

% Block interpretation
% size ith block
clear sz nend nstart
sz(1)=dat.blockind(1);
nstart(1)=1;
nend(1)=1+sz(1)^2-1;
for ij=2:dat.nblocks
    sz(ij)=dat.blockind(ij)-dat.blockind(ij-1);
    nstart(ij)=1+sum(sz(1:ij-1).^2);
    nend(ij)=nstart(ij)+sz(ij).^2-1;
end

% GRACE dataset
nmax=size(data.C,1)-1;
ntime=size(data.C,3);

% Initialization
dataDDK=data;


% CAS de l'ordre 0
ordre=0;
% block d'ordre 0, degree 2-> 120
ij=ordre+1;
block=reshape(dat.pack1(nstart(ij):nend(ij)),sz(ij),sz(ij));
% on garde degree 2-> nmax
block=block(1:nmax-1,1:nmax-1);

% GRACE ordre 0 degree 2->nmax
for klm=1:ntime
    coef=squeeze(data.C(3:end,ordre+1,klm));
    %coefF=block*coef;

    % remplacement
    dataDDK.C(3:end,ij,klm)=block*coef;
end

ordre=1;
while ordre<nmax+1
    % block d'ordre ordre, degree 2-> 120
    ij=2*(ordre);
    blockC=reshape(dat.pack1(nstart(ij):nend(ij)),sz(ij),sz(ij));
    blockS=reshape(dat.pack1(nstart(ij+1):nend(ij+1)),sz(ij+1),sz(ij+1));
    % on garde degree 2-> nmax
    fin=min(nmax-1,nmax-ordre+1);
    blockC=blockC(1:fin,1:fin);
    blockS=blockS(1:fin,1:fin);
    % GRACE ordre ordre degree 2->nmax
    deb=max(3,ordre+1);
    
    for klm=1:ntime
        coefC=squeeze(data.C(deb:end,ordre+1,klm));
        coefS=squeeze(data.S(deb:end,ordre+1,klm));
%         coefCF=blockC*coefC;
%         coefSF=blockS*coefS;

        % remplacement
        dataDDK.C(deb:end,ordre+1,klm)=blockC*coefC;
        dataDDK.S(deb:end,ordre+1,klm)=blockS*coefS;
    end

    % increment
    ordre=ordre+1;
end
% imagesc(log10(abs(data.C(:,:,1))))
% set(gca,'Clim',[-12 0])
% colorbar
% figure(2)
% imagesc(log10(abs(dataDDK.C(:,:,1))))
% set(gca,'Clim',[-12 0])
% colorbar


    