
function nouv=ddk(number,data);
% ddk filtering

switch number
    case 1 %strongest
        file='Wbd_2-120.a_1d14p_4';
    case 2
        file='Wbd_2-120.a_1d13p_4';
    case 3
        file='Wbd_2-120.a_1d12p_4';
    case 4
        file='Wbd_2-120.a_5d11p_4';
    case 5
        file='Wbd_2-120.a_1d11p_4';
end

dat=read_BIN(file);

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
    coefF=block*coef;

    % remplacement
    dataDDK.C(3:end,ij,klm)=coefF;
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
        coefCF=blockC*coefC;
        coefSF=blockS*coefS;

        % remplacement
        dataDDK.C(deb:end,ordre+1,klm)=coefCF;
        dataDDK.S(deb:end,ordre+1,klm)=coefSF;
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


    