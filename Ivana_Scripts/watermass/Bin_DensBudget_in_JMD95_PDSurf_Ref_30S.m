%  cd /data1/project/gyre/SOSE_ITER60/pd_surf_scripts
%  Bin_DensBudget_in_JMD95_PDSurf_Ref_30S


clear all
close all
warning off

rho0=1035;

IT1=input('FROM which iteration number 1, 74, 147....  ? ')
IT2=input('UNTIL which iteration number 73, 146, 219....  ? ')
lt=IT2-IT1+1;  % length of the time period considered

NY_LIM=288;   % input('What is the index of the northermost point until which you integrate?   YC=-30 ')

cd /data1/project/gyre/SOSE_ITER22/GRID
hFacC = rdmds('hFacC');
DRF=rdmds('DRF');
NZ=length(DRF);
RAC=rdmds('RAC');

% coordinates
ILON = [1:2159];        NX=length(ILON);
ILAT = [1:NY_LIM];
IDEP = [1:NZ];


vol=zeros(NX,NY_LIM,NZ);
display('Estimate the net volume of each numerical cell')
for kz=1:NZZ
vol(:,:,kz)=hFacC(ILON,ILAT,kz).*RAC(ILON,ILAT).*squeeze(DRF(kz));
end

% define binning intervals
deltaRHO=0.1
deltaRHO2=0.5*deltaRHO;
rhobin =[1025.25:deltaRHO:1028];
nbins=floor(length(rhobin)); nbins1=nbins-1;



ic=0;
for      it = IT1:IT2
ic=ic+1;


cd  /data2/project/gyre-data/SOSE_ITER60/STATE
T = rdmds('THETA','rec',it);
S = rdmds('SSALT','rec',it);
T=THETA(ILON,ILAT,IDEP);
S=SALT(ILON,ILAT,IDEP); 

cd /data1/project/gyre/SOSE_ITER60/pd_surf_scripts
PD = densjmd95(S,T,zeros(size(S)));

display('Load all the terms from potential density equation : tendency Rho...')
cd    /data2/project/gyre-data/SOSE_ITER60/PDbudget_RefSurf
tend=rdmds('TEND_R_RefSurf',it);
adv=rdmds('ADV_R_RefSurf',it);
surf=rdmds('DIV_R_SurfFLX_RefSurf',it);
DIFE=rdmds('DIFE_R_RefSurf',it);
DIFI=rdmds('DIFI_R_RefSurf',it);
DIFKP=rdmds('DIFKP_R_RefSurf',it);

tend=tend(ILON,ILAT,IDEP);
adv=adv(ILON,ILAT,IDEP);
surf=surf(ILON,ILAT,IDEP);
DIFE=DIFE(ILON,ILAT,IDEP);
DIFI=DIFI(ILON,ILAT,IDEP);
DIFKP=DIFKP(ILON,ILAT,IDEP);

diab=DIFE+DIFI+DIFKP+surf;
adv(:,:,1)=tend(:,:,1)-diab(:,:,1); % because "WVEL(1)" is missing for the free surface


figure(1); clf
plot(squeeze(sum(sum(tend))),'r'); hold on
plot(squeeze(sum(sum(diab))),'b'); plot(squeeze(sum(sum(adv))),'k');
plot(squeeze(sum(sum(tend))-sum(sum(adv))-sum(sum(diab))),'o--','linewidth',2,'color','k');  
legend('tend','diab','adv','res','sum');
title(['Potential Density budget in the Southern Ocean,   iter= ',num2str(it)]); pause(0.01)


[n m l]=size(PD);  nml=n*m*l;
bindRdt=zeros(n,m,nbins1); bin_tend=bindRdt;   bin_adv=bin_tend; bin_surf=bin_tend;
bin_DIFE=bin_tend;         bin_DIFI=bin_tend;  bin_DIFKP=bin_tend;


% Bin all individual terms of potential density budget
for nb=1:nbins1
bindepth=zeros(size(PD));
bindepthr=reshape(bindepth,1,nml);
PDr=reshape(PD,1,nml);

bindepthr(find(PDr>=rhobin(nb) & PDr<rhobin(nb+1)))=1;
bindepth=reshape(bindepthr,n,m,l);

bin_tend(:,:,nb)=sum(((tend.*bindepth).*vol),3);
bin_adv(:,:,nb)=sum(((adv.*bindepth).*vol),3);
bin_surf(:,:,nb)=sum(((surf.*bindepth).*vol),3);
bin_DIFE(:,:,nb)=sum(((DIFE.*bindepth).*vol),3);
bin_DIFI(:,:,nb)=sum(((DIFI.*bindepth).*vol),3);
bin_DIFKP(:,:,nb)=sum(((DIFKP.*bindepth).*vol),3);
end %

%
% estimate FORMATION RATES in Sverdrups
%

Bin_Tend=diff(bin_tend,1,3)/(1e6*deltaRHO);
Bin_Adv=diff(bin_adv,1,3)/(1e6*deltaRHO);
Bin_Surf_Form=diff(bin_surfRflx,1,3)/(1e6*deltaRHO);

Bin_DIFE=diff(bin_DIFE,1,3)/(1e6*deltaRHO);
Bin_DIFI=diff(bin_DIFI,1,3)/(1e6*deltaRHO);
Bin_DIFKP=diff(bin_DIFKP,1,3)/(1e6*deltaRHO);




