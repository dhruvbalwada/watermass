%  cd  /data1/project/gyre/SOSE_with_Ryan
%  Temp_Budget

clear all
close all

display('GET GRID DETAILS NEEDED')
cd /data2/project/gyre-data/SOSE_with_Ryan
RAC = rdmds('RAC');
hFacC = rdmds('hFacC');
hFacW = rdmds('hFacW');
hFacS = rdmds('hFacS');
DX = rdmds('DXG');
DY = rdmds('DYG');
DZ = rdmds('DRF');  NZZ=length(DZ);
DRF=rdmds('DRF');
DRF=squeeze(DRF);
mask=squeeze(hFacC(:,[1:end-1],:));

ZC = (cumsum(squeeze(DRF))+(cumsum([0 squeeze(DRF(1:NZZ-1))']))')./(2);

YC=rdmds('YC');
XC=rdmds('XC');
maskC = hFacC;maskC(maskC~=0)=1;



[NX NY NZ] = size(hFacC);
VOL = zeros(NX,NY,NZ+1,'single');
AREAWEST = VOL;AREASOUTH = VOL; AREACELL = VOL;

% added extra layer (NZ+1) of "ground" for convenience with vertical derivatives
for k = 1:NZ
  VOL(:,:,k) = hFacC(:,:,k).*RAC.*DZ(k);
end
RecipVol = VOL;RecipVol(RecipVol==0)=inf;RecipVol = 1./RecipVol;
ILON = 1;        nx=length(ILON);
ILAT = [1:NY-1]; ny=length(ILAT); 
IDEP = [1:NZ];   nz=length(IDEP); 

% ~~~~~~~~~~~     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

display(' READ IN HEAT BUDGET TERMS ')

% ADVECTIVE TERMS -- zonal, then meridional, then vertical
% Q=rdmds('DiagGAD-T.0009473760','rec',3);
%  tm1 = Q(ILON+1,ILAT,IDEP);% 'ADVx_TH '
%  tm2 = Q(ILON  ,ILAT,IDEP);% 'ADVx_TH '
Ax  = zeros(nx,ny,nz);

Q=rdmds('DiagGAD-T.0009473760','rec',4);
 tm1 = Q(ILON,ILAT+1,IDEP);% 'ADVy_TH '
 tm2 = Q(ILON,ILAT  ,IDEP);% 'ADVy_TH '
Ay  = (tm1 - tm2);


Q=rdmds('DiagGAD-T.0009473760','rec',2);
Q(:,:,nz+1)=zeros(nx,ny+1,1);
 tm1 = Q(ILON,ILAT,IDEP  );% 'ADVr_TH '
 tm2 = Q(ILON,ILAT,IDEP+1);% 'ADVr_TH '
Az = (tm1 - tm2);


% DIFFUSIVE TERMS -- zonal, then meridional, then vertical
% Q=rdmds('DiagGAD-T.0009473760','rec',6);
% tm1 = Q(ILON+1,ILAT,IDEP);% 'DFxE_TH '
% tm2 = Q(ILON  ,ILAT,IDEP);% 'DFxE_TH '
Dx  = zeros(nx,ny,nz);

Q=rdmds('DiagGAD-T.0009473760','rec',7);
 tm1 = Q(ILON,ILAT+1,IDEP);% 'DFyE_TH '
 tm2 = Q(ILON,ILAT  ,IDEP);% 'DFyE_TH '
Dy  = (tm1 - tm2);


Q=rdmds('DiagGAD-T.0009473760','rec',5);
Q(:,:,nz+1)=zeros(nx,ny+1,1);
 tm1 = Q(ILON,ILAT,IDEP  );% 'DFrE_TH '
 tm2 = Q(ILON,ILAT,IDEP+1);% 'DFrE_TH '
Dze = (tm1 - tm2);%explicit vertical diff



Q=rdmds('DiagGAD-T.0009473760','rec',8);
Q(:,:,nz+1)=zeros(nx,ny+1,1);
 tm1 = Q(ILON,ILAT,IDEP );% 'DFrI_TH '
 tm2 = Q(ILON,ILAT,IDEP+1);% 'DFrI_TH '
Dzi = (tm1 - tm2);%implicit vertical diff


day=24*60*60;
% AND TENDENCY dT/dt
Q=rdmds('DiagGAD-T.0009473760','rec',1);
TT = Q(ILON,ILAT,IDEP)/day;  %'TOTTTEND'

% forcing 
clear Q                              
Q=rdmds('DiagSurfT.0009473760','rec',1); 

% Q is surfaceForcingT *HeatCapacity_Cp*rUnit2mass;   gT is +surfaceForcingT/drf(1)
% 1 / (Cp*rho_nil*DRF(1) );

fac= 1./(3994*1035*DRF(1)); 
TFLX =  fac.*Q(ILON,ILAT);        % [deg/sec]  this is just to change the units on net heat flux and SW flux , these are all surface fields!



% FLUXES DETERMINED IN 
% gad_calc_rhs.F
%        gTracer(i,j,k,bi,bj)=gTracer(i,j,k,bi,bj)
%     &   -recip_hFacC(i,j,k,bi,bj)*recip_drF(k)
%     &   *recip_rA(i,j,bi,bj)*recip_deepFac2C(k)*recip_rhoFacC(k)
%     &   *( (fZon(i+1,j)-fZon(i,j))
%     &     +(fMer(i,j+1)-fMer(i,j))
%     &     +(fVerT(i,j,kDown)-fVerT(i,j,kUp))*rkSign
%     &     -localT(i,j)*( (uTrans(i+1,j)-uTrans(i,j))*advFac
%     &                   +(vTrans(i,j+1)-vTrans(i,j))*advFac
%     &                   +(rTransKp1(i,j)-rTrans(i,j))*rAdvFac
% WHERE
%  fVerT = KPPg_TH*maskUp*rhoFacF + DFrE_TH*maskUp + ADVr_TH;
%  fMer  = DFyE_TH*rhoFacF + ADVy_TH
%  fZon  = DFxE_TH*rhoFacF + ADVx_TH
% NOW GET TENDENCY FROM SUM OF RHS

gTracer = zeros(size(Ax));
gTracer = -RecipVol(ILON,ILAT,IDEP).*(Ax+Dx+Ay+Dy+Az+Dze+Dzi);

clear fac
gTracer(:,:,1) = gTracer(:,:,1) + TFLX;  

ADV=zeros(size(Ax)); DIFE=ADV;  DIFI=ADV; RES=ADV;
RecVol=RecipVol(ILON,ILAT,IDEP);


ADV = -RecVol.* (Ax+Ay+Az);     %   U dot grad T  = U \cdot \nabla T
DIFE =  -RecVol.*(Dx+Dy+Dze); %  Sum of diffusive terms ~ \nabla \kappa \nabla T
DIFI =  -RecVol.*Dzi;
RES =    TT - gTracer;                                 % is almost zero below top couple layers 

DIV_FLX=zeros(size(Ax));
DIV_FLX(:,:,1) = TFLX;  
DIAB=DIV_FLX+DIFE+DIFI;

%      TEST 

IIXX=ILON;
JJYY=ILAT;
KKZZ=IDEP;
maskC=hFacC;
maskC(find(hFacC~=1))=0;

mm=squeeze(maskC(IIXX,JJYY,KKZZ)); mm=double(mm); 

figure(1)
clf
plot(squeeze(sum(mm.*squeeze(double(TT(IIXX,JJYY,KKZZ))))),'ro-'); hold on
plot(squeeze(sum(mm.*squeeze(ADV(IIXX,JJYY,KKZZ)))),'x-','linewidth',2,'color','k'); 
plot(squeeze(sum(mm.*squeeze(DIAB(IIXX,JJYY,KKZZ)))),'bo-'); grid on; title(['TEMP budget'])
plot(squeeze(sum(-mm.*squeeze(RES(IIXX,JJYY,KKZZ)))),'o-','linewidth',2,'color','g'); grid on ; axis('tight')  ;
legend('T tend','div of T adv','T form','T sum')

%~~~~~~~~~~~~~~~~~
t=rdmds('TTtave.0009473760');
t=squeeze(t);
t=t([1:end-1],:);

yy=YC(IIXX,JJYY);
zz=ZC;

res=TT-ADV-DIAB;
lim=3e-3
surf_forc=zeros(size(squeeze(ADV)));
surf_forc(:,1)=TFLX;

figure(2)
clf
subplot(321); hold on 
[h]=pcolor(yy,zz,day*squeeze(double(TT))'); shading flat; [h]=colorbar; caxis([-1 1]*lim);   colormap('redblue'); 
set(h,'position',[0.94    0.709    0.015    0.2]);  axis('tight'); AXX=axis; title(['TOTTEND  [x day  ^{\circ}]  ']);
hold on; contour(yy,zz,mask',[0.1 0.1],'linewidth',2,'color','k'); contour(yy,zz,t',12,'color','k')
set(gca,'ydir','reverse'); 
subplot(322)
pcolor(yy,zz,day*squeeze(double(ADV))'); shading flat;  caxis([-1 1]*lim); title(['Tendency from advection  [x day] ']); 
axis(AXX); hold on;  contour(yy,zz,mask',[0.1 0.1],'linewidth',2,'color','k');  contour(yy,zz,t',12,'color','k')
set(gca,'ydir','reverse')
subplot(323)
pcolor(yy,zz,day*squeeze(double(DIAB))'); shading flat; caxis([-1 1]*lim);  hold on ; 
contour(yy,zz,mask',[0.1 0.1],'linewidth',2,'color','k'); contour(yy,zz,t',12,'color','k')
title(['Tendency from diffusion T  [x day] ']);  axis(AXX); 
set(gca,'ydir','reverse')
subplot(325); hold on 
pcolor(yy,zz,day*squeeze(double(ADV))'); shading flat; caxis([-1 1]*lim);  title(['ADV  [x day] ']);  axis(AXX); hold on
contour(yy,zz,mask',[0.1 0.1],'linewidth',2,'color','k'); contour(yy,zz,t',12,'color','k')
set(gca,'ydir','reverse')
subplot(324); hold on 
pcolor(yy,zz,day*squeeze(double(surf_forc))'); shading flat; title(['Surface Forcing [x day] ']);  axis(AXX); hold on
contour(yy,zz,mask',[0.1 0.1],'linewidth',2,'color','k'); contour(yy,zz,t',12,'color','k')
set(gca,'ydir','reverse'); caxis([-1 1]*lim); 
subplot(326)
pcolor(yy,zz,day*squeeze(double(res))'); shading flat; caxis([-1 1]*lim*1e-4); hold on;  contour(yy,zz,t',12,'color','k')
contour(yy,zz,mask',[0.1 0.1],'linewidth',2,'color','k');
 title([' Residual  [x day] - NOTE scale different by factor of 1e4']);  axis(AXX); [h2]=colorbar;  contour(yy,zz,t',12,'color','k')
set(h2,'position',[0.94    0.11    0.015    0.2]); set(gca,'ydir','reverse')

