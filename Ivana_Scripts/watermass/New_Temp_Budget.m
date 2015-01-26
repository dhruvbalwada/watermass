%  cd  /data1/project/gyre/SOSE_ITER60/temp_budget_scripts
%  New_Temp_Budget

display('Run for iteration 100 ')
clear all
close all
if_sav=0;


IT1=input('FROM which iteration number 1, 74, 147....  ? ')
IT2=input('UNTIL which iteration number 73, 146, 219....  ? ')

if_sav=1
rho0=1035

display('GET GRID DETAILS NEEDED')
% GET GRID DETAILS NEEDED
cd /data1/project/gyre/SOSE_ITER22/GRID
RAC = rdmds('RAC');
hFacC = rdmds('hFacC');
hFacW = rdmds('hFacW');
hFacS = rdmds('hFacS');
DX = rdmds('DXG');
DY = rdmds('DYG');
DZ = rdmds('DRF');
DRF = rdmds('DRF');

YC=rdmds('YC');
XC=rdmds('XC');


[NX NY NZ] = size(hFacC);
VOL = zeros(NX,NY,NZ+1,'single');
AREAWEST = VOL;AREASOUTH = VOL; AREACELL = VOL;
% added extra layer (NZ+1) of "ground" for convenience with vertical derivatives
for k = 1:NZ
  VOL(:,:,k) = hFacC(:,:,k).*RAC.*DZ(k);
  AREAWEST(:,:,k)  = DY.*DZ(k).*hFacW(:,:,k);
  AREASOUTH(:,:,k) = DX.*DZ(k).*hFacS(:,:,k);
  AREACELL(:,:,k)  = RAC.*hFacC(:,:,k);
end
RecipVol = VOL;RecipVol(RecipVol==0)=inf;RecipVol = 1./RecipVol;

ILON = [1:2159];
ILAT = [1:319];
IDEP = [1:42];


DRF0(1)=0; 
DRF0(2:NZ)=DRF(1:NZ-1); DRF=squeeze(DRF);  DRF0=squeeze(DRF0); 
ZC = (cumsum(DRF')+cumsum(DRF0))./(2);
ZC(1:6)

% ~~~~~~~~~~~    DO TIME LOOP   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for   it=IT1:IT2      % 3 yrs of data 
[it]
display(' READ IN HEAT BUDGET TERMS ')
display('CALCULATE DENSITY BUDGET')
fname_int=[ 'ADVx_TH ' 'ADVy_TH ' 'ADVr_TH '  'DFxE_TH ' 'DFyE_TH ' 'DFrE_TH ' 'DFrI_TH ' 'KPPg_TH ' 'TOTTTEND'];
fname_stt =[ 'THETA   '  'UVEL    ' 'VVEL    ' 'WVEL    '];
fname_frc=['TFLUX' 'EXFswnet' 'oceFreez'];


inm1=[1  9  17 25 33  41 49  57  65];
inm2=[7 15  23 31 39 47 55  63  72];
LR=length(inm1)   % number of fields to read in 


cd  /data2/project/gyre-data/SOSE_ITER60/TEMPbdgt
Q=zeros(NX,NY,NZ+1,LR);   %extra z lvl for vert der.
for ifields=1:LR
 [fname_int(inm1(ifields):inm2(ifields))]
 tmp = rdmds([fname_int(inm1(ifields):inm2(ifields))],'rec',it);
 Q(:,:,[1:NZ],ifields) =tmp; tmp=zeros(NX,NY,NZ);
end 


% ADVECTIVE TERMS -- zonal, then meridional, then vertical
 tm1 = Q(ILON+1,ILAT,IDEP,1);% 'ADVx_TH '
 tm2 = Q(ILON  ,ILAT,IDEP,1);% 'ADVx_TH '
Ax  = (tm1 - tm2);
 tm1 = Q(ILON,ILAT+1,IDEP,2);% 'ADVy_TH '
 tm2 = Q(ILON,ILAT  ,IDEP,2);% 'ADVy_TH '
Ay  = (tm1 - tm2);
 tm1 = Q(ILON,ILAT,IDEP  ,3);% 'ADVr_TH '
 tm2 = Q(ILON,ILAT,IDEP+1,3);% 'ADVr_TH '
Az = (tm1 - tm2);

% DIFFUSIVE TERMS -- zonal, then meridional, then vertical
 tm1 = Q(ILON+1,ILAT,IDEP,4);% 'DFxE_TH '
 tm2 = Q(ILON  ,ILAT,IDEP,4);% 'DFxE_TH '
Dx  = (tm1 - tm2);
 tm1 = Q(ILON,ILAT+1,IDEP,5);% 'DFyE_TH '
 tm2 = Q(ILON,ILAT  ,IDEP,5);% 'DFyE_TH '
Dy  = (tm1 - tm2);
 tm1 = Q(ILON,ILAT,IDEP  ,6);% 'DFrE_TH '
 tm2 = Q(ILON,ILAT,IDEP+1,6);% 'DFrE_TH '
Dze = (tm1 - tm2);%explicit vertical diff
 tm1 = Q(ILON,ILAT,IDEP  ,7);% 'DFrI_TH '
 tm2 = Q(ILON,ILAT,IDEP+1,7);% 'DFrI_TH '
Dzi = (tm1 - tm2);%implicit vertical diff
 tm1 = Q(ILON,ILAT,IDEP  ,8);% 'KPPg_TH '
 tm2 = Q(ILON,ILAT,IDEP+1,8);% 'KPPg_TH '
KP = (tm1- tm2);% KPP term


% AND TENDENCY dT/dt
TT = Q(ILON,ILAT,IDEP,9)./86400;  %'TOTTTEND'
clear Q

% ALSO NEED UVW TRANSPORT SO FROM STATE
cd  /data2/project/gyre-data/SOSE_ITER60/STATE
im1=[1 9  17 25 ];
im2=[5 12 20 28 ];
lns=length(im1);
Q = zeros(NX,NY,NZ+1,lns,'single');
for ifields=1:lns
 [fname_stt(im1(ifields):im2(ifields))]
 Q(:,:,[1:NZ],ifields) = rdmds([fname_stt(im1(ifields):im2(ifields))],'rec',it);
end

THETA = Q(ILON,ILAT,IDEP,1);%'THETA   '
% GET DIVERGENCE, HORIZONTAL THEN VERTICAL
 tm1 =   Q(ILON+1,ILAT,IDEP,2).*AREAWEST(ILON+1,ILAT,IDEP);% UVEL
 tm2 =   Q(ILON  ,ILAT,IDEP,2).*AREAWEST(ILON,ILAT,IDEP);% UVEL
DU  =   (tm1 - tm2);
 tm1 = Q(ILON,ILAT+1,IDEP,3).*AREASOUTH(ILON,ILAT+1,IDEP);% VVEL
 tm2 = Q(ILON,ILAT  ,IDEP,3).*AREASOUTH(ILON,ILAT,IDEP);% VVEL
DV = (tm1 - tm2);
 tm1 = Q(ILON,ILAT,IDEP  ,4).*AREACELL(ILON,ILAT,IDEP);% WVEL
 tm1(:,:,1)=zeros(NX-1,NY-1,1); % <<<<<<-------------- set w=0 at the surafce 
 tm2 = Q(ILON,ILAT,IDEP+1,4).*AREACELL(ILON,ILAT,IDEP+1);% WVEL
DW = (tm1 - tm2); 

%  if (min(IDEP)==1) % get surface forcing if surface in region
%  Q   = rdmds([fname_frc],'rec',[1]); % TFLUX


clear Q                              
cd  /data2/project/gyre-data/SOSE_ITER60/TFORCING
Q=rdmds('TFLUX.0000000060','rec',it); 
SW=rdmds('EXFswnet.0000000060','rec',it);  % 


% Q is surfaceForcingT *HeatCapacity_Cp*rUnit2mass;   gT is +surfaceForcingT/drf(1)


fac= 1./(3994*1035*DRF(1));       % 1 / (Cp*rho_nil*DRF(1) );
TFLX =  fac.*Q(ILON,ILAT);        % [deg/sec]  this is just to change the units on net heat flux and SW flux , these are all surface fields!
Qsw=    fac.*SW(ILON,ILAT); 


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


SWFRAC =   [   0.7695     0.0975      0.0600   0.0349       0.0192    0.0104   0.0050      0.0022   8.6887e-04   2.9287e-04 ...
      8.3983e-05    2.5737e-05  0      0 ];


gTracer = zeros(size(Ax));
gTracer = -RecipVol(ILON,ILAT,IDEP)   ...
          .*( (Ax+Dx+Ay+Dy+Az+Dze+KP) ... %
               + Dzi                  ...
              - THETA.*(DU+DV+DW)    );

gTracer(:,:,1) = gTracer(:,:,1) + TFLX-(1-SWFRAC(1))*Qsw(ILON,ILAT);  

for k = 2:12
gTracer(:,:,k) = gTracer(:,:,k) +SWFRAC(k)*Qsw(ILON,ILAT)*(DRF(1)/DRF(k));
end

ADV1 =zeros(size(Ax)); ADV2=zeros(size(Ax)); ADV=zeros(size(Ax)); DIFE=ADV; DIFKP=ADV; DIFI=ADV; RES=ADV;
RecVol=RecipVol(ILON,ILAT,IDEP);
ADV1x=zeros(size(Ax));  ADVx=ADV1x;   ADV1y=zeros(size(Ay));  ADVy=ADV1y; ADV1z=zeros(size(Az));  ADVz=ADV1z; 
ADV2x=zeros(NX,NY,NZ);  ADV2y=ADV2x;  ADV2z=ADV2x;


ADV1 = -RecVol.* (Ax+Ay+Az);     %   U dot grad T  = U \cdot \nabla T
ADV2 =  RecVol.*(THETA(ILON,ILAT,IDEP).*(DU+DV+DW));  %   T times div U = T \nabla \cdot U
ADV=ADV1+ADV2;

ADV1x = -RecVol.* Ax; 
ADV2x =  RecVol.*(THETA(ILON,ILAT,IDEP).*DU);  
ADVx=ADV1x+ADV2x;

ADV1y = -RecVol.* Ay; 
ADV2y =  RecVol.*(THETA(ILON,ILAT,IDEP).*DV);  
ADVy=ADV1y+ADV2y;

ADV1z = -RecVol.* Az; 
ADV2z =  RecVol.*(THETA(ILON,ILAT,IDEP).*DW);  
ADVz=ADV1z+ADV2z;


DIFE =  -RecVol.*(Dx+Dy+Dze); %  Sum of diffusive terms ~ \nabla \kappa \nabla T
DIFKP = -RecVol.* KP;
DIFI =  -RecVol.*Dzi;
RES =    TT - gTracer;                                 % is almost zero below top couple layers 

DIV_FLX=zeros(size(Ax));
DIV_FLX(:,:,1)=TFLX-(1-SWFRAC(1))*Qsw(ILON,ILAT);
for k=2:12
DIV_FLX(:,:,k)=SWFRAC(k)*Qsw(ILON,ILAT)*(DRF(1)/DRF(k));
end

display('Save!')
if(if_sav==1)
cd    /data2/project/gyre-data/SOSE_ITER60/Tbudget_NEWEST
for isav=1:10
[isav]
  if(isav==1); fname='TEND_T' ; end
  if(isav==2); fname='ADV_T' ; end
  if(isav==3); fname='DIFE_T'  ; end  
  if(isav==4); fname='DIFI_T'  ; end  
  if(isav==5); fname='DIFKP_T'  ; end      
  if(isav==6); fname='DIV_TFLX';  end            
  if(isav==7); fname='RES_T'; end           
  if(isav==8); fname='ADV_Tx' ; end
  if(isav==9); fname='ADV_Ty' ; end
  if(isav==10); fname='ADV_Tz' ; end

   
if(it<10)
  fidd=fopen([fname,'.000000000' num2str(it) '.data'],'w','b');
  fidm=fopen([fname,'.000000000' num2str(it) '.meta'],'w','b');
end
if(it<100 & it>9)
  fidd=fopen([fname,'.00000000' num2str(it) '.data'],'w','b');
  fidm=fopen([fname,'.00000000' num2str(it) '.meta'],'w','b');
end
if(it>99)
  fidd=fopen([fname,'.0000000' num2str(it) '.data'],'w','b');
  fidm=fopen([fname,'.0000000' num2str(it) '.meta'],'w','b');
end

  if(isav==1); fwrite(fidd,TT,'single') ; end
  if(isav==2); fwrite(fidd,ADV,'single') ; end
  if(isav==3); fwrite(fidd,DIFE,'single') ; end      
  if(isav==4); fwrite(fidd,DIFI,'single');  end         
  if(isav==5); fwrite(fidd,DIFKP,'single') ; end    
  if(isav==6); fwrite(fidd,DIV_FLX,'single') ; end
  if(isav==7); fwrite(fidd,RES,'single') ; end
  if(isav==8); fwrite(fidd,ADVx,'single') ; end
  if(isav==9); fwrite(fidd,ADVy,'single') ; end
  if(isav==10); fwrite(fidd,ADVz,'single') ; end
              
fclose(fidd);
  fprintf(fidm,'%s\n',' nDims = [   3 ];');
  fprintf(fidm,'%s\n',' dimList = [');
  fprintf(fidm,'%s\n','           2159,    1, 2159,');
  fprintf(fidm,'%s\n','            319,    1,  319,');
  fprintf(fidm,'%s\n','             42,    1,   42');
  fprintf(fidm,'%s\n','          ];');
  fprintf(fidm, '%s%s%s\n', [' format = [ '''] , ['float32'] ,[''' ];']);
  fprintf(fidm,'%s\n',' nrecords = [     1 ];');    %  REC#?
  fprintf(fidm,'%s%s%s\n',' timeStepNumber = [    ', num2str(it),  ' ];');
fclose(fidm);
end  % isave loop 
end % if sav 


end % time
return



%      TEST 

IIXX=[1:2159];
JJYY=[1:300];
KKZZ=IDEP;
mm=squeeze(hFacC(IIXX,JJYY,KKZZ)); mm=double(mm); 

figure(1)
clf
subplot(211)
plot(squeeze(sum(sum(mm.*squeeze(double(TT(IIXX,JJYY,KKZZ)))))),'ro-'); hold on
plot(squeeze(sum(sum(mm.*squeeze(ADV(IIXX,JJYY,KKZZ))))),'x-','linewidth',2,'color','k'); 
plot(squeeze(sum(sum(mm.*squeeze(DIV_FLX(IIXX,JJYY,KKZZ))))),'bo-'); grid on; title(['TEMP budget   iter= ',num2str(it)])
plot(squeeze(sum(sum(-mm.*squeeze(RES(IIXX,JJYY,KKZZ))))),'o-','linewidth',2,'color','g'); grid on ; axis('tight')  ;
legend('S tend','S adv','S form','S sum')
%~~~~~~~~~~~~~~~~~

subplot(212)
plot(squeeze(sum(sum(-mm.*squeeze(RES(IIXX,JJYY,KKZZ))))),'o-'); grid on ;  hold on;  title('sum')
plot(squeeze(sum(sum(mm.*squeeze(ADV1(IIXX,JJYY,KKZZ))))),'x-','linewidth',2,'color','g'); 
plot(squeeze(sum(sum(mm.*squeeze(ADV2(IIXX,JJYY,KKZZ))))),'x-','linewidth',2,'color','r'); legend('SUM of the budget','ADV1','ADV2'); axis('tight')


res(1:4)=sum(sum(squeeze(RES(:,[1:309],[1:4]))))
cd /data1/project/gyre/SOSE_ITER60/temp_budget_scripts
klev=1
lim=0.01
figure(15)
clf
subplot(411)
pcolor(1e6*squeeze(double(RES(:,:,klev)'))); shading flat; colorbar; caxis([-1 1].*lim); colormap('redblue'); 
title(['Residual (MISMATCH) x 1e6   level=1   iter= ',num2str(it),'  res= ',num2str(res(1))])
subplot(412)
pcolor(1e6*squeeze(double(RES(:,:,2)'))); shading flat; colorbar; caxis([-1 1].*lim); 
title(['Residual (MISMATCH) x 1e6   level=2  ,  res= ',num2str(res(1,2))])
subplot(413)
pcolor(1e6*squeeze(double(RES(:,:,3)'))); shading flat; colorbar; caxis([-1 1].*lim); 
title(['Residual (MISMATCH) x 1e6   level=3  ,  res= ',num2str(res(1,3))])
subplot(414)
pcolor(1e6*squeeze(double(RES(:,:,4)'))); shading flat; colorbar; caxis([-1 1].*lim); 
title(['Residual (MISMATCH) x 1e6   level=4   ,  res= ',num2str(res(1,4))])

figure(16)
clf
subplot(411)
pcolor(-squeeze(double(TFLX(:,:,1)'))); shading flat; colorbar; caxis([-1 1].*2e-5); colormap('redblue'); 
subplot(412)
pcolor(-squeeze(double(Qsw(:,:,1)'))); shading flat; colorbar; caxis([-1 1].*2e-5); colormap('redblue'); 


%~~~~~~~~~~~~~~~~~~~~~
ADV1x=zeros(size(Ax)); ADV1y=ADV1x; ADV1z=ADV1x;ADV2x=zeros(size(Ax)); ADV2y=ADV2x; ADV2z=ADV2x;
ADV1x = -(RecipVol(ILON,ILAT,IDEP).*Ax);
ADV1y = -(RecipVol(ILON,ILAT,IDEP).*Ay);
ADV1z = -(RecipVol(ILON,ILAT,IDEP).*Az);
ADV2x = ((RecipVol(ILON,ILAT,IDEP).*THETA(ILON,ILAT,IDEP)).*DU);
ADV2y = ((RecipVol(ILON,ILAT,IDEP).*THETA(ILON,ILAT,IDEP)).*DV);
ADV2z = ((RecipVol(ILON,ILAT,IDEP).*THETA(ILON,ILAT,IDEP)).*DW);

mm=squeeze(hFacC(ILON,ILAT,IDEP)); mm=double(mm); 
figure(5)
clf
subplot(211)
plot(squeeze(sum(sum(mm.*ADV1x(ILON,ILAT,IDEP)))),'x-','linewidth',2,'color','k'); hold on; grid on 
plot(squeeze(sum(sum(mm.*ADV1y(ILON,ILAT,IDEP)))),'x-','linewidth',2,'color','r'); title(['iter =  ',num2str(it)])
plot(squeeze(sum(sum(mm.*ADV1z(ILON,ILAT,IDEP)))),'x-','linewidth',2,'color','b'); 
plot(squeeze(sum(sum(mm.*ADV1x(ILON,ILAT,IDEP))))+squeeze(sum(sum(mm.*ADV1y(ILON,ILAT,IDEP))))+...
       squeeze(sum(sum(mm.*ADV1z(ILON,ILAT,IDEP)))),'linewidth',2,'color','g'); legend('ADV1x','ADV1y','ADV1z','sum'); axis('tight')
subplot(212)
plot(squeeze(sum(sum(mm.*ADV2x(ILON,ILAT,IDEP)))),'x-','linewidth',2,'color','k'); hold on 
plot(squeeze(sum(sum(mm.*ADV2y(ILON,ILAT,IDEP)))),'x-','linewidth',2,'color','r'); grid on 
plot(squeeze(sum(sum(mm.*ADV2z(ILON,ILAT,IDEP)))),'x-','linewidth',2,'color','b');
% plot(squeeze(sum(sum(mm.*ADV2x(ILON,ILAT,IDEP))))+squeeze(sum(sum(mm.*ADV2y(ILON,ILAT,IDEP))))+...
 





