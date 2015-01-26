% cd /data1/project/gyre/SOSE_ITER60/salt_budget_scripts
% New_Salt_Budget

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
display(' READ IN SALT BUDGET TERMS ')
fname_int=[ 'ADVx_SL ' 'ADVy_SL ' 'ADVr_SL '  'DFxE_SL ' 'DFyE_SL ' 'DFrE_SL ' 'DFrI_SL ' 'KPPg_SL ' 'TOTSTEND'];
fname_stt =[ 'SSALT   '  'UVEL    ' 'VVEL    ' 'WVEL    '];
fname_frc=['SFLUX'];


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
 tm1 = Q(ILON+1,ILAT,IDEP,1);% 'ADVx_SL '
 tm2 = Q(ILON  ,ILAT,IDEP,1);% 'ADVx_SL '
Ax  = (tm1 - tm2);
 tm1 = Q(ILON,ILAT+1,IDEP,2);% 'ADVy_SL '
 tm2 = Q(ILON,ILAT  ,IDEP,2);% 'ADVy_SL '
Ay  = (tm1 - tm2);
 tm1 = Q(ILON,ILAT,IDEP  ,3);% 'ADVr_SL '
 tm2 = Q(ILON,ILAT,IDEP+1,3);% 'ADVr_SL '
Az = (tm1 - tm2);

% DIFFUSIVE TERMS -- zonal, then meridional, then vertical
 tm1 = Q(ILON+1,ILAT,IDEP,4);% 'DFxE_SL '
 tm2 = Q(ILON  ,ILAT,IDEP,4);% 'DFxE_SL '
Dx  = (tm1 - tm2);
 tm1 = Q(ILON,ILAT+1,IDEP,5);% 'DFyE_SL '
 tm2 = Q(ILON,ILAT  ,IDEP,5);% 'DFyE_SL '
Dy  = (tm1 - tm2);
 tm1 = Q(ILON,ILAT,IDEP  ,6);% 'DFrE_SL '
 tm2 = Q(ILON,ILAT,IDEP+1,6);% 'DFrE_SL '
Dze = (tm1 - tm2);%explicit vertical diff
 tm1 = Q(ILON,ILAT,IDEP  ,7);% 'DFrI_SL '
 tm2 = Q(ILON,ILAT,IDEP+1,7);% 'DFrI_SL '
Dzi = (tm1 - tm2);%implicit vertical diff
 tm1 = Q(ILON,ILAT,IDEP  ,8);% 'KPPg_SL '
 tm2 = Q(ILON,ILAT,IDEP+1,8);% 'KPPg_SL '
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

SALT = Q(ILON,ILAT,IDEP,1);%'SALT   '
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


%  oceSflux 	net surface Salt flux into the ocean (+=down), >0 increases salinity  diagUnits = 'g/m^2/s  
clear Q
Q = zeros(NX,NY,'single');
cd  /data2/project/gyre-data/SOSE_ITER60/TFORCING
Q=rdmds('SFLUX.0000000060','rec',it);
fac = 1./(1035*DRF(1)); 
mask=squeeze(hFacC(ILON,ILAT,1));
mask(find(hFacC(ILON,ILAT,1)~=1))=0;                            % 1 / (rho_nil*DRF(1) );  [psu/sec]
SFLX =  1*fac.*Q(ILON,ILAT);     % Q is surfaceForcingS *rUnit2mass  % gS is +surfaceForcingS/drf(1)
SFLX=SFLX.*mask;
clear Q 



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
%  fVerT = KPPg_SL*maskUp*rhoFacF + DFrE_SL*maskUp + ADVr_SL;
%  fMer  = DFyE_SL*rhoFacF + ADVy_SL
%  fZon  = DFxE_SL*rhoFacF + ADVx_SL
% NOW GET TENDENCY FROM SUM OF RHS

gTracer = zeros(size(Ax));
gTracer = -RecipVol(ILON,ILAT,IDEP)   ...
          .*( (Ax+Dx+Ay+Dy+Az+Dze+KP) ... %
               + Dzi                );  ...
%              - SALT.*(DU+DV+DW)    );

gTracer(:,:,1) = gTracer(:,:,1) + SFLX(ILON,ILAT) - SALT(ILON,ILAT,1).*(DU(ILON,ILAT,1)+DV(ILON,ILAT,1)+DW(ILON,ILAT,1)); 
RES =    TT - gTracer;   


ADV1 =zeros(size(Ax)); ADV2=zeros(size(Ax)); ADV=zeros(size(Ax)); DIFE=ADV; DIFKP=ADV; DIFI=ADV; RES=ADV;
RecVol=RecipVol(ILON,ILAT,IDEP);
ADV1x=zeros(size(Ax));  ADVx=ADV1x; ADV1y=zeros(size(Ay));  ADVy=ADV1y; ADV1z=zeros(size(Az));  ADVz=ADV1z; 
ADV2x=zeros(NX,NY);  ADV2y=ADV2x;  ADV2z=ADV2x;

ADV1 = -RecVol.* (Ax+Ay+Az);     %   U dot grad T  = U \cdot \nabla T
ADV2 =  RecVol(ILON,ILAT,1).*(SALT(ILON,ILAT,1).*(DU(ILON,ILAT,1)+DV(ILON,ILAT,1)+DW(ILON,ILAT,1)));   %   T times div U = T \nabla \cdot U
ADV=ADV1;
ADV(:,:,1)=ADV1(:,:,1)+ADV2(:,:,1);

ADV1x = -RecVol.* Ax;     %   U dot grad T  = U \cdot \nabla T
ADV2x =  RecVol(ILON,ILAT,1).*(SALT(ILON,ILAT,1).*DU(ILON,ILAT,1));   
ADVx=ADV1x; ADVx(:,:,1)=ADV1x(:,:,1)+ADV2x(:,:,1);

ADV1y = -RecVol.* Ay;     
ADV2y =  RecVol(ILON,ILAT,1).*(SALT(ILON,ILAT,1).*DV(ILON,ILAT,1));   
ADVy=ADV1y; ADVy(:,:,1)=ADV1y(:,:,1)+ADV2y(:,:,1);

ADV1z = -RecVol.* Az;     
ADV2z =  RecVol(ILON,ILAT,1).*(SALT(ILON,ILAT,1).*DW(ILON,ILAT,1));   
ADVz=ADV1z; ADVz(:,:,1)=ADV1z(:,:,1)+ADV2z(:,:,1);


DIFE =  -RecVol.*(Dx+Dy+Dze); %  Sum of diffusive terms ~ \nabla \kappa \nabla T
DIFKP = -RecVol.* KP;
DIFI =  -RecVol.*Dzi;

if(if_sav==1)
cd    /data2/project/gyre-data/SOSE_ITER60/Sbudget_NEWEST
for isav=7:9
[isav]
  if(isav==1); fname='TEND_S' ; end
  if(isav==2); fname='ADV_S' ; end
  if(isav==3); fname='DIFE_S'  ; end  
  if(isav==4); fname='DIFI_S'  ; end  
  if(isav==5); fname='DIFKP_S'  ; end                 
  if(isav==6); fname='RES_S'; end      
  if(isav==7); fname='ADV_Sx' ; end
  if(isav==8); fname='ADV_Sy' ; end
  if(isav==9); fname='ADV_Sz' ; end
          
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
  if(isav==6); fwrite(fidd,RES,'single') ; end
  if(isav==7); fwrite(fidd,ADVx,'single') ; end
  if(isav==8); fwrite(fidd,ADVy,'single') ; end
  if(isav==9); fwrite(fidd,ADVz,'single') ; end
       
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

                              % is almost zero below top couple layers 
DIV_FLX=zeros(size(Ax));
DIV_FLX=DIFE+DIFI+DIFKP;
DIV_FLX(:,:,1)=DIFE(:,:,1)+DIFI(:,:,1)+DIFKP(:,:,1)+SFLX(ILON,ILAT);

IIXX=[1:2159];
JJYY=[10:300];
KKZZ=IDEP;
% mmm=squeeze(hFacC(IIXX,JJYY,KKZZ)); mmm=double(mmm); 

figure(1)
clf
subplot(211)
plot(squeeze(sum(sum(squeeze(double(TT(IIXX,JJYY,KKZZ)))))),'ro-'); hold on
plot(squeeze(sum(sum(squeeze(ADV(IIXX,JJYY,KKZZ))))),'x-','linewidth',2,'color','k'); 
plot(squeeze(sum(sum(squeeze(DIV_FLX(IIXX,JJYY,KKZZ))))),'bo-'); grid on; title(['SALT budget   iter= ',num2str(it)])
plot(squeeze(sum(sum(squeeze(RES(IIXX,JJYY,KKZZ))))),'o--','linewidth',2,'color','g'); grid on ; axis('tight')  ;
legend('S tend','S adv','S form','S sum')
%~~~~~~~~~~~~~~~~~

subplot(212)
plot(squeeze(sum(sum(squeeze(RES(IIXX,JJYY,KKZZ))))),'o-'); grid on ;  hold on;  title('sum')
plot(squeeze(sum(sum(squeeze(ADV1(IIXX,JJYY,KKZZ))))),'x-','linewidth',2,'color','g'); 
% plot(squeeze(sum(sum(squeeze(ADV2(IIXX,JJYY,KKZZ))))),'x--','linewidth',2,'color','r'); legend('SUM of the budget','ADV1','ADV2'); axis('tight')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res(1:4)=sum(sum(squeeze(RES(:,:,[1:4]))))
cd /data1/project/gyre/SOSE_ITER60/salt_budget_scripts
klev=1
lim=1e-10
figure(15)
clf
subplot(511)
pcolor(1e6*squeeze(double(RES(:,:,klev)'))); shading flat; colorbar; caxis([-1 1].*lim); colormap('redblue'); 
title(['Residual (MISMATCH) x 1e6   level=1   iter= ',num2str(it),'  res= ',num2str(res(1))])
subplot(512)
pcolor(1e6*squeeze(double(RES(:,:,2)'))); shading flat; colorbar; caxis([-1 1].*lim); 
title(['Residual (MISMATCH) x 1e6   level=2  ,  res= ',num2str(res(1,2))])
subplot(513)
pcolor(1e6*squeeze(double(RES(:,:,3)'))); shading flat; colorbar; caxis([-1 1].*lim); 
title(['Residual (MISMATCH) x 1e6   level=3  ,  res= ',num2str(res(1,3))])


subplot(514)
pcolor(-squeeze(double(ADV2(:,:,1))')); shading flat; colorbar; caxis([-1 1].*5e-6); 
title(['ADV2 in layer 1'])
subplot(515)
pcolor(-squeeze(double(SFLX'))); shading flat; colorbar; caxis([-1 1].*1e-6); 
title(['SFLX'])


%~~~~~~~~~~~~~~~~~~~~~
ADV1x=zeros(size(Ax)); ADV1y=ADV1x; ADV1z=ADV1x;ADV2x=zeros(size(Ax)); ADV2y=ADV2x; ADV2z=ADV2x;
ADV1x = -(RecipVol(ILON,ILAT,IDEP).*Ax);
ADV1y = -(RecipVol(ILON,ILAT,IDEP).*Ay);
ADV1z = -(RecipVol(ILON,ILAT,IDEP).*Az);
ADV2x = ((RecipVol(ILON,ILAT,IDEP).*SALT(ILON,ILAT,IDEP)).*DU);
ADV2y = ((RecipVol(ILON,ILAT,IDEP).*SALT(ILON,ILAT,IDEP)).*DV);
ADV2z = ((RecipVol(ILON,ILAT,IDEP).*SALT(ILON,ILAT,IDEP)).*DW);

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
plot(squeeze(sum(sum(mm.*ADV2x(ILON,ILAT,IDEP))))+squeeze(sum(sum(mm.*ADV2y(ILON,ILAT,IDEP))))+...
       squeeze(sum(sum(mm.*ADV2z(ILON,ILAT,IDEP)))),'linewidth',2,'color','g');  legend('ADV2x','ADV2y','ADV2z','sum'); axis('tight')
 
