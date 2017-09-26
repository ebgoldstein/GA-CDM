clear all 
close all
clc

FF1=load('rawdensecloud/FF1dense.txt');    %4/26
FF2=load('rawdensecloud/FF2dense.txt');    %9/12
FF3=load('rawdensecloud/FF3dense.txt');    %11/

elev1=FF1(:,3);
elev2=FF2(:,3);
elev3=FF3(:,3);

grey1=((FF1(:,4)+FF1(:,5)+FF1(:,6))/3)/256;
grey2=((FF2(:,4)+FF2(:,5)+FF2(:,6))/3)/256;
grey3=((FF3(:,4)+FF3(:,5)+FF3(:,6))/3)/256;

minx=min(FF1(:,1));
miny=min(FF1(:,2));

FF1(:,1)=FF1(:,1)-minx;
FF2(:,1)=FF2(:,1)-minx;
FF3(:,1)=FF3(:,1)-minx;

FF1(:,2)=FF1(:,2)-miny;
FF2(:,2)=FF2(:,2)-miny;
FF3(:,2)=FF3(:,2)-miny;


% https://www.mathworks.com/matlabcentral/answers/9554-rotate-clouds-of-points
% rotate by 30 clockwise around (0,0)
% http://en.wikipedia.org/wiki/Rotation_matrix
rotAngle = 17; %from google earth
rotRad = rotAngle *(pi/180);
FF1xRot     = FF1(:,1)*cos(rotRad) - FF1(:,2)*sin(rotRad);
FF1yRot     = FF1(:,1)*sin(rotRad) + FF1(:,2)*cos(rotRad);
%plot(FF1xRot,FF1yRot,'k*');
%hold on
FF2xRot     = FF2(:,1)*cos(rotRad) - FF2(:,2)*sin(rotRad);
FF2yRot     = FF2(:,1)*sin(rotRad) + FF2(:,2)*cos(rotRad);
%plot(FF2xRot,FF2yRot,'r*');
FF3xRot     = FF3(:,1)*cos(rotRad) - FF3(:,2)*sin(rotRad);
FF3yRot     = FF3(:,1)*sin(rotRad) + FF3(:,2)*cos(rotRad);


% https://www.mathworks.com/help/matlab/ref/griddata.html
%Grid from
%y axis of 400 to 500
ygrid=400:1:500;
%x axis from 0 to 150
xgrid=0:1:150;

%remove data that is not within boundaris of the griddata call
xindone=find(FF1xRot>=-10 & FF1xRot<=160);
FF1xRot=FF1xRot(xindone);
FF1yRot=FF1yRot(xindone);
elev1=elev1(xindone);
grey1=grey1(xindone);
yindone=find(FF1yRot>=390 & FF1yRot<=510);
FF1xRot=FF1xRot(yindone);
FF1yRot=FF1yRot(yindone);
elev1=elev1(yindone);
grey1=grey1(yindone);

xindtwo=find(FF2xRot>=-10 & FF2xRot<=160);
FF2xRot=FF2xRot(xindtwo);
FF2yRot=FF2yRot(xindtwo);
elev2=elev2(xindtwo);
grey2=grey2(xindtwo);
yindtwo=find(FF2yRot>=390 & FF2yRot<=510);
FF2xRot=FF2xRot(yindtwo);
FF2yRot=FF2yRot(yindtwo);
elev2=elev2(yindtwo);
grey2=grey2(yindtwo);

xindthree=find(FF3xRot>=-10 & FF3xRot<=160);
FF3xRot=FF3xRot(xindthree);
FF3yRot=FF3yRot(xindthree);
elev3=elev3(xindthree);
grey3=grey3(xindthree);
yindthree=find(FF3yRot>=390 & FF3yRot<=510);
FF3xRot=FF3xRot(yindthree);
FF3yRot=FF3yRot(yindthree);
elev3=elev3(yindthree);
grey3=grey3(yindthree);

%topo griddata
[xq,yq]= meshgrid(xgrid,ygrid);
FF1regrid = griddata(FF1xRot,FF1yRot ,elev1,xq,yq);
FF2regrid = griddata(FF2xRot,FF2yRot ,elev2,xq,yq);
FF3regrid = griddata(FF3xRot,FF3yRot ,elev3,xq,yq);

%veg griddata
FF1Vregrid = griddata(FF1xRot,FF1yRot ,grey1,xq,yq);
FF2Vregrid = griddata(FF2xRot,FF2yRot ,grey2,xq,yq);
FF3Vregrid = griddata(FF3xRot,FF3yRot ,grey3,xq,yq);

%%%%CHOP the grid at the CDM boundaries
CDMTgrid1=FF1regrid(1:64,27:150);
CDMTgrid2=FF2regrid(1:64,27:150);
CDMTgrid3=FF3regrid(1:64,27:150);

CDMVgrid1=FF1Vregrid(1:64,27:150);
CDMVgrid2=FF2Vregrid(1:64,27:150);
CDMVgrid3=FF3Vregrid(1:64,27:150);

%impose beach slope of 4 degrees at cell 89
beachstarts=89;
beachslope=4;
heighreduction=-tand(beachslope);
beachface=[2.6:heighreduction:0 0];

b=repmat(beachface,64,1);

CDMgrid1withbeach=[CDMTgrid1(:,1:beachstarts) b];
CDMgrid2withbeach=[CDMTgrid2(:,1:beachstarts) b];
CDMgrid3withbeach=[CDMTgrid3(:,1:beachstarts) b];

DEMdiff=CDMgrid3withbeach-CDMgrid1withbeach;

vegbeachface=beachface*0;
bveg=repmat(vegbeachface,64,1);
CDMVEGgrid1withbeach=[CDMVgrid1(:,1:beachstarts) bveg];
CDMVEGgrid1withbeach(CDMVEGgrid1withbeach==nan) = 0;
CDMVEGgrid2withbeach=[CDMVgrid2(:,1:beachstarts) bveg];
CDMVEGgrid2withbeach(CDMVEGgrid2withbeach==nan) = 0;
CDMVEGgrid3withbeach=[CDMVgrid3(:,1:beachstarts) bveg];
CDMVEGgrid3withbeach(CDMVEGgrid3withbeach==nan) = 0;

%make the April Veg grids perent cover using a rule that i found.
CDMVEGgrid1withbeach=(-(1/.45)*(CDMVEGgrid1withbeach))+1;
CDMVEGgrid1withbeach(CDMVEGgrid1withbeach>=1)=0;
CDMVEGgrid1withbeach(CDMVEGgrid1withbeach<0)=0;
CDMVEGgrid1withbeach(isnan(CDMVEGgrid1withbeach))=0;

%make the Sept Veg grids perent cover using a rule that i found.
CDMVEGgrid2withbeach=(-(1/.45)*(CDMVEGgrid2withbeach))+1.45;
CDMVEGgrid2withbeach(CDMVEGgrid2withbeach>=1)=0;
CDMVEGgrid2withbeach(CDMVEGgrid2withbeach<0)=0;
CDMVEGgrid2withbeach(isnan(CDMVEGgrid2withbeach))=0;


%%%%%%%%%%%
%fix a hole
CDMVEGgrid1withbeach(40,64)=CDMVEGgrid1withbeach(40,63);


%%%%%%%%%%%
figure

subplot(2,2,1),imagesc(CDMVEGgrid1withbeach)
title('Veg 4/26/2016')
xlabel('m')
ylabel('m')
caxis([0 1])
colorbar

subplot(2,2,2),imagesc(CDMgrid1withbeach)
title('4/26/2016')
xlabel('m')
ylabel('m')
colorbar
caxis([1 5])

subplot(2,2,3),imagesc(CDMVEGgrid3withbeach)
title('Veg 11/2016')
xlabel('m')
ylabel('m')
colorbar
caxis([0 1])

subplot(2,2,4),imagesc(CDMgrid3withbeach)
title('11/2016')
xlabel('m')
ylabel('m')
colorbar
caxis([1 5])

%%%%Make it ready for the model
Apr=flipud(CDMgrid1withbeach');
%Apr=round(Apr*10000)/10000;
dlmwrite('April.dat',Apr,'delimiter',' ')

AprV=flipud(CDMVEGgrid1withbeach');
%AprV=round(AprV*10000)/10000;
dlmwrite('AprilV.dat',AprV,'delimiter',' ')

%make a dummy for CDM veg Y
AprilVy=AprV*0;
dlmwrite('AprilVy.dat',AprilVy,'delimiter',' ')

%%%%
Sept=flipud(CDMgrid2withbeach');
%Sept=round(Sept*10000)/10000;
dlmwrite('September.dat',Sept,'delimiter',' ')

SeptV=flipud(CDMVEGgrid2withbeach');
%SeptV=round(SeptV*10000)/10000;
dlmwrite('SeptemberV.dat',SeptV,'delimiter',' ')
%%%%

Nov=flipud(CDMgrid3withbeach');
%Nov=round(Nov*10000)/10000;
dlmwrite('November.dat',Nov,'delimiter',' ')

NovV=flipud(CDMVEGgrid2withbeach');
%NovV=round(NovV*10000)/10000;
dlmwrite('NovemberV.dat',NovV,'delimiter',' ')


%%%%
dlmwrite('Change.dat',(Nov-Apr),'delimiter',' ')
dlmwrite('ChangeV.dat',(NovV-AprV),'delimiter',' ')

%%%%%
%%%%%%%%%%%
figure

subplot(3,1,1),imagesc(CDMgrid1withbeach)
%title('4/26/2016')
colormap('gray')
xlabel('m')
ylabel('m')
colorbar
axis xy
caxis([1 5])

subplot(3,1,2),imagesc(CDMgrid3withbeach)
%title('9/12/2016')
colormap('gray')
xlabel('m')
ylabel('m')
colorbar
axis xy
caxis([1 5])

subplot(3,1,3),imagesc(CDMgrid3withbeach-CDMgrid1withbeach)
%title('Change')
colormap('jet')
xlabel('m')
ylabel('m')
colorbar
axis xy
caxis([-0.5 0.5])


% dlmwrite('April.dat',CDMgrid1withbeach,'delimiter',' ')
% dlmwrite('September.dat',CDMgrid2withbeach,'delimiter',' ')



