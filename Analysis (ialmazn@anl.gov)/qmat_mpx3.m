function dataout = qmat_mpx3(Tthdet, Rdet, Gamdet, Xdet,Ydet,plotflag)

%---------------------------------------------
%PARAMETERS FOR MOMENTUM SPACE BINNING MEDIPIX IMAGES
% All vector notation in beam convention [+outboard +up +downstream]

Ekev = 10.0; % beam energy in kev
kb = 2*pi*Ekev/12.39842; % beam momentum 1/A
hpxsz = 55; %microns (MEDIPIX3)
vpxsz = 55;

if(nargin<6) 
    plotflag=0; % output images
end


if(nargin<5) 
    Xdet = 0;
    Ydet = 0;
end

if(nargin<3)
    Tthdet=22; % **Cylindrical** two theta angle to the center of the detector (deg)
    Rdet=850000;  % **Cylindrical** radius from sample to the center detector (um)
    Gamdet=0;   % **Cylindrical** out of diffraction plane angle to the center of the detector (deg)
end

scrsz = get(0,'ScreenSize');  % Plot size parameters
imd = scrsz(3)/5;
imb = scrsz(4)-imd;

vNdet = 516; 
vNcen = round(vNdet/2);
vimaxis = (1-vNcen):(vNdet-vNcen);
hNdet = 516; 
hNcen = round(hNdet/2);
himaxis = (1-hNcen):(hNdet-hNcen);

%Detector center and pixel plane vectors
tthd = pi * Tthdet / 180;
gamd = pi * Gamdet / 180;
Detcen = [Rdet*sin(tthd)*cos(gamd) Rdet*sin(gamd) Rdet*cos(tthd)*cos(gamd)];
Detunit = Detcen./(sqrt(dot(Detcen,Detcen)));
jvec = cross([0 1 0],Detunit); % horizontal pixels unit vector (low pixel number, low tth)
jvec = jvec./(sqrt(dot(jvec,jvec)));
ivec = cross(jvec,Detunit); % vertical pixels unit vector (low pixel number, more up)

kfmat=zeros(vNdet,hNdet,4);
kimat=zeros(vNdet,hNdet,4);
qmat=zeros(vNdet,hNdet,4);

% Assign magnitude and vector q to active detector pixels within q and gamma range
%display('assigning q vectors...');
%Generate pixel arrays (um)
pix_x = repmat((himaxis.*(hpxsz)),vNdet,1)+Xdet;
pix_y = repmat((vimaxis.*(vpxsz))',1,hNdet)+Ydet;
%Calculate real space pixel positions in mm as a vector sum
kfmat(:,:,1) = pix_x(:,:).*jvec(1)+pix_y(:,:).*ivec(1)+Detcen(1);
kfmat(:,:,2) = pix_x(:,:).*jvec(2)+pix_y(:,:).*ivec(2)+Detcen(2);
kfmat(:,:,3) = pix_x(:,:).*jvec(3)+pix_y(:,:).*ivec(3)+Detcen(3);
% if(plotflag)
%        figure(11);clf;set(11,'Position',[imd imb imd imd]);clf;imagesc(kfmat(:,:,1));axis square;
%        title('X real space position'); pause(0.5);
%        figure(16);set(16,'Position',[imd imb imd imd]);clf;imagesc(kfmat(:,:,2));axis square;
%        title('Y real space position'); 
%        figure(17);set(17,'Position',[imd imb imd imd]);clf;imagesc(kfmat(:,:,3));axis square;
%        title('Z real space position');
% end
kfmat(:,:,4) = sqrt(kfmat(:,:,1).^2+kfmat(:,:,2).^2+kfmat(:,:,3).^2);

%Assign kfinal unit vector to each pixel based on position
kfmat(:,:,1) = kfmat(:,:,1)./kfmat(:,:,4);
kfmat(:,:,2) = kfmat(:,:,2)./kfmat(:,:,4);
kfmat(:,:,3) = kfmat(:,:,3)./kfmat(:,:,4);

kimat(:,:,1)=zeros(vNdet,hNdet);
kimat(:,:,2)=zeros(vNdet,hNdet);
kimat(:,:,3)=ones(vNdet,hNdet);
kimat(:,:,4)=ones(vNdet,hNdet);

%Calculate q vectors
qmat(:,:,1) = kb*(kfmat(:,:,1)-kimat(:,:,1));
qmat(:,:,2) = kb*(kfmat(:,:,2)-kimat(:,:,2));
qmat(:,:,3) = kb*(kfmat(:,:,3)-kimat(:,:,3));
qmat(:,:,4) = sqrt(qmat(:,:,1).^2+qmat(:,:,2).^2+qmat(:,:,3).^2);
if(plotflag)
       figure(118);set(118,'Position',[imd imb imd imd]);clf;imagesc(acosd(kfmat(:,:,3)));axis image;
       title('Two theta per pixel map');colorbar;
       figure(119);set(119,'Position',[imd imb imd imd]);clf;imagesc(90-acosd(kfmat(:,:,2)));axis image;
       title('Gamma per pixel map');colorbar
       figure(120);set(120,'Position',[imd imb imd imd]);clf;imagesc(qmat(:,:,4));axis image;
       title('Momentum transfer (A^-1) per pixel map');colorbar
       figure(121);set(120,'Position',[imd imb imd imd]);clf;imagesc(2*pi./qmat(:,:,4));axis image;
       title('D spacing (A) per pixel map');colorbar
end

if (plotflag==1)
%     figure(11);set(11,'Position',[imd imb imd imd]);clf;imagesc(imq);axis square;
%     title('Magnitude q per pixel'); pause(0.5);
%     figure(16);set(16,'Position',[imd imb imd imd]);clf;imagesc(imqx);axis square;
%     title('Magnitude qx per pixel'); 
%     figure(17);set(17,'Position',[imd imb imd imd]);clf;imagesc(imqy);axis square;
%     title('Magnitude qy per pixel'); 
%     figure(18);set(18,'Position',[imd imb imd imd]);clf;imagesc(imqz);axis square;
%     title('Magnitude qz per pixel'); 

end

dataout.twotheta = acosd(kfmat(:,:,3));
dataout.qmat = qmat(:,:,4);
dataout.gamma = 90-acosd(kfmat(:,:,2));
end