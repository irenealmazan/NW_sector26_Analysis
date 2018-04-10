function data_out = get_twoth(data,plotflag)
% function to interpolate twotheta and gamma per pixel from centroids
Ekev = 10.0; % beam energy in kev
kb = 2*pi*Ekev/12.39842; % beam momentum 1/A
lambda = 2*pi/kb; % wavelength A

qmat = qmat_mpx3(data.twotheta,data.rdet*1000,data.gamma,0,0,0);
data_out = zeros(size(data.XCent,1),size(data.XCent,2),3);

for ii = 1:size(data.XCent,1)
    for jj = 1:size(data.XCent,2)
        data_out(ii,jj,1) = interp2(qmat.twotheta,data.XCent(ii,jj),data.YCent(ii,jj));
        data_out(ii,jj,3) = interp2(qmat.gamma,data.XCent(ii,jj),data.YCent(ii,jj));
    end
end
data_out(:,:,2) = lambda./(2*sind(data_out(:,:,1)/2));

XRF1 = data.scan(1).XRF;
if(plotflag)
%figure(121);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data_out(:,:,1));axis image;colormap hot; shading interp; title('Two theta');set(gca, 'YDir', 'normal');colorbar;
%figure(122);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data_out(:,:,2));axis image;colormap hot; shading interp; title('D spacing');set(gca, 'YDir', 'normal');colorbar;
%figure(123);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data_out(:,:,3));axis image;colormap hot; shading interp; title('Gamma');set(gca, 'YDir', 'normal');colorbar;
fluo = data.scan(1).XRF(:,:,1);
maxfluo = max(fluo(:));
cutoff = 0.3;
ind = find(fluo < cutoff*maxfluo);
tempm=data_out(:,:,1);tempm(ind)=NaN;
figure(124);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),tempm);axis image;colormap jet; shading interp; title('Two theta');set(gca, 'YDir', 'normal');colorbar;
tempm=data_out(:,:,2);tempm(ind)=NaN;
figure(125);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),tempm);axis image;colormap jet; shading interp; title('D spacing');set(gca, 'YDir', 'normal');colorbar;
tempm=data_out(:,:,3);tempm(ind)=NaN;
figure(126);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),tempm);axis image;colormap jet; shading interp; title('Gamma');set(gca, 'YDir', 'normal');colorbar;
tempm=data.thcen;tempm(ind)=NaN;
figure(127);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),tempm);axis image;colormap jet; shading interp; title('Theta');set(gca, 'YDir', 'normal');colorbar;
figure(128);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),fluo);axis image;colormap jet; shading interp; title('XRF');set(gca, 'YDir', 'normal');colorbar;

end
