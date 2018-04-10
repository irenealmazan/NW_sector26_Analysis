function display_rock(data)
XRF1 = data.scan(1).XRF;
   figure(102);plot(data.curve(1,:),data.curve(2,:),'o','MarkerSize',10,'MarkerFaceColor','k');axis square;   
   figure(103);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data.XCent);axis image;colormap hot; shading interp; title('ROI X Centroid');set(gca, 'YDir', 'normal');colorbar;
   figure(104);clf reset;hfig=imagesc(XRF1(1,:,3),XRF1(:,1,2),data.scan(1).XRF(:,:,1));axis image;colormap hot; shading interp; title('Zn XRF');set(gca, 'YDir', 'normal');colorbar;
   % Add CCD image numbers and pass structure to data cursor
            pass2click.data = data;
            pass2click.xaxis = [min(XRF1(1,:,3)) max(XRF1(1,:,3)) size(XRF1(1,:,3),2)];
            pass2click.yaxis = [min(XRF1(:,1,2)) max(XRF1(:,1,2)) size(XRF1(:,1,2),1)];
            set(hfig, 'UserData', pass2click);
            datacursormode on;
            dcm_obj = datacursormode(gcf);
            set(dcm_obj, 'DisplayStyle', 'window');
           set(dcm_obj, 'UpdateFcn', @click4rock);
   figure(105);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data.YCent);axis image;colormap hot; shading interp; title('ROI Y Centroid');set(gca, 'YDir', 'normal');colorbar;
   figure(106);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data.thcen);axis image;colormap hot; shading interp; title('Centroid Theta');set(gca, 'YDir', 'normal');colorbar;
for kk = 1:size(data.curve,2);
    figure(107);clf reset
    imagesc(XRF1(1,:,3),XRF1(:,1,2),data.scan(kk).XRF(:,:,1));
    pause(1);
end

end