function data_out=analyze_thscan_new(scannums,plotflag)

if(nargin<2) 
    plotflag = 0;
end
%Edit detector channels here
XRFchan = 7;
imagechan = 23;
thetachan = 56;
twothetachan = 55;
gammachan = 54;
rdetchan = 53;

%Edit ROI here

%ROIXstart = 158; ROIXwidth = 200; ROIYstart = 158; ROIYwidth = 200; %111
%central box, large, NW1
%ROIXstart = 255; ROIXwidth = 150; ROIYstart = 320; ROIYwidth = 190; 
%bottom box
ROIXstart = 150; ROIXwidth = 250; ROIYstart = 100; ROIYwidth = 300; 


data_out.thvals = zeros(max(size(scannums)),1);

XRF1 = loadmda(['mda/26idbSOFT_' num2str(scannums(1),'%4.4d') '.mda'],XRFchan,0,0);
twotheta = loadmda(['mda/26idbSOFT_' num2str(scannums(1),'%4.4d') '.mda'],twothetachan,0,0);
gamma = loadmda(['mda/26idbSOFT_' num2str(scannums(1),'%4.4d') '.mda'],gammachan,0,0);
rdet = loadmda(['mda/26idbSOFT_' num2str(scannums(1),'%4.4d') '.mda'],rdetchan,0,0);
for ii =1:size(XRF1,1)
    for jj=1:size(XRF1,2)
        data_out.ii(ii).jj(jj).im = zeros(ROIYwidth+1,ROIXwidth+1);
        data_out.ii(ii).jj(jj).immaxrc = zeros(ROIYwidth+1,ROIXwidth+1);
    end
end
chi2 = zeros(20,20);
data_out.curve=zeros(2,max(size(scannums)));
data_out.Intensity = zeros(size(XRF1,1),size(XRF1,2));
data_out.XCent = zeros(size(XRF1,1),size(XRF1,2));
data_out.YCent = zeros(size(XRF1,1),size(XRF1,2));
data_out.thmax = zeros(size(XRF1,1),size(XRF1,2));
data_out.thcen = zeros(size(XRF1,1),size(XRF1,2));
data_out.ROIXstart = ROIXstart;
data_out.ROIYstart = ROIYstart;
data_out.ROIXwidth = ROIXwidth;
data_out.ROIYwidth = ROIYwidth;
data_out.twotheta = twotheta(1,1,1);
data_out.gamma = gamma(1,1,1);
data_out.rdet = rdet(1,1,1);


%%

%%
h=waitbar(0,'Loading multiple scans');
if(plotflag)
    figure(101);
    clf reset;
    imagesc(XRF1(1,:,3),XRF1(:,1,2),XRF1(:,:,1));colormap hot;axis square tight;
    pause(0.5)
end
for mm=1:max(size(scannums))
    waitbar(mm/max(size(scannums,2)));
    mdanum=scannums(mm);
    th_temp = loadmda(['mda/26idbSOFT_' num2str(mdanum,'%4.4d') '.mda'],thetachan,0,0);
    data_out.thvals(mm) = th_temp(1,1,1);
    data_out.curve(1,mm) = th_temp(1,1,1);
    data_out.scan(mm).scannum = mdanum;
    ccdnums = loadmda(['mda/26idbSOFT_' num2str(mdanum,'%4.4d') '.mda'],imagechan,0,0);
    XRF2 = loadmda(['mda/26idbSOFT_' num2str(mdanum,'%4.4d') '.mda'],XRFchan,0,0);
   
    %%{
    % lineup with fluorescence centroid
    y_mass = sum(XRF2(:,:,1),1);
    x_mass = sum(XRF2(:,:,1),2);
    yline = ([1:1:numel(y_mass)]);
    xline = transpose([1:1:numel(x_mass)]);
    shiftx = round(round(numel(xline)/2) - sum((x_mass).*xline)/sum((x_mass)));
    shifty = round(round(numel(yline)/2) - sum((y_mass).*yline)/sum((y_mass)));
    data_out.scan(mm).XRF = XRF2;  %keep the same scanning axes
    data_out.scan(mm).XRF(:,:,1) = circshift(XRF2(:,:,1),[shiftx,shifty]);
    data_out.scan(mm).ccdnums = circshift(ccdnums,[shiftx,shifty]);
    %}
   %{
    % line up with difference minimization
     for jj = 1:size(chi2,1)
        for kk = 1:size(chi2,2)
            chi2(jj,kk) = sum(sum((XRF1(:,:,1) - circshift(XRF2(:,:,1),[jj-10,kk-10])).^2));
        end
    end
    [xcen,ycen,intemp] = find(chi2==min(min(chi2)));
    
    data_out.scan(mm).XRF = circshift(XRF2,[xcen-10,ycen-10]);
    data_out.scan(mm).ccdnums = circshift(ccdnums,[xcen-10,ycen-10]);
    %}
    
    for ii = 1:size(XRF1,1)
        for jj= 1:size(XRF1,2)
            filename=['Images/' num2str(mdanum) '/scan_' num2str(mdanum) '_img_' num2str(data_out.scan(mm).ccdnums(ii,jj), '%6.6d') '.tif'];
            
            %ccd = double(imread(filename));
            ccd = double(imread(filename, 'PixelRegion', {[ROIYstart, ROIYstart+ROIYwidth],[ROIXstart, ROIXstart+ROIXwidth]})); 
            ccd1=ccd;
            ccd=ccd.*(ccd1>0);

            % hot pixel removal UNCOMMENT HERE
            %%{
            ccd1 = zeros(size(ccd,1),size(ccd,2),4);
            ccd1(:,:,1) = circshift(ccd,[0 1]);
            ccd1(:,:,2) = circshift(ccd,[1 0]);
            ccd1(:,:,3) = circshift(ccd,[0 -1]);
            ccd1(:,:,4) = circshift(ccd,[-1 0]);
            ccd2 = median(ccd1,3);
            ccdmask = ccd>(ccd2+10);   %CHANGE HOT PIXEL THRESHOLD HERE
            ccd = ccd.*(1-ccdmask)+ccd2.*ccdmask;
            %}
            
            %ROI1 = ccd(ROIYstart:(ROIYstart+ROIYwidth),ROIXstart:(ROIXstart+ROIXwidth));
            ROI1 = ccd;

            %display(size(data_out.ii(ii).jj(jj).im));display(size(ROI1));
            data_out.ii(ii).jj(jj).im = data_out.ii(ii).jj(jj).im + ROI1;
            data_out.ii(ii).jj(jj).rc(mm) = sum(ROI1(:));
            

            %identify and save image from maximum of rocking curve at each
            %position
            if mm==1
                data_out.ii(ii).jj(jj).immaxrc = ROI1;
                data_out.ii(ii).jj(jj).th = th_temp(mm);
            elseif sum(ROI1(:))> sum(data_out.ii(ii).jj(jj).immaxrc(:)) %replace immaxrc image only if current ROI integration is greater than previous
                data_out.ii(ii).jj(jj).immaxrc = ROI1;
                data_out.ii(ii).jj(jj).th = th_temp(mm);
            end 
                        
            %data_out.curve(2,mm) = data_out.curve(2,mm) +sum(sum(ccd(ROIYstart:(ROIYstart+ROIYwidth),ROIXstart:(ROIXstart+ROIXwidth))));
            data_out.curve(2,mm) = data_out.curve(2,mm) +sum(sum(ROI1));
        end
    end
    if(plotflag)
    figure(101);
    clf reset;
    imagesc(XRF1(1,:,3),XRF1(:,1,2),data_out.scan(mm).XRF(:,:,1));colormap hot;axis square tight;
    pause(0.5)
    end
end
close(h);

fluo = data_out.scan(1).XRF(:,:,1);
maxfluo = max(fluo(:));
cutoff = 0.3;
ind = find(fluo < cutoff*maxfluo);

h=waitbar(0,'Calculating centroids');
for ii =1:size(XRF1,1)
    waitbar(ii/size(XRF1,1));
    for jj=1:size(XRF1,2)
        for kk = 1:max(size(scannums))
            data_out.thcen(ii,jj) = data_out.thcen(ii,jj)+data_out.thvals(kk)*data_out.ii(ii).jj(jj).rc(kk);
        end
        data_out.thcen(ii,jj) = data_out.thcen(ii,jj)/sum(sum(data_out.ii(ii).jj(jj).rc));
        %ROI1 = data_out.ii(ii).jj(jj).immaxrc; %use to calculate pixel centroids from MAX
        ROI1 = data_out.ii(ii).jj(jj).im; %use to calculate pixel centroids from SUM
        data_out.Intensity(ii,jj)=sum(sum(ROI1));
        line1=sum(ROI1,1);  % horizontal
        line2=sum(ROI1,2);  % vertical
        data_out.thmax(ii,jj) = data_out.ii(ii).jj(jj).th;
        for kk=1:size(line1,2)
            data_out.XCent(ii,jj)=data_out.XCent(ii,jj)+kk*line1(kk)/data_out.Intensity(ii,jj);
        end
        data_out.XCent(ii,jj)=data_out.XCent(ii,jj)+ROIXstart;
        for kk=1:size(line2,1)
            data_out.YCent(ii,jj)=data_out.YCent(ii,jj)+kk*line2(kk)/data_out.Intensity(ii,jj);
        end
        data_out.YCent(ii,jj)=data_out.YCent(ii,jj)+ROIYstart;
    end
end
close(h);

%data_out.XCent(ind) = NaN;
%data_out.YCent(ind) = NaN;
%data_out.thmax(ind) = NaN;

  if(plotflag)  
   figure(102);plot(data_out.curve(1,:),data_out.curve(2,:),'o','MarkerSize',10,'MarkerFaceColor','k');axis square;   
   figure(103);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data_out.XCent);axis image;colormap hot; shading interp; title('ROI X Centroid');set(gca, 'YDir', 'normal');colorbar;
   figure(104);clf;hfig=imagesc(XRF1(1,:,3),XRF1(:,1,2),data_out.scan(1).XRF(:,:,1));axis image;colormap hot; shading interp; title('Zn XRF');set(gca, 'YDir', 'normal');colorbar;
   % Add CCD image numbers and pass structure to data cursor
            pass2click.data = data_out;
            pass2click.xaxis = [min(XRF1(1,:,3)) max(XRF1(1,:,3)) size(XRF1(1,:,3),2)];
            pass2click.yaxis = [min(XRF1(:,1,2)) max(XRF1(:,1,2)) size(XRF1(:,1,2),1)];
            set(hfig, 'UserData', pass2click);
            datacursormode on;
            dcm_obj = datacursormode(gcf);
            set(dcm_obj, 'DisplayStyle', 'window');
           set(dcm_obj, 'UpdateFcn', @click4rock);
   figure(105);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data_out.YCent);axis image;colormap hot; shading interp; title('ROI Y Centroid');set(gca, 'YDir', 'normal');colorbar;
   figure(106);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data_out.thmax);axis image;colormap hot; shading interp; title('Maximum Theta');set(gca, 'YDir', 'normal');colorbar;
   figure(107);clf;imagesc(XRF1(1,:,3),XRF1(:,1,2),data_out.thcen);axis image;colormap hot; shading interp; title('Centroid Theta');set(gca, 'YDir', 'normal');colorbar;

  end
end




