% This scripts takes a dataout structure (generated by
% analyze_thscan_new.m) and displays:
% - the diffraction intensity maps for each angle of the rocking curve

% addpaths
addpath(genpath('Analysis (ialmazn@anl.gov)/'));

%% Diffraction intensity map:

% number of angles in the rocking curve scan
numthvals = size(data_rock.thvals,1);

% spatial coordinates (in micrometers) corresponding to the value of the
% hybridX,hybridY:
Xgrid = data_rock.scan(1).XRF(1,:,3);
Ygrid = data_rock.scan(1).XRF(:,1,2);

% Fluorescence map:
fluo = data_rock.scan(1).XRF(:,:,1);

figure(301);
clf;

for mmm = 1:numthvals
    for iii = 1:size(data_rock.ii,2)
        for jjj = 1:size(data_rock.ii(iii).jj,2)
            data_rock.thval(mmm).diffint(iii,jjj) = data_rock.ii(iii).jj(jjj).rc(mmm);                        
        end
    end
    
    subplot(ceil(sqrt(numthvals)),ceil(sqrt(numthvals)),mmm);
    imagesc(Xgrid,(Ygrid),(log10(data_rock.thval(mmm).diffint))); % the flipud is required to match the python figures
    set(gca, 'YDir', 'normal');
    axis image;
    colormap jet;
    colorbar;
    xlabel(['X(\mu m)']);
    ylabel(['Y(\mu m)']);
    title(['log(diff) \theta = ' num2str(data_rock.thvals(mmm))]);
end

figure(302);
clf;
imagesc(Xgrid,(Ygrid),(log10(fluo))); % the flipud is required to match the python figures
set(gca, 'YDir', 'normal');
axis image;
colormap jet;
colorbar;
xlabel(['X(\mu m)']);
ylabel(['Y(\mu m)']);
title(['log(fluo) ']);


%% Plot the D spacing, the theta angle of the maximum of the rocking curve and the gamma for each pixel and the contours

% data_out(:,:,1) is the two theta value per pixel
% data_out(:,:,2) is the D spacing value per pixel (using Bragg's law)
% data_out(:,:,3) is the gamma value per pixel
data_out = get_twoth(data_rock,1);

figure(124);
hold on;
[C,h] = contour(Xgrid,Ygrid,data_out(:,:,1));

