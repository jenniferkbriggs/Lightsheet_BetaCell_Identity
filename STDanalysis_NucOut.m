function[Correlation, pixelsx, pixelsy] = STDanalysis_NucOut(images,NucLoc,opts)
addpath('~/Documents/GitHub/UniversalCode/')
%Jennifer Briggs 03.2022
%This function takes a nucleus location and itteratively pulls more pixels
%around that location until the pixels on the outermost edge appear to not
%be apart of the cell. Whether or not the pixels are a part of the cell is
%based on the assumption that pixels from the same cell will have extremely
%similar behavior (minus noise). Therefore, once the behavior of pixels on
%the outermost radius of the cell deviates too much (based on the threshold
%(thr), the cell radius and corresponding timecourse is exported.

%images is the image file in matrix form (x,y,t). Currently only accepting
%3D matrices (e.g. no z stacks) so z stacks should be fed in separately.

%NucLoc is the [X,Y] position of the center of the nucleus (note, if this
%changes overtime, take the median for the most accurate result).

%thr is the threshold for how different the pixels can behave. The
%threshold is based on the standard deviation

%opts is a structure of options:
%opts.figs = 1 if want a gif to see the pixels be removed
    
    radius = 2; %inital radius 
    maxradius = 40;
    g = 0
    %pull in baseline oscillation
    [image_base,allpix] = getcal_radius(radius, images, NucLoc,1);
    try
        if opts.gif
        disp([opts.title '.gif'])
        giffig = figure, ti = tiledlayout(2,2), ti.TileSpacing = 'compact'
        caim = insertMarker(opts.ca, allpix,'color','r', 'Size', 1);
         nexttile(1), imshow(caim)
         xlim([NucLoc(1)-maxradius, NucLoc(1)+maxradius])
         ylim([NucLoc(2)-maxradius, NucLoc(2)+maxradius])
        %nexttile(1), plot(allpix(:,1),allpix(:,2),'o')
        nexttile(2), plot([1:size(image_base,1)], (image_base))
        sgtitle([opts.title])

        gif(['/Volumes/Briggs_10TB/Merrin/Confocal/AllGif' opts.title '.gif'])
        g = 1;
        end
        
    end
    
    Y = mean(image_base,2);
    p = 0; j = 1;
    check = 0;
    Correlation = zeros(maxradius^2, maxradius);
    pixelsx =  zeros(maxradius^2, maxradius);
    pixelsy =  zeros(maxradius^2, maxradius);

    for j = 1:maxradius %keep expanding until error is less than the threshold
        radius = radius+1;
        [image_new,allpix_new] = getcal_radius(radius, images, NucLoc, 0);
        
        X = mean(image_new,2);
        [r, p] = corr(image_new,Y);
        Correlation(1:length(r),j) = r;
        pixelsx(1:length(allpix_new), j) = allpix_new(:,1);
        pixelsy(1:length(allpix_new), j) = allpix_new(:,2);
%         nexttile(1), caim_new = insertMarker(opts.ca, allpix_new, 'color', 'g', 'size',1);
%         imshow(caim_new)
%         xlim([NucLoc(1)-maxradius, NucLoc(1)+maxradius])
%         ylim([NucLoc(2)-maxradius, NucLoc(2)+maxradius])
%         %plot(allpix_new(:,1),allpix_new(:,2),'o'),hold on,  plot(allpix(:,1), allpix(:,2),'o'), legend('Comparison','New'), title('Comparison Radius'), hold off
%         nexttile(2), plot([1:size(image_new,1)], X), hold on, plot([1:size(image_base,1)], Y), hold off, title('Calcium Fluoresence')
%         xlabel('Time')
%         legend('New', 'Baseline')
%         drawnow
%         %pause(1)
%         nexttile(3), plot(radius, r, 'o'), title('All Correlations'), hold on,
%         xlabel('Radius [pix]')
%         
%         nexttile(4), plot(radius, mean(r), 'o'), title('Correlation Average'), hold on,
%         xlabel('Radius [pix]')
        j = j+1;
        if g == 1
        gif
        end
        
        Y = mean(image_new,2);
    end
    %try
    %    saveas(gcf, ['/Volumes/Briggs_10TB/Merrin/Confocal/' opts.title '.fig'])
    %end
    %close(giffig)
    end

    
    function [casmooth,allpix] = getcal_radius(radius, images, NucLoc, filled)
    %function that finds pixels in the circle around the nucleus with a
    %given radius, extracts the calcium oscillations for these pixels and
    %then smooths the data using a gaussian smoother 
    if filled %filled determines whether we want just the outside layer or all of the pixels
    circ = @(radius, NucLoc) unique([reshape((round([1:radius]'.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round([1:radius]'.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');
    else
    circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');
    end
    allpix = circ(radius, NucLoc);
    
    for i = 1:length(allpix)
        calcium(:,i) = squeeze(images(allpix(i,1), allpix(i,2),:));
    end
     
    %Smooth data along time we assume that the noise that is smooted 
    %during this is random noise that contains no information about the
    %correlation between pixels. It may be good to try without
    %smoothing to see the change
   
    casmooth = smoothdata(calcium,'gaussian'); 
    end
