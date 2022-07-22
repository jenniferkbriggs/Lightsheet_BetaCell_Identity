function[Correlation, pixelsx, pixelsy] = STDanalysis_NucOut_Training(images,NucLoc,opts)
addpath('~/Documents/GitHub/UniversalCode/')
%Jennifer Briggs 03.2022
%This function calculates the correaltion of all pixels with in the
%'maxradius' of the nucleus. 

%images is the image file in matrix form (x,y,t). Currently only accepting
%3D matrices (e.g. no z stacks) so z stacks should be fed in separately.

%NucLoc is the [X,Y] position of the center of the nucleus (note, if this
%changes overtime, take the median for the most accurate result).

%thr is the threshold for how different the pixels can behave. The
%threshold is based on the standard deviation

%opts is a structure of options:
%opts.figs = 1 if want a gif to see the pixels be removed
    images = double(images)+0.01; %note that it rotates 90 degrees again. not sure why
    radius = 2; %inital radius 
    maxradius = 40;
    g = 0;
    %pull in baseline oscillation
    [image_base,allpix] = getcal_radius(radius, images, NucLoc,1);
    try
        if opts.gif
        disp([opts.title '.gif'])
        figure, ti = tiledlayout(2,2), ti.TileSpacing = 'compact'
        caim = insertMarker(opts.ca, allpix,'*', 'color','r', 'Size', 1);
         nexttile(1), imshow(caim)
         xlim([NucLoc(1)-maxradius-10, NucLoc(1)+maxradius+10])
         ylim([NucLoc(2)-maxradius-10, NucLoc(2)+maxradius+10])
        %nexttile(1), plot(allpix(:,1),allpix(:,2),'o')
        nexttile(2), plot([1:size(image_base,1)], (image_base))
        g = 1;
        end
        
    end
    
    Y = mean(image_base,2);
    p = 0; j = 1;
    check = 0;



        [image_new,allpix_new] = getcal_radius(maxradius, images, NucLoc, 1);
        
        [r, p] = corr(image_new,Y);
        if opts.normalized == 1
        r = r/mean(r);
        end
        Correlation = r;
        pixelsx = allpix_new(:,1);
        pixelsy = allpix_new(:,2);
       
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
