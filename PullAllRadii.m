%%Estimate radius of all cells in the nuclei
close all
clear all
clc

th_rad = 0.9
th_pix =  0.8183

addpath('~/Documents/GitHub/UniversalCode/');

%chose islet to analyze
filename = ["three","sample","five", "two","one"];
csvname = ["H2BmCherry Ucn3GCaMP-3_Detailed.csv", "H2BmCherry Ucn3GCaMP sample_Detailed.csv", "H2BmCherry Ucn3GCaMP-5_Detailed.csv", "H2BmCherry Ucn3GCaMP-2_Detailed.csv", "H2BmCherry Ucn3GCaMP-1_Detailed.csv"]
datapath = ['/Volumes/Briggs_10TB/Merrin/Confocal/'] 

maxradius = 45;
circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');

%initialize which time index to use nuclear locations
%timetouse = 243%three: 217 - somewhat wiggly. %sample: 194; %five: 243%two: 363; %one: 193
timetouse = [217, 194, 243, 363, 193];

for kt = 1:length(filename)
    isletfig = figure
    gif(char(strrep(strjoin(['/Volumes/Briggs_10TB/Merrin/Confocal/AllGif/' filename(kt) '.gif']),' ', '')))

  ca_files = dir(strrep(strjoin([datapath filename(kt) '/' '*C1*.tif']),' ',''));
nuc_files = dir(strrep(strjoin([datapath filename(kt) '/' '*C2*.tif']),' ',''));

    %% Import Images %%
    %import calcium
    for i = 1:length(ca_files)
        fulldatapath = strrep(strjoin([datapath filename(kt) '/' ca_files(i).name]),' ','');
        image = imread(fulldatapath);
        Islet_vid(:,:,i) = image; %calcium files are loaded as a X pixel x Y pixel x Time
        %imshow(imcomplement(image));, drawnow
        fullnucpath = strrep(strjoin([datapath filename(kt) '/' nuc_files(i).name]),' ', '');
        nuimage = imread(fullnucpath);
        Nuc_vid(:,:,i) = nuimage;
    end

    Islet_vid = rot90(Islet_vid, 2);%the image is reflected rotate it back
    Nuc_vid = rot90(Nuc_vid, 2);
    ca_im = rescale(mean(Islet_vid(:,:,timetouse(kt)-10:timetouse(kt)+10),3),0,1);%-double(min(Islet_vid,[],3)))./(double(max(Islet_vid,[],3))-double(min(Islet_vid,[],3)));%mean2([squeeze(mean(Islet_vid,2)),squeeze(mean(Islet_vid,1))]));


    imshow(ca_im)
   %% Load nucleus location 
        nucloc = readtable(strrep(strjoin([datapath csvname(kt)]),'/ ','/')); %import nucleus location csv
        cells_at_time = find(table2array(nucloc(:,7))==timetouse(kt));
        X = (table2array(nucloc(cells_at_time,1)));
        Y = table2array(nucloc(cells_at_time,2));

        nuimage = Nuc_vid(:,:,timetouse(kt));
        loc = [X,Y]; 

        imnew = insertMarker(nuimage, loc);
        figure, imshow(imnew)
        
        NucLoc = fliplr(loc);

    
    
    
    %Because of that rotation, my nucleus location is rotated: flip back
    imAv = mean(Nuc_vid, 3);
    imAv = imAv/max(imAv(:));
    for i = 2:length(NucLoc)
        opts.ca = ca_im;
        opts.gif = 0;
        opts.normalized = 1;
        opts.title = [filename(kt) 'Islet' num2str(i)]
        [Correlation, Pixelsx, Pixelsy] = STDanalysis_NucOut(Islet_vid,NucLoc(i,:),opts);
        
         
       radii = 1;
       score = 1;
       %set threshold
       while score > th_rad
           %if percent of pixels in radius is greater than threshold, we can increase radius
            [pix_belowth(radii)] = length(find(nonzeros(Correlation(:,radii)) > th_pix));
            [tot_pix] = length(nonzeros(Correlation(:,radii)));
            score = pix_belowth(radii)/tot_pix-radii/100;
            radii = radii+1;
       end
            
       nuclocforimages = loc;
        est_radius(kt,i) = radii; 
        figure(isletfig)
        [estlocations] = circ(est_radius(kt,i), loc(i,:));

        caim_new = insertMarker(opts.ca, estlocations, '*', 'color', 'g', 'size',1);
        imshow(caim_new)
        xlim([loc(i,1)-maxradius, loc(i,1)+maxradius])
        ylim([loc(i,2)-maxradius, loc(i,2)+maxradius])
        title(['Radius ' num2str(radii)])
        %est_radius_nores(kt,i) = input('Is there any calcium here? 1 if yes, 0 if no')

        hold off
        gif
        pause(0.5)
end
    close(isletfig)
    
end

save([datapath 'AllCellResults.mat'], 'est_radius')
