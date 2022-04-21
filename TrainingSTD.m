%%This script is how we define the training data and develop the ROC curve
close all
clear all
clc

addpath('~/Documents/GitHub/UniversalCode/');

%chose islet to analyze
filename = ["three","sample","five", "two","one"];
csvname = ["H2BmCherry Ucn3GCaMP-3_Detailed.csv", "H2BmCherry Ucn3GCaMP sample_Detailed.csv", "H2BmCherry Ucn3GCaMP-5_Detailed.csv", "H2BmCherry Ucn3GCaMP-2_Detailed.csv", "H2BmCherry Ucn3GCaMP-1_Detailed.csv"]
datapath = ['/Volumes/Briggs_10TB/Merrin/Confocal/'] 


%initialize which time index to use nuclear locations
%timetouse = 243%three: 217 - somewhat wiggly. %sample: 194; %five: 243%two: 363; %one: 193
timetouse = [217, 194, 243, 363, 193];
%Number of cells per islet for the training case
perislet = 5
rng(1) %set seed

for kt = 1:length(filename)
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


    try 
       load(strrep(strjoin(['/Volumes/Briggs_10TB/Merrin/Confocal/' filename(kt) 'RetrainAnalysis.mat']), ' ', ''))
    catch
        %% Load nucleus location 
        nucloc = readtable(strrep(strjoin([datapath csvname(kt)]),'/ ','/')); %import nucleus location csv
        cells_at_time = find(table2array(nucloc(:,7))==timetouse(kt));
        X = (table2array(nucloc(cells_at_time,1)));
        Y = table2array(nucloc(cells_at_time,2));

        nuimage = Nuc_vid(:,:,timetouse(kt));
        loc = [X,Y]; 

        imnew = insertMarker(nuimage, loc);
        figure, imshow(imnew)

        trainingcells = randi([1,length(X)], 1,perislet); %select 'perislet' number of training cells at random with generated seed

        %% Grab nucleus location and cell outline
        %create a color plot of nucleus and calcium
        caim_nuc = insertMarker(ca_im, loc(:,:));
        fig = figure('Name','Draw ROIs Around Cells of Interest');
        fig.Position = [1293 366 1076 918];
        fpos = get(fig,'Position')
        axOffset = (fpos(3:4)-[size(caim_nuc,2) size(caim_nuc,1)])/2;
        ha = axes('Parent',fig,'Units','pixels',...
                    'Position',[axOffset size(caim_nuc,2) size(caim_nuc,1)]);
        imshow(caim_nuc, 'Parent',ha);
        CellMask = double(zeros(size(ca_im)));
        for i = 1:perislet
            imshow(caim_nuc)
            title('Select Nucleus')
            training_loc = ginput(1);
            NucLoc_Ithink = find(sum(abs(loc - training_loc),2)<5);
            if length(NucLoc_Ithink)~= 1
                keyboard
            end
            NucLoc(i,:) = loc(NucLoc_Ithink,:);
            caim_onenuc = insertMarker(ca_im, NucLoc(i,:));
            imshow(caim_onenuc)
            title('Draw Around Cell')
            keyboard
            ROIMask = imfreehand(); %User draws region around cell
            ROIMask = createMask(ROIMask); %Mask is created from drawn region
            CellMask = CellMask + ROIMask.*i; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
            caim_nuc = insertMarker(caim_nuc, NucLoc(i,:), 'color','r');
        end
        imshow(caim_nuc)
        NucLoc = fliplr(NucLoc);

    end
    
    
    %Because of that rotation, my nucleus location is rotated: flip back

    for i = 1:perislet
        opts.ca = ca_im;
        opts.gif = 0;
        opts.normalized = 1;
        opts.title = [filename(kt) 'Islet' num2str(i)]
        [FinalCordata(i).Correlation, FinalCordata(i).Pixelsx, FinalCordata(i).Pixelsy] = STDanalysis_NucOut_Training(Islet_vid,NucLoc(i,:),opts);
    end
    save(strrep(strjoin(['/Volumes/Briggs_10TB/Merrin/Confocal/' filename(kt) 'NormalizeAnalysis.mat']), ' ', ''), 'CellMask','FinalCordata','NucLoc')


end
