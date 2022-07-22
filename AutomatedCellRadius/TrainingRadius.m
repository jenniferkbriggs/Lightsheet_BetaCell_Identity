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
    end

    Islet_vid = rot90(Islet_vid, 2);%the image is reflected rotate it back


    try 
       load(strrep(strjoin(['/Volumes/Briggs_10TB/Merrin/Confocal/' filename(kt) 'SecondAttempt_Analysis.mat']), ' ', ''))
    end
    
    
    %Because of that rotation, my nucleus location is rotated: flip back

    for i = 1:perislet
        opts.ca = Islet_vid(:,:,timetouse(kt));
        opts.gif = 0;
        opts.normalized = 1;
        opts.title = [filename(kt) 'Islet' num2str(i)]
        [FinalCordata(i).Correlation, FinalCordata(i).Pixelsx, FinalCordata(i).Pixelsy] = STDanalysis_NucOut(Islet_vid,NucLoc(i,:),opts);
    end
    save(strrep(strjoin(['/Volumes/Briggs_10TB/Merrin/Confocal/' filename(kt) 'SecondAttempt_radius.mat']), ' ', ''), 'CellMask','FinalCordata','NucLoc')
end
