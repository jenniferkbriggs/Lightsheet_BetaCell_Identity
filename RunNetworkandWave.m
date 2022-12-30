%% Analysis of 3D data
%code for network analysis of 3D data
% Jennifer Briggs 06/02/2022


close all 
clear all
clc
%changing defaults
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)
addpath('~/Documents/GitHub/Functional_and_Structural_Networks')
try
addpath('/Users/brigjenn/Documents/GitHub/UniversalCode')
end

%% Change these!!
pka = 0
fileloc = ('/Volumes/Briggs_10TB/Merrin/Ca_Courses/Ca_Courses_new/EJ129 10G/EJ129 10G F03_Statistics/')
filename = 'F03 10g'
fullfile = [fileloc filename]
savename = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/'

%% Here we load the files
    calcium = readmatrix([fullfile '_Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
    calcium(1:3,:) = []; %the CSV you have has the first 3 rows as NAN so we remove them
    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
    Locations = readmatrix([fullfile '_Pos.csv']);
    Locations = Locations(:,1:3);
    figure, plot(mean(calcium'))
    title('Is there photobleaching? 1 for yes, 0 for no')
    photobleaching= input('Is there photobleaching? 1 for yes, 0 for no')
    if pka
    figure, plot(mean(calcium'))
    title('Select beginning of second phase, beginning of pka administration, and end of time course')
    [cuttime, ~] = ginput(3)
    close(figure(1))
    end
    calstore = calcium;
    timestore = time;
    fignum = 1

    if pka
        numtrial = 2
    else
        numtrial = 1
        cuttime(1) = 1
        cuttime(2) = length(calcium)
    end

    
figure(7), plot(time,calcium) %plot calcium
disp("Resize the figure and then click continue after you are happy with it")
Opts.direction_hold = 1
Opts.figs = 0;
[start_indx, end_indx] = identify_oscillations(calcium, time, 0)
xline(start_indx, 'label','This is where we start')
xline(end_indx, 'label','This is where we end')

close figure 7
%% Run network analysis
network = Network_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, filename, Locations, calstore, timestore,start_indx, end_indx)
%% Run waveintiator analysis
wave = Wave_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, filename, Locations, calstore, timestore,start_indx, end_indx, time)