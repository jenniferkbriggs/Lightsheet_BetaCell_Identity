%% Analysis of 3D data
% In this code, we use network and wave analysis to compare results across
% 3D and 2D. 


% close all 
% clear all
% clc
%changing defaults
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)
addpath('~/Documents/GitHub/Functional_and_Structural_Networks')
try
addpath('/Users/brigjenn/Documents/GitHub/UniversalCode')
end

%% Change these!!
pka = 0 %1 if there is there PKA, 0 if control
fileloc = ('/Volumes/Briggs_10TB/Merrin/2023Data/EJ129 10G/')
files = dir(fileloc);
files_dir = {};
for i = 1:length(files)
    if isfolder([fileloc files(i).name]) & ~contains(files(i).name, '.')
        files_dir{end+1} = {files(i).name}
    end
end

for i = 4%1:length(files_dir)
fullfile = [fileloc files_dir{i}{1} '/']
name = files_dir{i}{1}
savename = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/'

%% Here we load the files
    calcium = readmatrix([fullfile 'Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
    Locations = readmatrix([fullfile 'Pos.csv']);
    Locations = Locations(:,1:3);
    figure(2), plot(mean(calcium'))
    title('Is there photobleaching? 1 for yes, 0 for no')
    photobleaching= 1%input('Is there photobleaching? 1 for yes, 0 for no')
    close(figure(2))
    if pka
    figure, plot(mean(calcium'))
    title('Select beginning of second phase, beginning of pka administration, and end of time course')
    [cuttime, ~] = ginput(3)
    close(figure(1))
    end
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
[start_indx, end_indx] = identify_oscillations(calcium, time,0)
xline(start_indx, 'label','This is where we start')
xline(end_indx, 'label','This is where we end')
close figure 7

%% Find the 2D slice (along Z-plane)
middle_Z = mean(Locations(:,3));
cells_2D = find(abs(middle_Z - Locations(:,3))< 3);
figure, plot3(Locations(:,1),Locations(:,2),Locations(:,3),'o'), hold on,  plot3(Locations(cells_2D,1),Locations(cells_2D,2),Locations(cells_2D,3),'ro', 'MarkerFaceColor','r')
calcium_2D = calcium(:,cells_2D);

%% Run network analysis
network2 = Network_3D(calcium_2D, cuttime, numtrial, photobleaching, [savename '2D'], fileloc, name, Locations(cells_2D,:), timestore,start_indx, end_indx, 1)
network3 = Network_3D(calcium, cuttime, numtrial, photobleaching, [savename '3D'], fileloc, name, Locations, timestore,start_indx, end_indx, 1)


network_all.(name(end-2:end)).twoD = network2;
network_all.(name(end-2:end)).threeD = network3;
%% Run waveintiator analysis
wave2 = Wave_3D(calcium_2D, cuttime, numtrial, photobleaching, [savename '2D'], fileloc, name, Locations(cells_2D,:), timestore,start_indx, end_indx, time, 1)
wave3 = Wave_3D(calcium, cuttime, numtrial, photobleaching, [savename '3D'], fileloc, name, Locations, timestore,start_indx, end_indx, time, 1)

wave_all.(name(end-2:end)).twoD = network2;
wave_all.(name(end-2:end)).threeD = network3;
end
%%
figure, bar(cells_2D, mean(wave2.phase));
hold on, plot(mean(wave3.phase),'r*')