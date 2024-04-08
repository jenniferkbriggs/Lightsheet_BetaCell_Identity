%% Analysis of 3D data
% In this code, we use network and wave analysis to compare results across
% 3D and 2D. 
% 
% close all 
% clear all
% clc
% %changing defaults
% load('/Users/brigjenn/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Anschutz/Islet/3DLightSheet/CORRECTEDWAVE_Analysis_After_ChamberChange20240326tenperc_allAnalysis.mat')
% wave3 = wave_all;
% network3 = network_all;
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)
addpath(genpath('~/Documents/GitHub/Functional_and_Structural_Networks'))
try
addpath('/Users/brigjenn/Documents/GitHub/UniversalCode')
end



savetime =  datestr(datetime('today'),'yyyymmdd');
%%Change these!!
pka = 1 %1 if there is there PKA, 0 if control
fileloc = ('/Volumes/Briggs_10TB/Merrin/EJ155 New data set/') %directory with folders for each experiment
files = dir(fileloc);
files_dir = {};
for i = 5:length(files)
    if isfolder([fileloc files(i).name]) & ~contains(files(i).name, '.')
        files_dir{end+1} = {files(i).name}
    end
end

%% Start analysis
for i = 1:length(files_dir)
fullfile = [fileloc files_dir{i}{1} '/']
name = files_dir{i}{1}
savename = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange/2Dvs3D'

%% Here we load the files

%%Here we load the files
try
    Locations = readmatrix([fullfile name(1:end-10) 'detailed.csv']);
catch
    Locations = readmatrix([fullfile name(1:end-10) 'Detailed.csv']);
end
    Locations = Locations(:,1:3);
    calcium = readmatrix([fullfile name(1:end-10) 'Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv

    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
   
    photobleaching= 0%input('Is there photobleaching? 1 for yes, 0 for no')
    calstore = calcium;
    timestore = time;
    fignum = 1

    load([fullfile 'waveindx_samenumber.mat']);

%% Find the 2D slice (along Z-plane)
minny = min(Locations(:,3));
maxy = max(Locations(:,3));
middle_Z = range(Locations(:,3))./2 + minny;


cells_middle = find(abs(middle_Z - Locations(:,3))< 3);
calcium_middle = calcium(:,cells_middle);

thirdsZ = 3*range(Locations(:,3))./4 + minny;
cells_thirds = find(abs(thirdsZ - Locations(:,3))< 3);
figure, p1 = plot3(Locations(:,1),Locations(:,2),Locations(:,3),'o'), hold on,  p2 = plot3(Locations(cells_thirds,1),Locations(cells_thirds,2),Locations(cells_thirds,3),'bo', 'MarkerFaceColor','b')
p3 = plot3(Locations(cells_middle,1),Locations(cells_middle,2),Locations(cells_middle,3),'ro', 'MarkerFaceColor','r')
legend([p2, p3], '1/4 Z-stack','1/2 Z-stack')
calcium_third = calcium(:,cells_thirds);

if i == 29
    middle_Z = -50;


    cells_middle = find(abs(middle_Z - Locations(:,3))< 3);
    calcium_middle = calcium(:,cells_middle);
    
    thirdsZ = -30;
    cells_thirds = find(abs(thirdsZ - Locations(:,3))< 3);
    figure, p1 = plot3(Locations(:,1),Locations(:,2),Locations(:,3),'o'), hold on,  p2 = plot3(Locations(cells_thirds,1),Locations(cells_thirds,2),Locations(cells_thirds,3),'bo', 'MarkerFaceColor','b')
    p3 = plot3(Locations(cells_middle,1),Locations(cells_middle,2),Locations(cells_middle,3),'ro', 'MarkerFaceColor','r')
    legend([p2, p3], '1/4 Z-stack','1/2 Z-stack')
    calcium_third = calcium(:,cells_thirds);
end

name = name(find(~isspace(name)));
%remove dashes:
toremove = strfind(name, '-');
name(toremove) = [];

saveas(gcf, [savename '/' name '2Dvs3D.fig'])
saveas(gcf, [savename '/' name '2Dvs3D.png'])
end
%% Run network analysis
% network_threequarter = Network_3D(calcium_third, cuttime, numtrial, photobleaching, [savename '/threequarter'], fileloc, name, Locations(cells_thirds,:), timestore,start_indx, end_indx);
% network_half = Network_3D(calcium_middle, cuttime, numtrial, photobleaching, [savename '/half'], fileloc, name, Locations(cells_middle,:), timestore,start_indx, end_indx);
% 
% 
% network_all.(name).half2D = network_half;
% network_all.(name).half2D.cellindex = cells_middle;
% network_all.(name).threequarter2D = network_threequarter;
% network_all.(name).threequarter2D.cellindex = cells_thirds;
% network_all.(name).threeD = network3.(name);
% 
% close all

%% Run waveintiator analysis

wave_threequarter = Wave_3D(calcium_third, cuttime, numtrial, photobleaching, [savename '/threequarter'], fileloc, name, Locations(cells_thirds,:), timestore,start_indx, end_indx, time,0, 1,1);
wave_half = Wave_3D(calcium_middle, cuttime, numtrial, photobleaching, [savename '/half'], fileloc, name,  Locations(cells_middle,:), timestore,start_indx, end_indx, time,0, 1,1);
wave_3 = Wave_3D(calcium, cuttime, numtrial, photobleaching, [savename '/ThreeD'], fileloc, name,  Locations, timestore,start_indx, end_indx, time,0, 1,1);



wave_all.(name).half2D = wave_half;
wave_all.(name).threequarter2D = wave_threequarter;
wave_all.(name).threeD = wave_3;

disp(i)
end

save([savename '/3Dvs2D_betternormalizedwave.mat'], 'wave_all','network_all')
%%
figure, bar(cells_2D, mean(wave2.phase));
hold on, plot(mean(wave3.phase),'r*')