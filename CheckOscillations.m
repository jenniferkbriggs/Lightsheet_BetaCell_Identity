%% Analysis of 3D data
%code for network analysis of 3D data
% Jennifer Briggs 06/02/2022
% 
% 
close all 
% clear all
clc
%changing defaults
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)
addpath('~/Documents/GitHub/Functional_and_Structural_Networks')
try
addpath('/Users/brigjenn/Documents/GitHub/UniversalCode')
end


%%Change these!!
pka = 1 %1 if there is there PKA, 0 if control
fileloc = ('/Volumes/Briggs_10TB/Merrin/EJ155 New data set/') %directory with folders for each experiment
files = dir(fileloc);
files_dir = {};
for i = 1:length(files)
    if isfolder([fileloc files(i).name]) & ~contains(files(i).name, '.')
        files_dir{end+1} = {files(i).name}
    end
end

%% Start analysis
for i = 2:length(files_dir)
fullfile = [fileloc files_dir{i}{1} '/']
name = files_dir{i}{1}


savename = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange'

%%Here we load the files
    calcium = readmatrix([fullfile name(1:end-10) 'Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv

    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
   
    photobleaching= 0%input('Is there photobleaching? 1 for yes, 0 for no')
    calstore = calcium;
    timestore = time;
    fignum = 1
   
       load([fullfile 'waveindx_samenumber.mat']);
        %[start_indx_repol, end_indx_repol] = identify_oscillations_end(end_indx, calcium, time)

   
        figure
        plot(mean(calcium')) %plot calcium
      
        xline(cuttime, 'linewidth',3, 'color', 'b')

        for i = 1:length(start_indx_repol)
        hold on, xline(start_indx_repol(i), 'label',['start: ' num2str(i)]), xline(end_indx_repol(i), 'color','r', 'label', ['end:' num2str(i)])
        end
        title(name)
        keyboard


        figure
        plot(mean(calcium')) %plot calcium
      
        xline(cuttime, 'linewidth',3, 'color', 'b')

        for i = 1:length(start_indx)
            hold on, xline((start_indx(i)), 'label',['start: ' num2str(i)]), xline((end_indx(i)), 'color','r', 'label', ['end:' num2str(i)])
        end
        keyboard


        save([fullfile, 'waveindx_samenumber.mat'], 'start_indx','end_indx', 'start_indx_repol', 'end_indx_repol','cuttime','numtrial')
        close all

  end
