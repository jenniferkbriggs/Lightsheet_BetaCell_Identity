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


savename = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange'

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
 try 
       load([fullfile 'waveindx_samenumber.mat']);

  catch

        if pka
            figure, plot(mean(calcium')) 
            title('Select beginning of second phase, beginning of pka administration, and end of time course')
            admin = find(abs(time - 1800)<0.2); %when did the administration occur
            xline(admin(1))
            [cuttime, ~] = ginput(3) % if there is a bump or blip at administration, select BEFORE the blip - it will be excluded in analysis
            close(figure(1))
        end

        if pka
            numtrial = 2
        else
            numtrial = 1
            cuttime(1) = 1
            cuttime(2) = length(calcium)
        end
        figure, plot(time, mean(calcium'))
        disp("Resize the figure and then click continue after you are happy with it")
        Opts.direction_hold = 0
        Opts.figs = 0;
        xline(time(round(cuttime(2))), 'linewidth',3)

        %identify depolarization
        [start_indx, end_indx] = identify_oscillations(calcium, time, 0)
if ~exist('end_indx_repol')
        %identify repolarization
        [start_indx_repol, end_indx_repol] = identify_oscillations_end(end_indx, calcium, time)
end
        numosc = length(end_indx_repol)./2;
        if end_indx_repol(numosc) > cuttime(2) %make the end phase analysis isn't including the activator
            end_indx_repol(numosc) ;
        end
        keyboard
        save([fullfile, 'waveindx_samenumber.mat'], 'start_indx','end_indx', 'start_indx_repol', 'end_indx_repol','cuttime','numtrial')
  end


  %quantify peak activation: 
  

%start_indx - start_indx (2)

for j = 1:length(start_indx)/2
    if j == length(start_indx)/2
        wave_end = cuttime(2);
    else
        wave_end = start_indx(j+1);
    end  

cal_max_pre(i, j) = mean(max(calcium(start_indx(j):wave_end,:)));
end


for j = length(start_indx)/2+1:length(start_indx)
    if j == length(start_indx)
        wave_end = cuttime(3);
    else
        wave_end = start_indx(j+1);
    end  

cal_max_post(i, j-length(start_indx)/2) = mean(max(calcium(start_indx(j):wave_end,:)));
end
cal_norm(i) = range(calcium(:));
end

Pre_Control = cal_max_pre(Control,:)

%fit a photobleaching line: 
for k = 1:length(Control)
    x = 1:length(nonzeros(cal_max_pre(Control(k),:)));
    P1 = polyfit(x, nonzeros(cal_max_pre(Control(k),:))./cal_norm(Control(k)), 1);
    P2 = polyfit(x, nonzeros(cal_max_post(Control(k),:))./cal_norm(Control(k)), 1);
    %P(k,:) = polyfit([cal_max_pre(Control(k),:), cal_max_post(Control(k),:)],1);
    Control_slope(k,:) = [P1(1), P2(1)]; %slope
end

%fit a photobleaching line: 
for k = 1:length(GKa)
    x = 1:length(nonzeros(cal_max_pre(GKa(k),:)));
    P1 = polyfit(x, nonzeros(cal_max_pre(GKa(k),:))./cal_norm(GKa(k)), 1);
    P2 = polyfit(x, nonzeros(cal_max_post(GKa(k),:))./cal_norm(GKa(k)), 1);
    %P(k,:) = polyfit([cal_max_pre(Control(k),:), cal_max_post(Control(k),:)],1);
    GKa_slope(k,:) = [P1(1), P2(1)]; %slope
end

%fit a photobleaching line: 
for k = 1:length(PKa)
    x = 1:length(nonzeros(cal_max_pre(PKa(k),:)));
    P1 = polyfit(x, nonzeros(cal_max_pre(PKa(k),:))./cal_norm(PKa(k)), 1);
    P2 = polyfit(x, nonzeros(cal_max_post(PKa(k),:))./cal_norm(PKa(k)), 1);
    %P(k,:) = polyfit([cal_max_pre(Control(k),:), cal_max_post(Control(k),:)],1);
    PKa_slope(k,:) = [P1(1), P2(1)]; %slope
end
