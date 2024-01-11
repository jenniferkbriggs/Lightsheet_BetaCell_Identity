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

ct = 1


%% Change these!!
pka = 1 %1 if there is there PKA, 0 if control
fileloc = ('/Volumes/Briggs_10TB/Merrin/2023Data/EJ106-113 10GPKa/')
files = dir(fileloc);
files_dir = {};
for i = 1:length(files)
    if isfolder([fileloc files(i).name]) & ~contains(files(i).name, '.')
        files_dir{end+1} = {files(i).name}
    end
end

for i = 1:length(files_dir)
fullfile = [fileloc files_dir{i}{1} '/']
newfiles = dir(fullfile);
newfiles_a = {};
for ll = 1:length(newfiles)
    if isfolder([fullfile newfiles(ll).name]) & ~contains(newfiles(ll).name, '.')
        newfiles_a{end+1} = {newfiles(ll).name}
    end
end

for j = 1:length(newfiles_a)
fullfile = [fileloc files_dir{i}{1} '/' newfiles_a{j}{1} '/']
savename = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/'
name = newfiles_a{j}{1};
%%Here we load the files
    %calcium = readmatrix([fullfile 'Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
    calcium = readmatrix([fullfile 'Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
    figure, plot(time, calcium)
    title(num2str(j))

    try
    Locations = readmatrix([fullfile 'Pos.csv']);
    catch
        Locations = readmatrix([fullfile ' pos_Detailed.csv']);
    end
    Locations = Locations(:,1:3);
%     figure, plot(mean(calcium'))
%     title('Is there photobleaching? 1 for yes, 0 for no')
    photobleaching= 0;%input('Is there photobleaching? 1 for yes, 0 for no')

    calstore = calcium;
    timestore = time;
    fignum = 1

   try 
       load([fullfile, 'waveindx_samenumber.mat']); 
   catch
           
        if pka
            figure, plot(mean(calcium'))
            title('Select beginning of second phase, beginning of pka administration, and end of time course')
            xline(1600*2)
            [cuttime, ~] = ginput(3)
            close(figure(1))
        end

        if pka
            numtrial = 2
        else
            numtrial = 1
            cuttime(1) = 1
            cuttime(2) = length(calcium)
        end
    
    
        figure(7), plot(time, mean(calcium')) %plot calcium
        disp("Resize the figure and then click continue after you are happy with it")
        Opts.direction_hold = 1
        Opts.figs = 0;
        xline(time(round(cuttime(2))), 'linewidth',3)
        [start_indx, end_indx] = identify_oscillations(calcium, time, 0)
        xline(start_indx, 'label','This is where we start')
        xline(end_indx, 'label','This is where we end')
         keyboard
        save([fullfile, 'waveindx_samenumber.mat'], 'start_indx','end_indx', 'cuttime','numtrial')
   end
   

%close figure 7
%%Run network analysisclose all
for u = 1:length(start_indx)
    if start_indx(u) < 0
        start_indx(u) = 1;
    end
end
%network = Network_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, name, Locations, timestore,start_indx, end_indx)
%%Run waveintiator analysis
wave = Wave_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, name, Locations, timestore,start_indx, end_indx, time, 1, 1)

try
wave_all.(name(end-3:end))= wave;
%network_all.(name(end-3:end))= network;
wave_all.(name(end-3:end)).Locations= normalize(Locations, "range");
catch
wave_all.(name(end-2:end))= wave;
%network_all.(name(end-2:end))= network;
wave_all.(name(end-2:end)).Locations=  normalize(Locations, "range");;

end
ct = ct+1
end
end

save([savename 'PKAAnalyses.mat'], 'network_all','wave_all')


%% 

fnames = fieldnames(wave_all);

% Range of phases:


for i = 1:length(fnames)
    numosc = length(wave_all.(fnames{i}).phaserange)/2
    pre(i) = mean(wave_all.(fnames{i}).phaserange(1:numosc));
    post(i) = mean(wave_all.(fnames{i}).phaserange(numosc+1:end));
end

phaseranget = array2table([pre', post'])
writetable(phaseranget, [savename 'PKARange_of_phases.csv'])



wave_cell = array2table([mean(nonzeros(wave_all.(fnames{1}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{1}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{2}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{2}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{3}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{3}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{4}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{4}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{5}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{5}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{6}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{6}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{7}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{7}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{8}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{8}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{9}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{9}).Degree_KL_wave(2).trial)) 
mean(nonzeros(wave_all.(fnames{10}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{10}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{11}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{11}).Degree_KL_wave(2).trial)) 
mean(nonzeros(wave_all.(fnames{12}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{12}).Degree_KL_wave(2).trial))
], 'VariableNames',{'Pre','PKA'})

writetable(wave_cell, [savename 'PKAWaveCell.csv'])

network_cell = array2table([mean(nonzeros(network_all.(fnames{1}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{1}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{2}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{2}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{3}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{3}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{4}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{4}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{5}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{5}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{6}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{6}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{7}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{7}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{8}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{8}).Degree_KL_net(2).trial)) 
mean(nonzeros(network_all.(fnames{9}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{9}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{10}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{10}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{11}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{11}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{12}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{12}).Degree_KL_net(2).trial))], 'VariableNames',{'Vehicle','PKA'})
writetable(network_cell, [savename 'PKANetCell.csv'])

wave_region = array2table([mean(nonzeros(wave_all.(fnames{1}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{1}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{2}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{2}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{3}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{3}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{4}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{4}).COG_KL_wave(2).trial)) 
mean(nonzeros(wave_all.(fnames{5}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{5}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{6}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{6}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{7}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{7}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{8}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{8}).COG_KL_wave(2).trial)) 
mean(nonzeros(wave_all.(fnames{9}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{9}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{10}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{10}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{11}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{11}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{12}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{12}).COG_KL_wave(2).trial))
], 'VariableNames',{'Pre','PKA'})
writetable(wave_region, [savename 'PKAWaveRegion.csv'])

network_region = array2table([mean(nonzeros(network_all.(fnames{1}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{1}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{2}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{2}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{3}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{3}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{4}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{4}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{5}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{5}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{6}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{6}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{7}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{7}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{8}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{8}).COG_KL_net(2).trial)) 
mean(nonzeros(network_all.(fnames{9}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{9}).COG_KL_net(2).trial)) 
mean(nonzeros(network_all.(fnames{10}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{10}).COG_KL_net(2).trial)) 
mean(nonzeros(network_all.(fnames{11}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{11}).COG_KL_net(2).trial)) 
mean(nonzeros(network_all.(fnames{12}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{12}).COG_KL_net(2).trial))], 'VariableNames',{'Vehicle','PKA'})
writetable(network_region, [savename 'PKANetRegion.csv'])



maintainednet = array2table([network_all.(fnames{1}).mainted_perc;...
    network_all.(fnames{2}).mainted_perc; network_all.(fnames{3}).mainted_perc;...
    network_all.(fnames{4}).mainted_perc; network_all.(fnames{5}).mainted_perc;...
    network_all.(fnames{6}).mainted_perc; network_all.(fnames{7}).mainted_perc;...
    network_all.(fnames{8}).mainted_perc;network_all.(fnames{9}).mainted_perc;...
    network_all.(fnames{10}).mainted_perc;network_all.(fnames{11}).mainted_perc;network_all.(fnames{12}).mainted_perc])
writetable(maintainednet, [savename 'PKANetmaintained.csv'])





pkavar = array2table([wave_all.(fnames{1}).intravar, wave_all.(fnames{1}).intravar_pka, wave_all.(fnames{1}).intervar; 
    wave_all.(fnames{2}).intravar, wave_all.(fnames{2}).intravar_pka, wave_all.(fnames{2}).intervar; 
    wave_all.(fnames{3}).intravar, wave_all.(fnames{3}).intravar_pka, wave_all.(fnames{3}).intervar; 
    wave_all.(fnames{4}).intravar, wave_all.(fnames{4}).intravar_pka, wave_all.(fnames{4}).intervar; 
    wave_all.(fnames{5}).intravar, wave_all.(fnames{5}).intravar_pka, wave_all.(fnames{5}).intervar;
    wave_all.(fnames{6}).intravar, wave_all.(fnames{6}).intravar_pka, wave_all.(fnames{6}).intervar;
    wave_all.(fnames{7}).intravar, wave_all.(fnames{7}).intravar_pka, wave_all.(fnames{7}).intervar; 
    wave_all.(fnames{8}).intravar, wave_all.(fnames{8}).intravar_pka, wave_all.(fnames{8}).intervar; 
    wave_all.(fnames{9}).intravar, wave_all.(fnames{9}).intravar_pka, wave_all.(fnames{9}).intervar; 
    wave_all.(fnames{10}).intravar, wave_all.(fnames{10}).intravar_pka, wave_all.(fnames{10}).intervar; 
    wave_all.(fnames{11}).intravar, wave_all.(fnames{11}).intravar_pka, wave_all.(fnames{11}).intervar; 
    wave_all.(fnames{12}).intravar, wave_all.(fnames{12}).intravar_pka, wave_all.(fnames{12}).intervar; 
], 'VariableNames',{'Vehicle Intravariability','PKA Intravariability', 'Intervariability'})
writetable(pkavar, [savename 'PKAvariability.csv'])

Netpkavar = array2table([network_all.(fnames{1}).intravar, network_all.(fnames{1}).intravar_pka, network_all.(fnames{1}).intervar; 
    network_all.(fnames{2}).intravar, network_all.(fnames{2}).intravar_pka, network_all.(fnames{2}).intervar; 
    network_all.(fnames{3}).intravar, network_all.(fnames{3}).intravar_pka, network_all.(fnames{3}).intervar; 
    network_all.(fnames{4}).intravar, network_all.(fnames{4}).intravar_pka, network_all.(fnames{4}).intervar; 
    network_all.(fnames{5}).intravar, network_all.(fnames{5}).intravar_pka, network_all.(fnames{5}).intervar;
    network_all.(fnames{6}).intravar, network_all.(fnames{6}).intravar_pka, network_all.(fnames{6}).intervar;
    network_all.(fnames{7}).intravar, network_all.(fnames{7}).intravar_pka, network_all.(fnames{7}).intervar; 
    network_all.(fnames{8}).intravar, network_all.(fnames{8}).intravar_pka, network_all.(fnames{8}).intervar; 
    network_all.(fnames{9}).intravar, network_all.(fnames{9}).intravar_pka, network_all.(fnames{9}).intervar; 
    network_all.(fnames{10}).intravar, network_all.(fnames{10}).intravar_pka, network_all.(fnames{10}).intervar; 
    network_all.(fnames{11}).intravar, network_all.(fnames{11}).intravar_pka, network_all.(fnames{11}).intervar; 
    network_all.(fnames{12}).intravar, network_all.(fnames{12}).intravar_pka, network_all.(fnames{12}).intervar; 
], 'VariableNames',{'Vehicle Intravariability','PKA Intravariability', 'Intervariability'})
writetable(Netpkavar, [savename 'Network_PKAvariability.csv'])

%Locations:
longestosc = 0
for i = 1:12
    longestosc = max(longestosc, length(network_all.(fnames{i}).meanLocation(:,1)));
end
for i = 1:12
    ControlLocRadius(i,1:length(network_all.(fnames{i}).meanLocation(:,1))) = network_all.(fnames{i}).meanLocation(:,1);
    ControlLocRadius(i,length(network_all.(fnames{i}).meanLocation(:,1))+1:longestosc) = NaN;
    PKALocRadius(i,1:length(network_all.(fnames{i}).meanLocation_pka(:,1))) = network_all.(fnames{i}).meanLocation_pka(:,1);    
    PKALocRadius(i,length(network_all.(fnames{i}).meanLocation_pka(:,1))+1:longestosc) = NaN;

    ControlLocTheta(i,1:length(network_all.(fnames{i}).meanLocation(:,2))) = network_all.(fnames{i}).meanLocation(:,2);
    ControlLocTheta(i,length(network_all.(fnames{i}).meanLocation(:,2))+1:longestosc) = NaN;
    PKALocTheta(i,1:length(network_all.(fnames{i}).meanLocation_pka(:,2))) = network_all.(fnames{i}).meanLocation_pka(:,2);    
    PKALocTheta(i,length(network_all.(fnames{i}).meanLocation_pka(:,2))+1:longestosc) = NaN;

    ControlLocPhi(i,1:length(network_all.(fnames{i}).meanLocation(:,3))) = network_all.(fnames{i}).meanLocation(:,3);
    ControlLocPhi(i,length(network_all.(fnames{i}).meanLocation(:,3))+1:longestosc) = NaN;
    PKALocPhi(i,1:length(network_all.(fnames{i}).meanLocation_pka(:,3))) = network_all.(fnames{i}).meanLocation_pka(:,3);    
    PKALocPhi(i,length(network_all.(fnames{i}).meanLocation_pka(:,3))+1:longestosc) = NaN;
end

writematrix([ControlLocRadius, PKALocRadius], [savename, 'RadiusStability.csv'])
writematrix([ControlLocTheta, PKALocTheta], [savename, 'ThetaStability.csv'])
writematrix([ControlLocPhi, PKALocPhi], [savename, 'PhiStability.csv'])

for i = 1:12
    ControlLocRadiusw(i,1:length(wave_all.(fnames{i}).MeanLocation(:,1))) = wave_all.(fnames{i}).MeanLocation(:,1);
    ControlLocRadiusw(i,length(wave_all.(fnames{i}).MeanLocation(:,1))+1:longestosc) = NaN;
    PKALocRadiusw(i,1:length(wave_all.(fnames{i}).meanLocation_pka(:,1))) = wave_all.(fnames{i}).meanLocation_pka(:,1);    
    PKALocRadiusw(i,length(wave_all.(fnames{i}).meanLocation_pka(:,1))+1:longestosc) = NaN;

    ControlLocThetaw(i,1:length(wave_all.(fnames{i}).MeanLocation(:,2))) = wave_all.(fnames{i}).MeanLocation(:,2);
    ControlLocThetaw(i,length(wave_all.(fnames{i}).MeanLocation(:,2))+1:longestosc) = NaN;
    PKALocThetaw(i,1:length(wave_all.(fnames{i}).meanLocation_pka(:,2))) = wave_all.(fnames{i}).meanLocation_pka(:,2);    
    PKALocThetaw(i,length(wave_all.(fnames{i}).meanLocation_pka(:,2))+1:longestosc) = NaN;

    ControlLocPhiw(i,1:length(wave_all.(fnames{i}).MeanLocation(:,3))) = wave_all.(fnames{i}).MeanLocation(:,3);
    ControlLocPhiw(i,length(wave_all.(fnames{i}).MeanLocation(:,3))+1:longestosc) = NaN;
    PKALocPhiw(i,1:length(wave_all.(fnames{i}).meanLocation_pka(:,3))) = wave_all.(fnames{i}).meanLocation_pka(:,3);    
    PKALocPhiw(i,length(wave_all.(fnames{i}).meanLocation_pka(:,3))+1:longestosc) = NaN;
end

writematrix([ControlLocRadiusw, PKALocRadiusw], [savename, 'waveRadiusStability.csv'])
writematrix([ControlLocThetaw, PKALocThetaw], [savename, 'waveThetaStability.csv'])
writematrix([ControlLocPhiw, PKALocPhiw], [savename, 'wavePhiStability.csv'])

%ST location:
for i = 1:12
    ControlLocRadiusST(i,1:length(network_all.(fnames{i}).STLocation(:,1))) = network_all.(fnames{i}).STLocation(:,1);
    ControlLocRadiusST(i,length(network_all.(fnames{i}).STLocation(:,1))+1:longestosc) = NaN;
    PKALocRadiusST(i,1:length(network_all.(fnames{i}).STLocation_pka(:,1))) = network_all.(fnames{i}).STLocation_pka(:,1);    
    PKALocRadiusST(i,length(network_all.(fnames{i}).STLocation_pka(:,1))+1:longestosc) = NaN;

    ControlLocThetaST(i,1:length(network_all.(fnames{i}).STLocation(:,2))) = network_all.(fnames{i}).STLocation(:,2);
    ControlLocThetaST(i,length(network_all.(fnames{i}).STLocation(:,2))+1:longestosc) = NaN;
    PKALocThetaST(i,1:length(network_all.(fnames{i}).STLocation_pka(:,2))) = network_all.(fnames{i}).STLocation_pka(:,2);    
    PKALocThetaST(i,length(network_all.(fnames{i}).STLocation_pka(:,2))+1:longestosc) = NaN;

    ControlLocPhiST(i,1:length(network_all.(fnames{i}).STLocation(:,3))) = network_all.(fnames{i}).STLocation(:,3);
    ControlLocPhiST(i,length(network_all.(fnames{i}).STLocation(:,3))+1:longestosc) = NaN;
    PKALocPhiST(i,1:length(network_all.(fnames{i}).STLocation_pka(:,3))) = network_all.(fnames{i}).STLocation_pka(:,3);    
    PKALocPhiST(i,length(network_all.(fnames{i}).STLocation_pka(:,3))+1:longestosc) = NaN;
end
writematrix([ControlLocRadiusST, PKALocRadiusST], [savename, 'RadiusStabilityST.csv'])
writematrix([ControlLocThetaST, PKALocThetaST], [savename, 'ThetaStabilityST.csv'])
writematrix([ControlLocPhiST, PKALocPhiST], [savename, 'PhiStabilityST.csv'])

for i = 1:12
    ControlLocRadiuswST(i,1:length(wave_all.(fnames{i}).STLocation(:,1))) = wave_all.(fnames{i}).STLocation(:,1);
    ControlLocRadiuswST(i,length(wave_all.(fnames{i}).STLocation(:,1))+1:longestosc) = NaN;
    PKALocRadiuswST(i,1:length(wave_all.(fnames{i}).STLocation_pka(:,1))) = wave_all.(fnames{i}).STLocation_pka(:,1);    
    PKALocRadiuswST(i,length(wave_all.(fnames{i}).STLocation_pka(:,1))+1:longestosc) = NaN;

    ControlLocThetawST(i,1:length(wave_all.(fnames{i}).STLocation(:,2))) = wave_all.(fnames{i}).STLocation(:,2);
    ControlLocThetawST(i,length(wave_all.(fnames{i}).STLocation(:,2))+1:longestosc) = NaN;
    PKALocThetawST(i,1:length(wave_all.(fnames{i}).STLocation_pka(:,2))) = wave_all.(fnames{i}).STLocation_pka(:,2);    
    PKALocThetawST(i,length(wave_all.(fnames{i}).STLocation_pka(:,2))+1:longestosc) = NaN;

    ControlLocPhiwST(i,1:length(wave_all.(fnames{i}).STLocation(:,3))) = wave_all.(fnames{i}).STLocation(:,3);
    ControlLocPhiwST(i,length(wave_all.(fnames{i}).STLocation(:,3))+1:longestosc) = NaN;
    PKALocPhiwST(i,1:length(wave_all.(fnames{i}).STLocation_pka(:,3))) = wave_all.(fnames{i}).STLocation_pka(:,3);    
    PKALocPhiwST(i,length(wave_all.(fnames{i}).STLocation_pka(:,3))+1:longestosc) = NaN;
end

writematrix([ControlLocRadiuswST, PKALocRadiuswST], [savename, 'waveRadiusStabilityST.csv'])
writematrix([ControlLocThetawST, PKALocThetawST], [savename, 'waveThetaStabilityST.csv'])
writematrix([ControlLocPhiwST, PKALocPhiwST], [savename, 'wavePhiStabilityST.csv'])


%look at "lenghth" traveled:
locchange = [];
for i = 1:12
locchange = [locchange sqrt((mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10,1)) - mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10_pka,1)))^2 + ...
    (mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10,2)) - mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10_pka,2)))^2 + ...
    (mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10,3)) - mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10_pka,3)))^2)];
end
writematrix(locchange, [savename 'locchange_wave_PKA.csv'])

locchangenet = [];
for i = 1:12
locchangenet = [locchangenet, sqrt((mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells,1)) - mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells_pka,1)))^2 + ...
    (mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells,2)) - mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells_pka,2)))^2 + ...
    (mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells,3)) - mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells_pka,3)))^2)];
end
writematrix(locchangenet, [savename 'locchange_net_PKA.csv'])

% look at correlation
figure, 
for i = 1:length(fnames)
    hold on
    plot([network_all.(fnames{i}).correlation, network_all.(fnames{i}).correlation_pka], 'linewidth',3)
end