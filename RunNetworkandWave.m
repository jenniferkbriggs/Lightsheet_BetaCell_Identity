%% Analysis of 3D data
%code for network analysis of 3D data
% Jennifer Briggs 06/02/2022
% 
% 
close all 
%clear all
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
for i = 1:length(files_dir)
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
end

for i = 1:14
fullfile = [fileloc files_dir{i}{1} '/']
name = files_dir{i}{1}

%%Here we load the files
calcium = readmatrix([fullfile name(1:end-10) 'Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
time = calcium(:,1);   %time is in the first column so pull this out;
calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
Locations = readmatrix([fullfile name(1:end-10) 'Detailed.csv']);
   
Locations = Locations(:,1:3);
cut = find(isnan(Locations(1,:)))
Locations(cut, :) = [];
    calstore = calcium;
    timestore = time;
load([fullfile 'waveindx_samenumber.mat']);

%% 09j
%close figure 7
%%Run network analysisclose all

network = Network_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, name, Locations, timestore,start_indx, end_indx)
%%Run waveintiator analysis
wave = Wave_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, name, Locations, timestore,start_indx, end_indx, time,1, 1)

%remove spaces from string:
name = name(find(~isspace(name)));
%remove dashes:
toremove = strfind(name, '-');
name(toremove) = [];

wave_all.(name)= wave;
network_all.(name)= network;
wave_all.(name).Locations=  normalize(Locations, "range");
calcium_all.(name) = calcium;
time_all.(name) = time;
disp(i)

end


%save all results to matlab
save([savename 'control_PKAanalyses.mat'], 'network_all','wave_all', 'time_all', 'calcium_all')

%% 

%% Write results to csv. 

if pka
fnames = fieldnames(wave_all);
wave_cell = array2table([mean(nonzeros(wave_all.(fnames{1}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{1}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{2}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{2}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{3}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{3}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{4}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{4}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{5}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{5}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{6}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{6}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{7}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{7}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{8}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{8}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{9}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{9}).Degree_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{10}).Degree_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{10}).Degree_KL_wave(2).trial))], 'VariableNames',{'Pre','Post'})

writetable(wave_cell, [savename 'controlWaveCell.csv'])


network_cell = array2table([mean(nonzeros(network_all.(fnames{1}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{1}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{2}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{2}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{3}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{3}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{4}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{4}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{5}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{5}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{6}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{6}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{7}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{7}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{8}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{8}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{9}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{9}).Degree_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{10}).Degree_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{10}).Degree_KL_net(2).trial))], 'VariableNames',{'Pre','Post'})
writetable(network_cell, [savename 'controlNetCell.csv'])


wave_region = array2table([mean(nonzeros(wave_all.(fnames{1}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{1}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{2}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{2}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{3}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{3}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{4}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{4}).COG_KL_wave(2).trial)) 
mean(nonzeros(wave_all.(fnames{5}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{5}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{6}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{6}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{7}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{7}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{8}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{8}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{9}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{9}).COG_KL_wave(2).trial))
mean(nonzeros(wave_all.(fnames{10}).COG_KL_wave(1).trial)), mean(nonzeros(wave_all.(fnames{10}).COG_KL_wave(2).trial))], 'VariableNames',{'Pre','Post'})
writetable(wave_region, [savename 'controlWaveRegion.csv'])

network_region = array2table([mean(nonzeros(network_all.(fnames{1}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{1}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{2}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{2}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{3}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{3}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{4}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{4}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{5}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{5}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{6}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{6}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{7}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{7}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{8}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{8}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{9}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{9}).COG_KL_net(2).trial))
mean(nonzeros(network_all.(fnames{10}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{10}).COG_KL_net(2).trial))], 'VariableNames',{'Pre','Post'})
writetable(network_region, [savename 'controlNetRegion.csv'])



netregion = [mean(nonzeros([network_all.(fnames{1}).COG_KL_net(1).trial; network_all.(fnames{1}).COG_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{2}).COG_KL_net(1).trial; network_all.(fnames{2}).COG_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{3}).COG_KL_net(1).trial; network_all.(fnames{3}).COG_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{4}).COG_KL_net(1).trial; network_all.(fnames{4}).COG_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{5}).COG_KL_net(1).trial; network_all.(fnames{5}).COG_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{6}).COG_KL_net(1).trial; network_all.(fnames{6}).COG_KL_net(2).trial])),...
    mean(nonzeros([network_all.(fnames{7}).COG_KL_net(1).trial; network_all.(fnames{7}).COG_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{8}).COG_KL_net(1).trial; network_all.(fnames{8}).COG_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{9}).COG_KL_net(1).trial; network_all.(fnames{9}).COG_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{10}).COG_KL_net(1).trial; network_all.(fnames{10}).COG_KL_net(2).trial]))];
%    mean(nonzeros([network_all.(fnames{12}).COG_KL_net(1).trial; network_all.(fnames{12}).COG_KL_net(2).trial]))];
waveregion = [mean(nonzeros([wave_all.(fnames{1}).COG_KL_wave(1).trial; wave_all.(fnames{1}).COG_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{2}).COG_KL_wave(1).trial; wave_all.(fnames{2}).COG_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{3}).COG_KL_wave(1).trial; wave_all.(fnames{3}).COG_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{4}).COG_KL_wave(1).trial; wave_all.(fnames{4}).COG_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{5}).COG_KL_wave(1).trial; wave_all.(fnames{5}).COG_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{6}).COG_KL_wave(1).trial; wave_all.(fnames{6}).COG_KL_wave(2).trial])),...
    mean(nonzeros([wave_all.(fnames{7}).COG_KL_wave(1).trial; wave_all.(fnames{7}).COG_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{8}).COG_KL_wave(1).trial; wave_all.(fnames{8}).COG_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{9}).COG_KL_wave(1).trial; wave_all.(fnames{9}).COG_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{10}).COG_KL_wave(1).trial; wave_all.(fnames{10}).COG_KL_wave(2).trial]))];
%    mean(nonzeros([wave_all.(fnames{12}).COG_KL_wave(1).trial; wave_all.(fnames{12}).COG_KL_wave(2).trial]))];
netcell = [mean(nonzeros([network_all.(fnames{1}).Degree_KL_net(1).trial; network_all.(fnames{1}).Degree_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{2}).Degree_KL_net(1).trial; network_all.(fnames{2}).Degree_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{3}).Degree_KL_net(1).trial; network_all.(fnames{3}).Degree_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{4}).Degree_KL_net(1).trial; network_all.(fnames{4}).Degree_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{5}).Degree_KL_net(1).trial; network_all.(fnames{5}).Degree_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{6}).Degree_KL_net(1).trial; network_all.(fnames{6}).Degree_KL_net(2).trial])),...
    mean(nonzeros([network_all.(fnames{7}).Degree_KL_net(1).trial; network_all.(fnames{7}).Degree_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{8}).Degree_KL_net(1).trial; network_all.(fnames{8}).Degree_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{9}).Degree_KL_net(1).trial; network_all.(fnames{9}).Degree_KL_net(2).trial])), ...
    mean(nonzeros([network_all.(fnames{10}).Degree_KL_net(1).trial; network_all.(fnames{10}).Degree_KL_net(2).trial]))];
%    mean(nonzeros([network_all.(fnames{12}).Degree_KL_net(1).trial; network_all.(fnames{12}).Degree_KL_net(2).trial]))];
wavecell =[mean(nonzeros([wave_all.(fnames{1}).Degree_KL_wave(1).trial; wave_all.(fnames{1}).Degree_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{2}).Degree_KL_wave(1).trial; wave_all.(fnames{2}).Degree_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{3}).Degree_KL_wave(1).trial; wave_all.(fnames{3}).Degree_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{4}).Degree_KL_wave(1).trial; wave_all.(fnames{4}).Degree_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{5}).Degree_KL_wave(1).trial; wave_all.(fnames{5}).Degree_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{6}).Degree_KL_wave(1).trial; wave_all.(fnames{6}).Degree_KL_wave(2).trial])),...
    mean(nonzeros([wave_all.(fnames{7}).Degree_KL_wave(1).trial; wave_all.(fnames{7}).Degree_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{8}).Degree_KL_wave(1).trial; wave_all.(fnames{8}).Degree_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{9}).Degree_KL_wave(1).trial; wave_all.(fnames{9}).Degree_KL_wave(2).trial])), ...
    mean(nonzeros([wave_all.(fnames{10}).Degree_KL_wave(1).trial; wave_all.(fnames{10}).Degree_KL_wave(2).trial]))];
%    mean(nonzeros([wave_all.(fnames{12}).Degree_KL_wave(1).trial; wave_all.(fnames{12}).Degree_KL_wave(2).trial]))];


%control - cellular vs. region:
cellvsregiont = array2table([waveregion, netregion; wavecell, netcell]); 
writetable(cellvsregiont, [savename 'cellvsregional.csv'])

writematrix([waveregion',wavecell'], [savename 'waveregionvscell.csv'])
writematrix([netregion',netcell'], [savename 'netregionvscell.csv'])



maintainednet = array2table([network_all.(fnames{1}).mainted_perc;
    network_all.(fnames{2}).mainted_perc; network_all.(fnames{3}).mainted_perc;
    network_all.(fnames{4}).mainted_perc; network_all.(fnames{5}).mainted_perc;
    network_all.(fnames{6}).mainted_perc; network_all.(fnames{7}).mainted_perc;
    network_all.(fnames{8}).mainted_perc; network_all.(fnames{9}).mainted_perc; 
    network_all.(fnames{10}).mainted_perc]);
    %network_all.(fnames{12}).mainted_perc])
writetable(wave_cell, [savename 'controlNetmaintained.csv'])


varnames = {'Pre- Intravariability','Post- Intravariability', 'Intervariability'}

pkavar = array2table([wave_all.(fnames{1}).intravar, wave_all.(fnames{1}).intravar_pka, wave_all.(fnames{1}).intervar; 
    wave_all.(fnames{2}).intravar, wave_all.(fnames{2}).intravar_pka, wave_all.(fnames{2}).intervar; 
    wave_all.(fnames{3}).intravar, wave_all.(fnames{3}).intravar_pka, wave_all.(fnames{3}).intervar; 
    wave_all.(fnames{4}).intravar, wave_all.(fnames{4}).intravar_pka, wave_all.(fnames{4}).intervar; 
    wave_all.(fnames{5}).intravar, wave_all.(fnames{5}).intravar_pka, wave_all.(fnames{5}).intervar;
    wave_all.(fnames{6}).intravar, wave_all.(fnames{6}).intravar_pka, wave_all.(fnames{6}).intervar;
    wave_all.(fnames{7}).intravar, wave_all.(fnames{7}).intravar_pka, wave_all.(fnames{7}).intervar; 
    wave_all.(fnames{8}).intravar, wave_all.(fnames{8}).intravar_pka, wave_all.(fnames{8}).intervar; 
    wave_all.(fnames{9}).intravar, wave_all.(fnames{9}).intravar_pka, wave_all.(fnames{9}).intervar; 
    wave_all.(fnames{10}).intravar, wave_all.(fnames{10}).intravar_pka, wave_all.(fnames{10}).intervar],'VariableNames',{'Pre- Intravariability','Post- Intravariability',  'Intervaraibility'})
writetable(pkavar, [savename 'controlvariability.csv'])

Netpkavar = array2table([network_all.(fnames{1}).intravar, network_all.(fnames{1}).intravar_pka, network_all.(fnames{1}).intervar; 
    network_all.(fnames{2}).intravar, network_all.(fnames{2}).intravar_pka, network_all.(fnames{2}).intervar; 
    network_all.(fnames{3}).intravar, network_all.(fnames{3}).intravar_pka, network_all.(fnames{3}).intervar; 
    network_all.(fnames{4}).intravar, network_all.(fnames{4}).intravar_pka, network_all.(fnames{4}).intervar; 
    network_all.(fnames{5}).intravar, network_all.(fnames{5}).intravar_pka, network_all.(fnames{5}).intervar;
    network_all.(fnames{6}).intravar, network_all.(fnames{6}).intravar_pka, network_all.(fnames{6}).intervar;
    network_all.(fnames{7}).intravar, network_all.(fnames{7}).intravar_pka, network_all.(fnames{7}).intervar; 
    network_all.(fnames{8}).intravar, network_all.(fnames{8}).intravar_pka, network_all.(fnames{8}).intervar; 
    network_all.(fnames{9}).intravar, network_all.(fnames{9}).intravar_pka, network_all.(fnames{9}).intervar; 
    network_all.(fnames{10}).intravar, network_all.(fnames{10}).intravar_pka, network_all.(fnames{10}).intervar],... 
    'VariableNames',{'Pre- Intravariability','Post- Intravariability',  'Intervaraibility'}) 
writetable(Netpkavar, [savename 'Network_controlvariability.csv'])



for i = 1:length(fnames)
    numosc = length(wave_all.(fnames{i}).phaserange)/2
    pre(i) = mean(wave_all.(fnames{i}).phaserange(1:numosc));
    post(i) = mean(wave_all.(fnames{i}).phaserange(numosc+1:end));
end

phaseranget = array2table([pre', post'])
writetable(phaseranget, [savename 'controlRange_of_phases.csv'])

locchange = [];
for i = 1:10
locchange = [locchange sqrt((mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10,1)) - mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10_pka,1)))^2 + ...
    (mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10,2)) - mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10_pka,2)))^2 + ...
    (mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10,3)) - mean(wave_all.(fnames{i}).Locations(wave_all.(fnames{i}).top10_pka,3)))^2)];
end
writematrix(locchange, [savename 'locchange_wave_control.csv'])



locchangenet = [];
for i = 1:10
locchangenet = [locchangenet, sqrt((mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells,1)) - mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells_pka,1)))^2 + ...
    (mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells,2)) - mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells_pka,2)))^2 + ...
    (mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells,3)) - mean(wave_all.(fnames{i}).Locations(network_all.(fnames{i}).top10cells_pka,3)))^2)];
end
writematrix(locchangenet, [savename 'locchange_net_control.csv'])


else
    fnames = fieldnames(wave_all);
    wave = array2table([mean(nonzeros(wave_all.(fnames{1}).COG_KL_wave.trial)), mean(nonzeros(wave_all.(fnames{1}).Degree_KL_wave.trial))
    mean(nonzeros(wave_all.(fnames{2}).COG_KL_wave.trial)), mean(nonzeros(wave_all.(fnames{2}).Degree_KL_wave.trial))
    mean(nonzeros(wave_all.(fnames{3}).COG_KL_wave.trial)), mean(nonzeros(wave_all.(fnames{3}).Degree_KL_wave.trial))
    mean(nonzeros(wave_all.(fnames{4}).COG_KL_wave.trial)), mean(nonzeros(wave_all.(fnames{4}).Degree_KL_wave.trial))], 'VariableNames',{'Regional','Cellular'})
    writetable(wave, [savename 'control_long_Wave_regionalvscellular.csv'])
    
    
    network = array2table([mean(nonzeros(network_all.(fnames{1}).COG_KL_net(1).trial)), mean(nonzeros(network_all.(fnames{1}).Degree_KL_net.trial))
    mean(nonzeros(network_all.(fnames{2}).COG_KL_net.trial)), mean(nonzeros(network_all.(fnames{2}).Degree_KL_net.trial))
    mean(nonzeros(network_all.(fnames{3}).COG_KL_net.trial)), mean(nonzeros(network_all.(fnames{3}).Degree_KL_net.trial))
    mean(nonzeros(network_all.(fnames{4}).COG_KL_net.trial)), mean(nonzeros(network_all.(fnames{4}).Degree_KL_net.trial))], 'VariableNames',{'Regional','Cellular'})
    writetable(network, [savename 'control_long_Net_regionalvscellular.csv'])

end
