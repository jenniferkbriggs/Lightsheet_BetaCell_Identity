%% Run script for analyzing most of the 3D consistency data

close all  % Close all open figures
clc        % Clear the command window

% Set default figure properties for visualization
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)

% Add necessary directories to the MATLAB path
addpath(genpath('~/Documents/GitHub/Functional_and_Structural_Networks'))
try
    addpath('/Users/brigjenn/Documents/GitHub/UniversalCode')
end

% Save current date for file naming
savetime = datestr(datetime('today'),'yyyymmdd');

%% User-defined parameters
pka = 1; % 1 if PKA is present, 0 if control
fileloc = '/Volumes/Briggs_10TB/Merrin/EJ155 New data set/'; % Directory with experimental data

% Collect folder names of experiments
files = dir(fileloc);
files_dir = {};
for i = 5:length(files)
    if isfolder([fileloc files(i).name]) && ~contains(files(i).name, '.')
        files_dir{end+1} = {files(i).name};
    end
end

%% Start analysis loop for each experiment
for i = 2:length(files_dir)
    fullfile = [fileloc files_dir{i}{1} '/'];
    name = files_dir{i}{1};
    
    % Define save path
    savename = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange';
    
    %% Load data files
    try
        Locations = readmatrix([fullfile name(1:end-10) 'detailed.csv']);
    catch
        Locations = readmatrix([fullfile name(1:end-10) 'Detailed.csv']);
    end
    Locations = Locations(:,1:3); % Extract XYZ coordinates
    calcium = readmatrix([fullfile name(1:end-10) 'Plot.csv']); % Load calcium intensity data
    
    time = calcium(:,1);  % Extract time vector
    calcium(:,1) = [];    % Remove time column from calcium data
    
    photobleaching = 0;   % Flag for photobleaching correction
    calstore = calcium;
    timestore = time;
    fignum = 1;
    
    %% Identify oscillations
    try
        load([fullfile 'waveindx_samenumber.mat']);
    catch
        if pka
            figure, plot(mean(calcium'));
            title('Select start of second phase, PKA administration, and end of recording');
            admin = find(abs(time - 1800)<0.2);
            xline(admin(1));
            [cuttime, ~] = ginput(3);
            close;
        end
        
        numtrial = pka + 1;
        cuttime(1) = 1;
        cuttime(2) = length(calcium);
        
        figure, plot(time, mean(calcium'));
        xline(time(round(cuttime(2))), 'linewidth',3);
        
        %Here we identify the oscillations
        [start_indx, end_indx] = identify_oscillations(calcium, time, 0);
        if ~exist('end_indx_repol')
            [start_indx_repol, end_indx_repol] = identify_oscillations_end(end_indx, calcium, time);
        end
        save([fullfile 'waveindx_samenumber.mat'], 'start_indx','end_indx', 'start_indx_repol', 'end_indx_repol','cuttime','numtrial');
    end
    
    %% Perform 3D Wave and Network Analysis
    % Run code to assess 3D wave consistency
    wave = Wave_3D(calcium, cuttime, numtrial, photobleaching, [savename '/AllImages'], fileloc, name, Locations, timestore, start_indx, end_indx, time,1, 1,1);

    % Run code to assess 3D wave consistency
    network = Network_3D(calcium, cuttime, numtrial, photobleaching, [savename '/AllImages'], fileloc, name, Locations, timestore, start_indx, 0);
    
    %% Clean experiment name for structure field
    name = regexprep(name, '[\s-]', '');
    network_all.(name) = network;
    disp(i);
    clear end_indx_repol;
end

%% Save the analysis results
% save([savename, savetime, 'tenperc_all' 'Analysis', '.mat'], 'network_all');
