%% This script verifies the correctness of oscillations generated in RunNetworkandWave
% Author: Jennifer Briggs
% Date: 06/02/2022

close all          % Close all open figure windows
% clear all        % (Commented out) Would clear all variables from workspace
clc                % Clear command window

% Set default figure properties for consistent visualization
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)

% Add necessary code directories to the MATLAB path
addpath('~/Documents/GitHub/Functional_and_Structural_Networks')
try
    addpath('/Users/brigjenn/Documents/GitHub/UniversalCode')
end

%% Define experimental conditions and data directory
pka = 1;  % Set to 1 if PKA is present, 0 for control condition
fileloc = '/Volumes/Briggs_10TB/Merrin/EJ155 New data set/';  % Path to experiment folders

% Retrieve list of experiment directories, ignoring hidden/system files
files = dir(fileloc);
files_dir = {};
for i = 1:length(files)
    if isfolder([fileloc files(i).name]) && ~contains(files(i).name, '.')
        files_dir{end+1} = {files(i).name};
    end
end

%% Begin analysis for each experiment
for i = 2:length(files_dir)  % Start from index 2 to skip '.' and '..'
    fullfile = [fileloc files_dir{i}{1} '/'];  % Full path to current experiment folder
    name = files_dir{i}{1};                   % Extract folder name

    % Define save location for analysis results
    savename = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange';

    %% Load calcium signal data
    calcium = readmatrix([fullfile name(1:end-10) 'Plot.csv']);  % Load calcium signal data
    time = calcium(:,1);         % Extract time from the first column
    calcium(:,1) = [];           % Remove time column, leaving only calcium intensity data

    photobleaching = 0;          % Flag for photobleaching correction (0 = no correction)
    calstore = calcium;          % Backup of calcium data
    timestore = time;            % Backup of time data
    fignum = 1;                  % Initialize figure counter

    % Load pre-identified oscillation indices
    load([fullfile 'waveindx_samenumber.mat']); 

    %% Plot calcium data with repolarization oscillation markers
    figure
    plot(mean(calcium'))  % Plot average calcium signal across cells
    xline(cuttime, 'linewidth', 3, 'color', 'b')  % Mark cuttime

    % Overlay start and end markers for repolarization oscillations
    for i = 1:length(start_indx_repol)
        hold on
        xline(start_indx_repol(i), 'label', ['start: ' num2str(i)])
        xline(end_indx_repol(i), 'color', 'r', 'label', ['end: ' num2str(i)])
    end
    title(name)
    keyboard  % Pause execution for interactive inspection

    %% Plot calcium data with full oscillation markers
    figure
    plot(mean(calcium'))  % Plot average calcium signal
    xline(cuttime, 'linewidth', 3, 'color', 'b')  % Mark cuttime

    % Overlay start and end markers for full oscillations
    for i = 1:length(start_indx)
        hold on
        xline(start_indx(i), 'label', ['start: ' num2str(i)])
        xline(end_indx(i), 'color', 'r', 'label', ['end: ' num2str(i)])
    end
    keyboard  % Pause for inspection

    % Save updated oscillation indices back to file
    save([fullfile, 'waveindx_samenumber.mat'], 'start_indx', 'end_indx', 'start_indx_repol', 'end_indx_repol', 'cuttime', 'numtrial')

    close all  
