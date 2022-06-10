%% Network Analysis of 3D data
%code for network analysis of 3D data
% Jennifer Briggs 06/02/2022


close all 
clear all
clc
%changing defaults
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)

addpath('~/Documents/GitHub/Functional_and_Structural_Networks')


fileloc = ('/Volumes/Briggs_10TB/Merrin/Ca_Courses/Singlecelltraces/EJ106/')
filename = 'F03 10G'
fullfile = [fileloc filename]
%Here we load the files
    calcium = readmatrix([fullfile '_Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
    calcium(1:3,:) = []; %the CSV you have has the first 3 rows as NAN so we remove them
    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
    Locations = readmatrix([fullfile ' pos_Detailed.csv']);
    Locations = Locations(:,1:3);

    
%find best threshold
Threshold = findoptRth(calcium);
%% Run Network Analysis

Opts.printasSubplot = 0;
Opts.Subplotnum = 0;
Opts.figs = 1;
Opts.multiseed = 0;
Opts.multiseednum = 1;




[N, adj, kperc, histArrayPercShort, Rij] = links(calcium, 0.999,Opts,1); %%This is where the network is built
[sorted, cellsor]= sort(N);


%network analysis is performed
[L, EGlob, CClosed, ELocClosed, COpen, ELocOpen, nopath]  = graphProperties(adj);

G = graph(adj);
figure, plot(G)

Hubs = cellsor((sorted - min(sorted))/range(sorted)>.60);

figure, plot(graph(sparse(adj)), 'Xdata',Locations(:,1), 'Ydata',Locations(:,2), 'Zdata',Locations(:,3),'NodeColor','k', 'MarkerSize',10)