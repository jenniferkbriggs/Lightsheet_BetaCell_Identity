%% Network Analysis of 3D data
%code for network analysis of 3D data
% Jennifer Briggs 06/02/2022

%changing defaults
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)

fileloc = ('/Volumes/Briggs_10TB/Merrin/Ca_Courses/')
filename = 'F03 10G'
fullfile = [fileloc filename]
%Here we load the files
    calcium = readmatrix([fullfile '_Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
    calcium(1:3,:) = []; %the CSV you have has the first 3 rows as NAN so we remove them
    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
    Locations = readmatrix([fullfile ' pos_Detailed.csv']);


%% Run Network Analysis

Opts.printasSubplot = 1;
Opts.Subplotnum = 1;
Opts.figs = 0;
Opts.multiseed = 0;
Opts.multiseednum = 1;

Threshold = 0.999; %will need to play with this

[N, adj, kperc, histArrayPercShort, Rij] = links(calcium, Threshold,Opts,1); %%This is where the network is built
[sorted, cellsor]= sort(N);


%network analysis is performed
[L, EGlob, CClosed, ELocClosed, COpen, ELocOpen, nopath]  = graphProperties(adj);



