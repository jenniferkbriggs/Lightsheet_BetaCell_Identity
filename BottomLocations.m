%% Just get distance from center of bottom cells: 

cd ~/'OneDrive - The University of Colorado Denver'/Anschutz/Islet/3DLightSheet/
load FinalAnalysis_04.05.2024.mat


%extract data
timeseries = fieldnames(wave_all)
GKa = find(contains(timeseries, 'GKa'));
PKa = find(contains(timeseries, 'PKa'));
Control = setxor([1:length(timeseries)], [GKa; PKa]);
Control(5) = [];
GKa(1) = [];
% PKa(3) = [];

savename  = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange/'


%this is what I want to extract: 

Varnames ={'Control Before','Control After', 'GKa Before', 'GKa After', 'PKa Before','PKa After'}

% PKa_corrbefore = NaN(1,14);
% PKa_corrafter = NaN(1,14);
% PKa_corrboth = NaN(1,14);
% GKa_corrbefore = NaN(1,14);
% GKa_corrafter = NaN(1,14);
% GKa_corrboth = NaN(1,14);
% Control_corrbefore = NaN(1,14);
% Control_corrafter = NaN(1,14);
% Control_corrboth = NaN(1,14);


%% location top 10

    %find average distance from center:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).dist_from_center));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).dist_from_center));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).dist_from_center));
    end
    average = [Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore'];

    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).distcenter_top10));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).distcenter_top10));

        bottomcells_w = wave_all.(timeseries{GKa(i)}).bottom10;
        bottomcells_n = network_all.(timeseries{GKa(i)}).bottom10cells;

        GKa_wlow(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).dist_from_center(bottomcells_w)));
        GKa_nlow(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).dist_from_center(nonzeros(bottomcells_n))));

    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).distcenter_top10));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).distcenter_top10));

        bottomcells_w = wave_all.(timeseries{PKa(i)}).bottom10;
        bottomcells_n = network_all.(timeseries{PKa(i)}).bottom10cells;

        PKa_wlow(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).dist_from_center(bottomcells_w)));
        PKa_nlow(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).dist_from_center(nonzeros(bottomcells_n))));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).distcenter_top10));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).distcenter_top10));

        bottomcells_w = wave_all.(timeseries{Control(i)}).bottom10;
        bottomcells_n = network_all.(timeseries{Control(i)}).bottom10cells;

        Control_wlow(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).dist_from_center(bottomcells_w)));
        Control_nlow(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).dist_from_center(nonzeros(bottomcells_n))));
    end
 

    Network_loc = [Control_corrbefore, GKa_corrbefore, PKa_corrbefore]';
    Wave_loc = [Control_corrafter,  GKa_corrafter, PKa_corrafter]';

    Network_low = [Control_nlow,  GKa_nlow, PKa_nlow]';
    Wave_low = [Control_wlow,  GKa_wlow, PKa_wlow]';
   % Wave_locend = [Control_corrboth,  GKa_corrboth, PKa_corrboth]';

    CorrelationTable = array2table([average, Network_loc, Network_low], 'VariableNames', {'Average Cell','High Degree', '0 Degree'})
        writetable(CorrelationTable, [savename, 'NetworkvsAverageLoc_withlow.csv'] )

    CorrelationTable = array2table([average, Wave_loc, Wave_low], 'VariableNames', {'Average Cell','Early Phase', 'Late Phase'})
        writetable(CorrelationTable, [savename, 'WavevsAverageLoc_withlow.csv'] )


    % CorrelationTable = array2table([average, Wave_locend], 'VariableNames', {'Average Cell','High Phase End'})
    %     writetable(CorrelationTable, [savename, 'WavevsAverageLoc_end.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

   