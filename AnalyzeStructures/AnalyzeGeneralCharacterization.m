%% Analyze network and wave from 3D

% Extract data field names from 'wave_all' structure
timeseries = fieldnames(wave_all);

% Identify indices for GKa and PKa conditions in timeseries
GKa = find(contains(timeseries, 'GKa'));
PKa = find(contains(timeseries, 'PKa'));

% Identify control group by excluding GKa and PKa
Control = setxor(1:length(timeseries), [GKa; PKa]);

% Remove specific entries from Control and GKa groups
Control(5) = []; % Removes the 5th entry from Control
GKa(1) = [];     % Removes the 1st entry from GKa

% Define directory to save analysis results
savename = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange/';

%% Correlation Analysis

% Compute average non-zero correlations for GKa group before and after treatment
for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).correlation));
    GKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).correlation_pka));
end

% Compute average non-zero correlations for PKa group before and after treatment
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).correlation));
    PKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).correlation_pka));
end

% Compute average non-zero correlations for Control group before and after treatment
for i = 1:length(Control)
    Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).correlation));
    Control_corrafter(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).correlation_pka));
end

% Save correlation results to CSV files
CorrelationTable = array2table([Control_corrbefore; Control_corrafter]', 'VariableNames', {'Pre-', 'Vehicle'});
writetable(CorrelationTable, [savename, 'Correlation', 'Control.csv']);

CorrelationTable = array2table([GKa_corrbefore; GKa_corrafter]', 'VariableNames', {'Pre-', 'GKa'});
writetable(CorrelationTable, [savename, 'Correlation', 'GKa.csv']);

CorrelationTable = array2table([PKa_corrbefore; PKa_corrafter]', 'VariableNames', {'Pre-', 'PKa'});
writetable(CorrelationTable, [savename, 'Correlation', 'PKa.csv']);

% Initialize arrays for further analysis
[Control_corrafter, Control_corrbefore, GKa_corrbefore, GKa_corrafter, PKa_corrbefore, PKa_corrafter, Control_corrboth, PKa_corrboth, GKa_corrboth] = mkarrays();

%% Frequency Analysis

% Compute average non-zero frequencies for GKa group before and after treatment
for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).freq));
    GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).freq_pka));
end

% Compute average non-zero frequencies for PKa group before and after treatment
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).freq));
    PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).freq_pka));
end

% Compute average non-zero frequencies for Control group before and after treatment
for i = 1:length(Control)
    Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).freq));
    Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).freq_pka));
end

% Save frequency results to CSV files
CorrelationTable = array2table([Control_corrbefore; Control_corrafter]', 'VariableNames', {'Pre-', 'Vehicle'});
writetable(CorrelationTable, [savename, 'Freq', 'Control.csv']);

CorrelationTable = array2table([GKa_corrbefore; GKa_corrafter]', 'VariableNames', {'Pre-', 'GKa'});
writetable(CorrelationTable, [savename, 'Freq', 'GKa.csv']);

CorrelationTable = array2table([PKa_corrbefore; PKa_corrafter]', 'VariableNames', {'Pre-', 'PKa'});
writetable(CorrelationTable, [savename, 'Freq', 'PKa.csv']);

[Control_corrafter, Control_corrbefore, GKa_corrbefore, GKa_corrafter, PKa_corrbefore, PKa_corrafter, Control_corrboth, PKa_corrboth, GKa_corrboth] = mkarrays();

%% Duty Cycle Analysis

% Compute average non-zero duty cycles for GKa group before and after treatment
for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).duty_cycle));
    GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).duty_cycle_pka));
end

% Compute average non-zero duty cycles for PKa group before and after treatment
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).duty_cycle));
    PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).duty_cycle_pka));
end

% Compute average non-zero duty cycles for Control group before and after treatment
for i = 1:length(Control)
    Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).duty_cycle));
    Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).duty_cycle_pka));
end

% Save duty cycle results to CSV files
CorrelationTable = array2table([Control_corrbefore; Control_corrafter]', 'VariableNames', {'Pre-', 'Vehicle'});
writetable(CorrelationTable, [savename, 'DutyCycle', 'Control.csv']);

CorrelationTable = array2table([GKa_corrbefore; GKa_corrafter]', 'VariableNames', {'Pre-', 'GKa'});
writetable(CorrelationTable, [savename, 'DutyCycle', 'GKa.csv']);

CorrelationTable = array2table([PKa_corrbefore; PKa_corrafter]', 'VariableNames', {'Pre-', 'PKa'});
writetable(CorrelationTable, [savename, 'DutyCycle', 'PKa.csv']);

[Control_corrafter, Control_corrbefore, GKa_corrbefore, GKa_corrafter, PKa_corrbefore, PKa_corrafter, Control_corrboth, PKa_corrboth, GKa_corrboth] = mkarrays();

%% Helper Function: Initialize Arrays
function [Control_corrafter, Control_corrbefore, GKa_corrbefore, GKa_corrafter, PKa_corrbefore, PKa_corrafter, Control_corrboth, PKa_corrboth, GKa_corrboth] = mkarrays()
    % Initialize NaN arrays for storing computed values
    PKa_corrbefore = NaN(1,14);
    PKa_corrafter = NaN(1,14);
    PKa_corrboth = NaN(1,14);
    GKa_corrbefore = NaN(1,14);
    GKa_corrafter = NaN(1,14);
    GKa_corrboth = NaN(1,14);
    Control_corrbefore = NaN(1,14);
    Control_corrafter = NaN(1,14);
    Control_corrboth = NaN(1,14);
end
