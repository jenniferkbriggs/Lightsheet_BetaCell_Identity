
%% Analyze network and wave from 3D

%extract data
timeseries = fieldnames(wave_all)
GKa = find(contains(timeseries, 'GKa'));
PKa = find(contains(timeseries, 'PKa'));
Control = setxor([1:length(timeseries)], [GKa; PKa]);
Control(5) = [];
GKa(1) = [];
%GKa(end-1:end) = [];
% PKa(3) = [];

savename  = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange/'


% Correlation

    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).correlation));
        GKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).correlation_pka));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).correlation));
        PKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).correlation_pka));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).correlation));
        Control_corrafter(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).correlation_pka))
    end

        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'Correlation', 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'Correlation', 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'Correlation', 'PKa.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

    %Frequency

    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).freq));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).freq_pka));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).freq));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).freq_pka));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).freq));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).freq_pka))
    end
 
        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'Freq', 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'Freq', 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'Freq', 'PKa.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


    %Duty Cycle
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).duty_cycle));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).duty_cycle_pka));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).duty_cycle));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).duty_cycle_pka));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).duty_cycle));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).duty_cycle_pka))
    end
 
        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'DutyCycle', 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'DutyCycle', 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'DutyCycle', 'PKa.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


function [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays()
PKa_corrbefore = NaN(1,14);
PKa_corrafter = NaN(1,14);
PKa_corrboth = NaN(1,14);
GKa_corrbefore = NaN(1,14);
GKa_corrafter = NaN(1,14);
GKa_corrboth = NaN(1,14);
Control_corrbefore = NaN(1,14);
Control_corrafter = NaN(1,14);
Control_corrboth = NaN(1,14);
Control_corrboth = NaN(1,14);
PKa_corrboth = NaN(1,14);
GKa_corrboth = NaN(1,14);
end
