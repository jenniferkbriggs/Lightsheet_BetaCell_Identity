%% Analyze network and wave from 3D

%extract data
timeseries = fieldnames(network_all)
GKa = find(contains(timeseries, 'GKa'));
PKa = find(contains(timeseries, 'PKa'));
Control = setxor([1:length(timeseries)], [GKa; PKa]);
Control(5) = [];
GKa(1) = [];
% PKa(3) = [];

savename  = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange/'


%this is what I want to extract: 

Varnames ={'Control Before','Control After', 'GKa Before', 'GKa After', 'PKa Before','PKa After'}

PKa_corrbefore = NaN(1,14);
PKa_corrafter = NaN(1,14);
PKa_corrboth = NaN(1,14);
GKa_corrbefore = NaN(1,14);
GKa_corrafter = NaN(1,14);
GKa_corrboth = NaN(1,14);
Control_corrbefore = NaN(1,14);
Control_corrafter = NaN(1,14);
Control_corrboth = NaN(1,14);

%look at the location spread of the waves:
for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(wave_all.(timeseries{GKa(i)}).STLocation);
    GKa_corrafter(i) = mean(wave_all.(timeseries{GKa(i)}).STLocation_pka);
end
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(wave_all.(timeseries{PKa(i)}).STLocation);
    PKa_corrafter(i) = mean(wave_all.(timeseries{PKa(i)}).STLocation_pka);
end
for i = 1:length(Control)
    Control_corrbefore(i) = mean(wave_all.(timeseries{Control(i)}).STLocation);
    Control_corrafter(i) = mean(wave_all.(timeseries{Control(i)}).STLocation_pka);
end
        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'LocationSpread' , 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'LocationSpread' , 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'LocationSpread' ,'PKa.csv'] ) 
%reset arrays to nan
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();



%look at the location spread of the waves:
for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(network_all.(timeseries{GKa(i)}).STLocation);
    GKa_corrafter(i) = mean(network_all.(timeseries{GKa(i)}).STLocation_pka);
end
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(network_all.(timeseries{PKa(i)}).STLocation);
    PKa_corrafter(i) = mean(network_all.(timeseries{PKa(i)}).STLocation_pka);
end
for i = 1:length(Control)
    Control_corrbefore(i) = mean(network_all.(timeseries{Control(i)}).STLocation);
    Control_corrafter(i) = mean(network_all.(timeseries{Control(i)}).STLocation_pka);
end
        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'LocationSpread_Network' , 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'LocationSpread_Network' , 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'LocationSpread_Network' ,'PKa.csv'] ) 
%reset arrays to nan
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


%compare location spread between network and wave: 

%look at the location spread of the waves:
for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(wave_all.(timeseries{GKa(i)}).STLocation);
    GKa_corrafter(i) = mean(network_all.(timeseries{GKa(i)}).STLocation);
end
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(wave_all.(timeseries{PKa(i)}).STLocation);
    PKa_corrafter(i) = mean(network_all.(timeseries{PKa(i)}).STLocation);
end
for i = 1:length(Control)
    Control_corrbefore(i) = mean(wave_all.(timeseries{Control(i)}).STLocation);
    Control_corrafter(i) = mean(network_all.(timeseries{Control(i)}).STLocation);
end
       

    CorrelationTable = array2table([Control_corrbefore, GKa_corrbefore, PKa_corrbefore; ...
        Control_corrafter,  GKa_corrafter, PKa_corrafter]', 'VariableNames', {'Wave', 'Network'})
    writetable(CorrelationTable, [savename, 'NetworkVsPhase_locationSpread.csv'] )

%reset arrays to nan
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();
% 

% %look at how to COG moves
% for i = 1:length(GKa)
%     GKa_corrbefore(i) = mean(mean(wave_all.(timeseries{GKa(i)}).cog_mov, 'omitnan'), 'omitnan');
%     GKa_corrafter(i) = mean(mean(network_all.(timeseries{GKa(i)}).cog_mov, 'omitnan'), 'omitnan');
% end
% for i = 1:length(PKa)
%     PKa_corrbefore(i) = mean(mean(wave_all.(timeseries{PKa(i)}).cog_mov, 'omitnan'), 'omitnan');
%     PKa_corrafter(i) = mean(mean(network_all.(timeseries{PKa(i)}).cog_mov, 'omitnan'), 'omitnan');
% end
% for i = 1:length(Control)
%     Control_corrbefore(i) = mean(mean(wave_all.(timeseries{Control(i)}).cog_mov, 'omitnan'), 'omitnan');
%     Control_corrafter(i) = mean(mean(network_all.(timeseries{Control(i)}).cog_mov, 'omitnan'), 'omitnan');
% end

    % 
    % CorrelationTable = array2table([Control_corrbefore, GKa_corrbefore, PKa_corrbefore; ...
    %     Control_corrafter,  GKa_corrafter, PKa_corrafter]', 'VariableNames', {'Wave', 'Network'})
    % writetable(CorrelationTable, [savename, 'NetworkVsPhase_locationMov.csv'] )

%reset arrays to nan
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


%% location top 10

    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).distcenter_top10));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).distcenter_top10_pka));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).distcenter_top10));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).distcenter_top10_pka));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).distcenter_top10));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).distcenter_top10_pka))
    end
 
        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'Wave_distancetop10', 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'Wave_distancetop10', 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'Wave_distancetop10', 'PKa.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

   

    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).distcenter_top10));
        GKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).distcenter_top10_pka));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).distcenter_top10));
        PKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).distcenter_top10_pka));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).distcenter_top10));
        Control_corrafter(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).distcenter_top10_pka))
    end
 
        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'Network_distancetop10', 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'Network_distancetop10', 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'Network_distancetop10', 'PKa.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

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

    %wave velocity
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).Velocity));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).Velocity_pka));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).Velocity));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).Velocity_pka));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).Velocity));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).Velocity_pka))
    end
 
        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'Velocity', 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'Velocity', 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'Velocity', 'PKa.csv'] )
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


    %Phase spread
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).phasespread_10));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).phasespread_10_pka));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).phasespread_10));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).phasespread_10_pka));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).phasespread_10));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).phasespread_10_pka))
    end
 
        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'PhaseSpread', 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'PhaseSpread', 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'PhaseSpread', 'PKa.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


%no intervention comparing network and wave

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
%        GKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).distcenter_top10));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).distcenter_top10));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).distcenter_top10));
        %PKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).distcenter_top10));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).distcenter_top10));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).distcenter_top10));
       % Control_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).distcenter_top10))
    end
 

    CorrelationTable = array2table([Control_corrbefore, GKa_corrbefore, PKa_corrbefore; ...
        Control_corrafter,  GKa_corrafter, PKa_corrafter]', 'VariableNames', {'Network', 'Wave'})
    writetable(CorrelationTable, [savename, 'NetworkWaveLocation.csv'] )

    Network_loc = [Control_corrbefore, GKa_corrbefore, PKa_corrbefore]';
    Wave_loc = [Control_corrafter,  GKa_corrafter, PKa_corrafter]';
   % Wave_locend = [Control_corrboth,  GKa_corrboth, PKa_corrboth]';

    CorrelationTable = array2table([average, Network_loc], 'VariableNames', {'Average Cell','High Degree'})
        writetable(CorrelationTable, [savename, 'NetworkvsAverageLoc.csv'] )

    CorrelationTable = array2table([average, Wave_loc], 'VariableNames', {'Average Cell','High Phase'})
        writetable(CorrelationTable, [savename, 'WavevsAverageLoc.csv'] )


    % CorrelationTable = array2table([average, Wave_locend], 'VariableNames', {'Average Cell','High Phase End'})
    %     writetable(CorrelationTable, [savename, 'WavevsAverageLoc_end.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


        %only 10%
    for i = 1:length(GKa)
        GKa_corrbefore(i) = ((network_all.(timeseries{GKa(i)}).top10percent_avg));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = ((network_all.(timeseries{PKa(i)}).top10percent_avg));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = ((network_all.(timeseries{Control(i)}).top10percent_avg));
    end
    average = [Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore'];

    for i = 1:length(GKa)
        GKa_corrbefore(i) = ((network_all.(timeseries{GKa(i)}).top10percent_avg));
        GKa_corrafter(i) = ((wave_all.(timeseries{GKa(i)}).top10percent_avg));
      %  GKa_corrafter(i) = ((wave_all_end.(timeseries{GKa(i)}).top10percent_avg));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = ((network_all.(timeseries{PKa(i)}).top10percent_avg));
        PKa_corrafter(i) = ((wave_all.(timeseries{PKa(i)}).top10percent_avg));
      %  PKa_corrboth(i) = ((wave_all_end.(timeseries{PKa(i)}).top10percent_avg));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = ((network_all.(timeseries{Control(i)}).top10percent_avg));
        Control_corrafter(i) = ((wave_all.(timeseries{Control(i)}).top10percent_avg));
       % Control_corrboth(i) = wave_all_end.(timeseries{Control(i)}).top10percent_avg;
    end
 

    %Waveloc_end = [Control_corrboth; GKa_corrafter; PKa_corrboth]';
    Network_loc = [Control_corrbefore; GKa_corrbefore; PKa_corrbefore]';
    Wave_loc = [Control_corrafter;  GKa_corrafter; PKa_corrafter]';
    CorrelationTable = array2table([Network_loc], 'VariableNames', {'Control','GKa','PKa'})
        writetable(CorrelationTable, [savename, 'Network_top10percentchange.csv'] )

    CorrelationTable = array2table([Wave_loc], 'VariableNames',  {'Control','GKa','PKa'})
        writetable(CorrelationTable, [savename, 'Wave_top10percentchange.csv'] )

    % CorrelationTable = array2table([Waveloc_end], 'VariableNames',  {'Control','GKa','PKa'})
    %     writetable(CorrelationTable, [savename, 'WaveEnd_top10percentchange.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

%% KL Consistency %% 

% Wave ------- 
Variables = {'COG_KL_wave', 'Degree_KL_wave'}
for j = 1:length(Variables)
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).(Variables{j})(1).trial));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).(Variables{j})(2).trial));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).(Variables{j})(1).trial));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).(Variables{j})(2).trial));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).(Variables{j})(1).trial));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).(Variables{j})(2).trial))
    end
    
    CorrelationTable = array2table([Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
    Control_corrafter', GKa_corrafter', PKa_corrafter']')
    
    writetable(CorrelationTable, [savename, Variables{j}, '.csv'] )

        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, Variables{j}, 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, Variables{j}, 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, Variables{j}, 'PKa.csv'] )
end
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


%KL Network!!
Variables = {'COG_KL_net', 'Degree_KL_net'}
for j = 1:length(Variables)
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).(Variables{j})(1).trial));
        GKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).(Variables{j})(2).trial));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).(Variables{j})(1).trial));
        PKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).(Variables{j})(2).trial));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).(Variables{j})(1).trial));
        Control_corrafter(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).(Variables{j})(2).trial))
    end
    
    CorrelationTable = array2table([Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
    Control_corrafter', GKa_corrafter', PKa_corrafter']');
    writetable(CorrelationTable, [savename, Variables{j}, '.csv'] )
        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, Variables{j}, 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, Variables{j}, 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, Variables{j}, 'PKa.csv'] )
end
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();



% Cell vs region
Variables = {'COG_KL_wave', 'Degree_KL_wave'}
Variables2 = {'COG_KL_net', 'Degree_KL_net'}

GKa_corrafter = [GKa_corrboth; GKa_corrboth];
GKa_corrbefore = [GKa_corrboth; GKa_corrboth];
Control_corrboth = [GKa_corrboth; GKa_corrboth];
Control_corrafter = [GKa_corrboth; GKa_corrboth];
Control_corrbefore = [GKa_corrboth; GKa_corrboth];
PKa_corrboth = [GKa_corrboth; GKa_corrboth];
PKa_corrafter = [GKa_corrboth; GKa_corrboth];
PKa_corrbefore = [GKa_corrboth; GKa_corrboth];
GKa_corrboth = [GKa_corrboth; GKa_corrboth];


for j = 1:length(Variables)
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).(Variables{j})(1).trial));
        GKa_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).(Variables2{j})(1).trial));
%       GKa_corrboth(j,i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).(Variables{j})(1).trial));

    end
    for i = 1:length(PKa)
        PKa_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).(Variables{j})(1).trial));
        PKa_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).(Variables2{j})(1).trial));
         %PKa_corrboth(j,i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).(Variables{j})(1).trial));
    end
    for i = 1:length(Control)
        Control_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).(Variables{j})(1).trial));
        Control_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{Control(i)}).(Variables2{j})(1).trial));
       % Control_corrboth(j,i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).(Variables{j})(1).trial));
    end
    
     
end
   CorrelationTable = array2table(fliplr([Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore']), 'VariableNames',{'Cellular', 'Regional'})
    writetable(CorrelationTable, [savename, 'CellvsRegionPhase.csv'] )

   CorrelationTable = array2table(fliplr([Control_corrafter'; GKa_corrafter'; PKa_corrafter']), 'VariableNames',{'Cellular', 'Regional'})
    writetable(CorrelationTable, [savename, 'CellvsRegionNetwork.csv'] )   
    
    % CorrelationTable = array2table(fliplr([Control_corrboth'; GKa_corrboth'; PKa_corrboth']), 'VariableNames',{'Cellular', 'Regional'})
    % writetable(CorrelationTable, [savename, 'CellvsRegionWaveEnd.csv'] )


[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();



    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = (network_all.(timeseries{GKa(i)}).change);
        GKa_corrafter(i) = ((wave_all.(timeseries{GKa(i)}).change));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = ((network_all.(timeseries{PKa(i)}).change));
        PKa_corrafter(i) = ((wave_all.(timeseries{PKa(i)}).change));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = (network_all.(timeseries{Control(i)}).change);
        Control_corrafter(i) = (wave_all.(timeseries{Control(i)}).change)
    end
    
    CorrelationTable = array2table([Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'], 'VariableNames', {'Control','GKa','PKa'});
    writetable(CorrelationTable, [savename, 'NetworkChange.csv'] )

    CorrelationTable = array2table([Control_corrafter',  GKa_corrafter', PKa_corrafter'], 'VariableNames', {'Control','GKa','PKa'});
    writetable(CorrelationTable, [savename, 'WaveChange.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();




%KL Phase!!
Variables = {'COG_KL_wave', 'Degree_KL_wave'}
for j = 1:length(Variables)
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).(Variables{j})(1).trial));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).(Variables{j})(2).trial));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).(Variables{j})(1).trial));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).(Variables{j})(2).trial));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).(Variables{j})(1).trial));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).(Variables{j})(2).trial))
    end

    CorrelationTable = array2table([Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
    Control_corrafter', GKa_corrafter', PKa_corrafter']')

    writetable(CorrelationTable, [savename, Variables{j}, '.csv'] )

        CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, Variables{j}, 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, Variables{j}, 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, Variables{j}, 'PKa.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

end

% 
% 
% %KL Phase!!
% Variables = {'COG_KL_wave', 'Degree_KL_wave'}
% for j = 1:length(Variables)
%     % Get regional consistency:
%     for i = 1:length(GKa)
%         GKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).(Variables{j})(1).trial));
%         GKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).(Variables{j})(2).trial));
%     end
%     for i = 1:length(PKa)
%         PKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).(Variables{j})(1).trial));
%         PKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).(Variables{j})(2).trial));
%     end
%     for i = 1:length(Control)
%         Control_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).(Variables{j})(1).trial));
%         Control_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).(Variables{j})(2).trial))
%     end
% 
%     CorrelationTable = array2table([Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
%     Control_corrafter', GKa_corrafter', PKa_corrafter']')
% 
%     writetable(CorrelationTable, [savename, Variables{j}, '.csv'] )
% 
%         CorrelationTable = array2table([Control_corrbefore; ...
%         Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
%     writetable(CorrelationTable, [savename, Variables{j}, 'depolarization_Control.csv'] )
%         CorrelationTable = array2table([GKa_corrbefore; ...
%         GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
%     writetable(CorrelationTable, [savename, Variables{j}, 'depolarization_GKa.csv'] )
%         CorrelationTable = array2table([PKa_corrbefore; ...
%         PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
%     writetable(CorrelationTable, [savename, Variables{j}, 'depolarization_PKa.csv'] )
% [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();
% end

%KL Phase vs Network!!
Variables = {'COG_KL_', 'Degree_KL_'}
for j = 1:length(Variables)
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).([Variables{j} 'net'])(1).trial));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).([Variables{j} 'wave'])(1).trial));
%        GKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).([Variables{j} 'wave'])(1).trial));

    end
    
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).([Variables{j} 'net'])(1).trial));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).([Variables{j} 'wave'])(1).trial));
      %  PKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).([Variables{j} 'wave'])(1).trial));

    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).([Variables{j} 'net'])(1).trial));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).([Variables{j} 'wave'])(1).trial));
      %  Control_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).([Variables{j} 'wave'])(1).trial))

    end
    
    CorrelationTable = array2table([Control_corrbefore, GKa_corrbefore, PKa_corrbefore; ...
        Control_corrafter,  GKa_corrafter, PKa_corrafter; ...
        Control_corrboth,  GKa_corrboth, PKa_corrboth]', 'VariableNames', {'Network', 'Wave','Wave End'})
    writetable(CorrelationTable, [savename, Variables{j}, 'NetworkVsDegree.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

end


%Difference in KL degree and COG
Variables = {'COG_KL_', 'Degree_KL_'}
for j = 1:length(Variables)
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(j,i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).([Variables{j} 'net'])(1).trial));
        GKa_corrafter(j,i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).([Variables{j} 'wave'])(1).trial));
        %GKa_corrboth(j,i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).([Variables{j} 'wave'])(1).trial));

    end
    for i = 1:length(PKa)
        PKa_corrbefore(j,i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).([Variables{j} 'net'])(1).trial));
        PKa_corrafter(j,i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).([Variables{j} 'wave'])(1).trial));
      %  PKa_corrboth(j,i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).([Variables{j} 'wave'])(1).trial));

    end
    for i = 1:length(Control)
        Control_corrbefore(j,i) = mean(nonzeros(network_all.(timeseries{Control(i)}).([Variables{j} 'net'])(1).trial));
        Control_corrafter(j,i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).([Variables{j} 'wave'])(1).trial));
    %    Control_corrboth(j,i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).([Variables{j} 'wave'])(1).trial))

    end
end
myarray = [(Control_corrbefore(2,:)-Control_corrbefore(1,:)), (GKa_corrbefore(2,:)-GKa_corrbefore(1,:)), (PKa_corrbefore(2,:)-PKa_corrbefore(1,:)); ...
        (Control_corrafter(2,:)-Control_corrafter(1,:)),  (GKa_corrafter(2,:)-GKa_corrafter(1,:)), (PKa_corrafter(2,:)-PKa_corrafter(1,:))]';
       % (Control_corrboth(2,:)-Control_corrboth(1,:)),  (GKa_corrboth(2,:)-GKa_corrboth(1,:)), (PKa_corrboth(2,:)-PKa_corrboth(1,:))]';
    CorrelationTable = array2table(myarray./max(myarray), 'VariableNames', {'Network', 'Wave'})
    writetable(CorrelationTable, [savename, Variables{j}, 'KL_COGminusdegree.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();
%% --- Lowest 10%
% Intravaraibility Network
VarNames = {'Pre- intravariability', 'Vehicle intravariability', 'Intervariability'}
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).intravar_bottom));
        GKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).intravarbottom_pka));
        GKa_corrboth(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).intervarbottom));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).intravar_bottom));
        PKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).intravarbottom_pka));
        PKa_corrboth(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).intervarbottom));

    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).intravar_bottom));
        Control_corrafter(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).intravarbottom_pka));
        Control_corrboth(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).intervarbottom));

    end
    
    CorrelationTable = array2table([Control_corrafter', GKa_corrafter', PKa_corrafter'; ...
        Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
        Control_corrboth', GKa_corrboth', PKa_corrboth']')
    writetable(CorrelationTable, [savename,'networkIntravariability_bottom.csv'] )

    CorrelationTable = array2table([Control_corrbefore; 
        Control_corrafter; 
        Control_corrboth]', 'VariableNames',VarNames)
    writetable(CorrelationTable, [savename,'Net_ControlIntravariability_bottom.csv'] )

    CorrelationTable = array2table([GKa_corrbefore; 
        GKa_corrafter; 
        GKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'GKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'Net_GKaIntravariability_bottom.csv'] )

        CorrelationTable = array2table([PKa_corrbefore; 
        PKa_corrafter; 
        PKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'PKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'Net_PKaIntravariability_bottom.csv'] )

[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();
 

% Intravaraibility Wave
VarNames = {'Pre- intravariability', 'Vehicle intravariability', 'Intervariability'}
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).intravar_bottom));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).intravar_bottom_pka));
        GKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).intervar_bottom));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).intravar_bottom));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).intravar_bottom_pka));
        PKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).intervar_bottom));

    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).intravar_bottom));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).intravar_bottom_pka));
        Control_corrboth(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).intervar_bottom));

    end
    
    CorrelationTable = array2table([Control_corrafter', GKa_corrafter', PKa_corrafter'; ...
        Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
        Control_corrboth', GKa_corrboth', PKa_corrboth']')
    writetable(CorrelationTable, [savename,'wave_Intravariabilitybottom.csv'] )

    CorrelationTable = array2table([Control_corrbefore; 
        Control_corrafter; 
        Control_corrboth]', 'VariableNames',VarNames)
    writetable(CorrelationTable, [savename,'wave_ControlIntravariabilitybottom.csv'] )

    CorrelationTable = array2table([GKa_corrbefore; 
        GKa_corrafter; 
        GKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'GKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'wave_GKaIntravariabilitybottom.csv'] )

        CorrelationTable = array2table([PKa_corrbefore; 
        PKa_corrafter; 
        PKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'PKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'wave_PKaIntravariabilitybottom.csv'] )

[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

 
% 
% % Intravaraibility Wave
% VarNames = {'Pre- intravariability', 'Vehicle intravariability', 'Intervariability'}
%     % Get regional consistency:
%     for i = 1:length(GKa)
%         GKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).intravar_bottom));
%         GKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).intravar_bottom_pka));
%         GKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).intervar_bottom));
%     end
%     for i = 1:length(PKa)
%         PKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).intravar_bottom));
%         PKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).intravar_bottom_pka));
%         PKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).intervar_bottom));
% 
%     end
%     for i = 1:length(Control)
%         Control_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).intravar_bottom));
%         Control_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).intravar_bottom_pka));
%         Control_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).intervar_bottom));

%     end
% 
%     CorrelationTable = array2table([Control_corrafter', GKa_corrafter', PKa_corrafter'; ...
%         Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
%         Control_corrboth', GKa_corrboth', PKa_corrboth']')
%     writetable(CorrelationTable, [savename,'Repolarization_Intravariabilitybottom.csv'] )
% 
%     CorrelationTable = array2table([Control_corrbefore; 
%         Control_corrafter; 
%         Control_corrboth]', 'VariableNames',VarNames)
%     writetable(CorrelationTable, [savename,'Repolarization_ControlIntravariabilitybottom.csv'] )
% 
%     CorrelationTable = array2table([GKa_corrbefore; 
%         GKa_corrafter; 
%         GKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'GKa intravariability', 'Intervariability'})
%     writetable(CorrelationTable, [savename,'Repolarization_GKaIntravariabilitybottom.csv'] )
% 
%         CorrelationTable = array2table([PKa_corrbefore; 
%         PKa_corrafter; 
%         PKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'PKa intravariability', 'Intervariability'})
%     writetable(CorrelationTable, [savename,'Repolarization_PKaIntravariabilitybottom.csv'] )
% 
% [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


%%  Network vs Wave - control
VarNames = {'Pre- intravariability', 'Vehicle intravariability', 'Intervariability'}
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).intravar_bottom));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).intravar_bottom));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).intravar_bottom));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).intravar_bottom));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).intravar_bottom));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).intravar_bottom));
    end
    
    CorrelationTable = array2table([Control_corrbefore', Control_corrafter'; GKa_corrbefore', GKa_corrafter'; PKa_corrbefore', PKa_corrafter'], 'VariableNames', {'Network', 'Wave'})
    writetable(CorrelationTable, [savename,'networkVsWaveIntravariability_bottom.csv'] )


        % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).intravar));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).intravar));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).intravar));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).intravar));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).intravar));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).intravar));
    end
    
    CorrelationTable = array2table([Control_corrbefore', Control_corrafter'; GKa_corrbefore', GKa_corrafter'; PKa_corrbefore', PKa_corrafter'], 'VariableNames', {'Network', 'Wave'})
    writetable(CorrelationTable, [savename,'networkVsWaveIntravariability.csv'] )


%% --- Calcium Wave Axis %% 
% Intravaraibility Network
VarNames = {'Pre- intravariability', 'Vehicle intravariability', 'Intervariability'}
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).rot));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).rot_pka));
        GKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).rot_intra));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).rot));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).rot_pka));
        PKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).rot_intra));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).rot));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).rot_pka));
        Control_corrboth(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).rot_intra));
    end

    CorrelationTable = array2table([Control_corrbefore; 
        Control_corrafter; 
        Control_corrboth]', 'VariableNames',VarNames)
    writetable(CorrelationTable, [savename,'CalciumWave_control.csv'] )

    CorrelationTable = array2table([GKa_corrbefore; 
        GKa_corrafter; 
        GKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'GKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'CalciumWave_GKal.csv'] )

        CorrelationTable = array2table([PKa_corrbefore; 
        PKa_corrafter; 
        PKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'PKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'CalciumWave_Pka.csv'] )

[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

% Just compare Before and After
VarNames = {'Pre- intravariability', 'Vehicle intravariability'}
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).rot));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).rot_pka));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).rot));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).rot_pka));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).rot));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).rot_pka));
    end

    %the normalization is wrong - multiply by 3 again. 
    CorrelationTable = array2table([Control_corrbefore; 
        Control_corrafter]', 'VariableNames',VarNames)
    writetable(CorrelationTable, [savename,'CalciumWave_control_nointer.csv'] )

    CorrelationTable = array2table([GKa_corrbefore; 
        GKa_corrafter]', 'VariableNames',{'Pre- intravariability', 'GKa intravariability'})
    writetable(CorrelationTable, [savename,'CalciumWave_GKa_nointer.csv'] )

        CorrelationTable = array2table([PKa_corrbefore; 
        PKa_corrafter]', 'VariableNames',{'Pre- intravariability', 'PKa intravariability'})
    writetable(CorrelationTable, [savename,'CalciumWave_Pka_nointer.csv'] )

[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();
% 

% 
% for i = 1:length(GKa)
%         GKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).rot));
%         GKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).rot_pka));
%         GKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).rot_intra));
%     end
%     for i = 1:length(PKa)
%         PKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).rot));
%         PKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).rot_pka));
%         PKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).rot_intra));
%     end
%     for i = 1:length(Control)
%         Control_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).rot));
%         Control_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).rot_pka));
%         Control_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).rot_intra));
%     end
% 
%     CorrelationTable = array2table([Control_corrbefore; 
%         Control_corrafter; 
%         Control_corrboth]', 'VariableNames',VarNames)
%     writetable(CorrelationTable, [savename,'Depol_CalciumWave_control.csv'] )
% 
%     CorrelationTable = array2table([GKa_corrbefore; 
%         GKa_corrafter; 
%         GKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'GKa intravariability', 'Intervariability'})
%     writetable(CorrelationTable, [savename,'Depol_CalciumWave_GKal.csv'] )
% 
%         CorrelationTable = array2table([PKa_corrbefore; 
%         PKa_corrafter; 
%         PKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'PKa intravariability', 'Intervariability'})
%     writetable(CorrelationTable, [savename,'Depol_CalciumWave_Pka.csv'] )
% 
% [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


%% --- Top 10%
% Intravaraibility Network
VarNames = {'Pre- intravariability', 'Vehicle intravariability', 'Intervariability'}
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).intravar));
        GKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).intravar_pka));
        GKa_corrboth(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).intervar));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).intravar));
        PKa_corrafter(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).intravar_pka));
        PKa_corrboth(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).intervar));

    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).intravar));
        Control_corrafter(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).intravar_pka));
        Control_corrboth(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).intervar));

    end
    
    CorrelationTable = array2table([Control_corrafter', GKa_corrafter', PKa_corrafter'; ...
        Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
        Control_corrboth', GKa_corrboth', PKa_corrboth']')
    writetable(CorrelationTable, [savename,'networkIntravariability.csv'] )

    CorrelationTable = array2table([Control_corrbefore; 
        Control_corrafter; 
        Control_corrboth]', 'VariableNames',VarNames)
    writetable(CorrelationTable, [savename,'Net_ControlIntravariability.csv'] )

    CorrelationTable = array2table([GKa_corrbefore; 
        GKa_corrafter; 
        GKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'GKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'Net_GKaIntravariabilitycsv'] )

        CorrelationTable = array2table([PKa_corrbefore; 
        PKa_corrafter; 
        PKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'PKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'Net_PKaIntravariability.csv'] )

[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


% Intravaraibility Wave
VarNames = {'Pre- intravariability', 'Vehicle intravariability', 'Intervariability'}
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).intravar));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).intravar_pka));
        GKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).intervar));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).intravar));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).intravar_pka));
        PKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).intervar));

    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).intravar));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).intravar_pka));
        Control_corrboth(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).intervar));

    end
    
    CorrelationTable = array2table([Control_corrafter', GKa_corrafter', PKa_corrafter'; ...
        Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
        Control_corrboth', GKa_corrboth', PKa_corrboth']')
    writetable(CorrelationTable, [savename,'Phase_Intravariability.csv'] )

    CorrelationTable = array2table([Control_corrbefore; 
        Control_corrafter; 
        Control_corrboth]', 'VariableNames',VarNames)
    writetable(CorrelationTable, [savename,'Phase_ControlIntravariability.csv'] )

    CorrelationTable = array2table([GKa_corrbefore; 
        GKa_corrafter; 
        GKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'GKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'Phase_GKaIntravariability.csv'] )

        CorrelationTable = array2table([PKa_corrbefore; 
        PKa_corrafter; 
        PKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'PKa intravariability', 'Intervariability'})
    writetable(CorrelationTable, [savename,'Phase_PKaIntravariability.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();



% 
% % Intravaraibility Wave
% VarNames = {'Pre- intravariability', 'Vehicle intravariability', 'Intervariability'}
%     % Get regional consistency:
%     for i = 1:length(GKa)
%         GKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).intravar));
%         GKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).intravar_pka));
%         GKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).intervar));
%     end
%     for i = 1:length(PKa)
%         PKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).intravar));
%         PKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).intravar_pka));
%         PKa_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).intervar));
% 
%     end
%     for i = 1:length(Control)
%         Control_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).intravar));
%         Control_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).intravar_pka));
%         Control_corrboth(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).intervar));
% 
%     end
% 
%     CorrelationTable = array2table([Control_corrafter', GKa_corrafter', PKa_corrafter'; ...
%         Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
%         Control_corrboth', GKa_corrboth', PKa_corrboth']')
%     writetable(CorrelationTable, [savename,'Depol_Phase_Intravariability.csv'] )
% 
%     CorrelationTable = array2table([Control_corrbefore; 
%         Control_corrafter; 
%         Control_corrboth]', 'VariableNames',VarNames)
%     writetable(CorrelationTable, [savename,'Depol_Phase_ControlIntravariability.csv'] )
% 
%     CorrelationTable = array2table([GKa_corrbefore; 
%         GKa_corrafter; 
%         GKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'GKa intravariability', 'Intervariability'})
%     writetable(CorrelationTable, [savename,'Depol_Phase_GKaIntravariability.csv'] )
% 
%         CorrelationTable = array2table([PKa_corrbefore; 
%         PKa_corrafter; 
%         PKa_corrboth]', 'VariableNames',{'Pre- intravariability', 'PKa intravariability', 'Intervariability'})
%     writetable(CorrelationTable, [savename,'Depol_Phase_PKaIntravariability.csv'] )
% [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


%not using

Variables = {'r','el','ez'};

for k = 1:3
% Get lcoation:
for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(network_all.(timeseries{GKa(i)}).meanLocation(:,k));
    GKa_corrafter(i) = mean(network_all.(timeseries{GKa(i)}).meanLocation_pka(:,k));
end
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(network_all.(timeseries{PKa(i)}).meanLocation(:,k));
    PKa_corrafter(i) = mean(network_all.(timeseries{PKa(i)}).meanLocation_pka(:,k));
end
for i = 1:length(Control)
    Control_corrbefore(i) = mean(network_all.(timeseries{Control(i)}).meanLocation(:,k));
    Control_corrafter(i) = mean(network_all.(timeseries{Control(i)}).meanLocation_pka(:,k));
end
CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'NetworkLoc_' Variables{k}, 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'NetworkLoc_' Variables{k}, 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'NetworkLoc_' Variables{k}, 'PKa.csv'] )

    [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

end



for k = 1:3
% Get lcoation:
for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(wave_all.(timeseries{GKa(i)}).MeanLocation(:,k));
    GKa_corrafter(i) = mean(network_all.(timeseries{GKa(i)}).meanLocation(:,k));
end
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(wave_all.(timeseries{PKa(i)}).MeanLocation(:,k));
    PKa_corrafter(i) = mean(network_all.(timeseries{PKa(i)}).meanLocation(:,k));
end
for i = 1:length(Control)
    Control_corrbefore(i) = mean(wave_all.(timeseries{Control(i)}).MeanLocation(:,k));
    Control_corrafter(i) = mean(network_all.(timeseries{Control(i)}).meanLocation(:,k));
end
CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'WaveLoc_' Variables{k}, 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'WaveLoc_' Variables{k}, 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'WaveLoc_' Variables{k}, 'PKa.csv'] )
    [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

end


%location wave vs network



for k = 1:3
% Get lcoation:
for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(wave_all.(timeseries{GKa(i)}).MeanLocation(:,k));
    GKa_corrafter(i) = mean(wave_all.(timeseries{GKa(i)}).meanLocation_pka(:,k));
end
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(wave_all.(timeseries{PKa(i)}).MeanLocation(:,k));
    PKa_corrafter(i) = mean(wave_all.(timeseries{PKa(i)}).meanLocation_pka(:,k));
end
for i = 1:length(Control)
    Control_corrbefore(i) = mean(wave_all.(timeseries{Control(i)}).MeanLocation(:,k));
    Control_corrafter(i) = mean(wave_all.(timeseries{Control(i)}).meanLocation_pka(:,k));
end
  CorrelationTable = array2table([Control_corrbefore, GKa_corrbefore, PKa_corrbefore; ...
        Control_corrafter,  GKa_corrafter, PKa_corrafter]', 'VariableNames', {'Wave', 'Network'})
    writetable(CorrelationTable, [savename,  Variables{k}, 'NetworkVsPhase_location.csv'] )
    [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

end



for i = 1:length(GKa)
    GKa_corrbefore(i) = mean(wave_all.(timeseries{GKa(i)}).MeanLocation(:,k));
    GKa_corrafter(i) = mean(network_all.(timeseries{GKa(i)}).meanLocation(:,k));
end
for i = 1:length(PKa)
    PKa_corrbefore(i) = mean(wave_all.(timeseries{PKa(i)}).MeanLocation(:,k));
    PKa_corrafter(i) = mean(network_all.(timeseries{PKa(i)}).meanLocation(:,k));
end
for i = 1:length(Control)
    Control_corrbefore(i) = mean(wave_all.(timeseries{Control(i)}).MeanLocation(:,k));
    Control_corrafter(i) = mean(network_all.(timeseries{Control(i)}).meanLocation(:,k));
end
CorrelationTable = array2table([Control_corrbefore; ...
        Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
    writetable(CorrelationTable, [savename, 'WaveLoc_' Variables{k}, 'Control.csv'] )
        CorrelationTable = array2table([GKa_corrbefore; ...
        GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
    writetable(CorrelationTable, [savename, 'WaveLoc_' Variables{k}, 'GKa.csv'] )
        CorrelationTable = array2table([PKa_corrbefore; ...
        PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
    writetable(CorrelationTable, [savename, 'WaveLoc_' Variables{k}, 'PKa.csv'] )
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

%% - wave end attic (now we call it depol and repol -- 
%  for i = 1:length(GKa)
%         GKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).distcenter_top10));
%         GKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).distcenter_top10_pka));
%     end
%     for i = 1:length(PKa)
%         PKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).distcenter_top10));
%         PKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).distcenter_top10_pka));
%     end
%     for i = 1:length(Control)
%         Control_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).distcenter_top10));
%         Control_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).distcenter_top10_pka))
%     end
% 
%         CorrelationTable = array2table([Control_corrbefore; ...
%         Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
%     writetable(CorrelationTable, [savename, 'Wave_distancetop10', 'depol_Control.csv'] )
%         CorrelationTable = array2table([GKa_corrbefore; ...
%         GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
%     writetable(CorrelationTable, [savename, 'Wave_distancetop10', 'depol_GKa.csv'] )
%         CorrelationTable = array2table([PKa_corrbefore; ...
%         PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
%     writetable(CorrelationTable, [savename, 'Wave_distancetop10', 'depol_PKa.csv'] )
% [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

% 
% 
% %look at the location spread of the waves:
% for i = 1:length(GKa)
%     GKa_corrbefore(i) = mean(wave_all_end.(timeseries{GKa(i)}).loc_spread);
%     GKa_corrafter(i) = mean(wave_all_end.(timeseries{GKa(i)}).loc_spread_pka);
% end
% for i = 1:length(PKa)
%     PKa_corrbefore(i) = mean(wave_all_end.(timeseries{PKa(i)}).loc_spread);
%     PKa_corrafter(i) = mean(wave_all_end.(timeseries{PKa(i)}).loc_spread_pka);
% end
% for i = 1:length(Control)
%     Control_corrbefore(i) = mean(wave_all_end.(timeseries{Control(i)}).loc_spread);
%     Control_corrafter(i) = mean(wave_all_end.(timeseries{Control(i)}).loc_spread_pka);
% end
%         CorrelationTable = array2table([Control_corrbefore; ...
%         Control_corrbefore]', 'VariableNames',{'Pre-', 'Vehicle'})
%     writetable(CorrelationTable, [savename, 'LocationSpread' , 'depol_Control.csv'] )
%         CorrelationTable = array2table([GKa_corrbefore; ...
%         GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
%     writetable(CorrelationTable, [savename, 'LocationSpread' , 'depol_GKa.csv'] )
%         CorrelationTable = array2table([PKa_corrbefore; ...
%         PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
%     writetable(CorrelationTable, [savename, 'LocationSpread' ,'depol_PKa.csv'] ) 
% %reset arrays to nan
% [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


%     %wave velocity
%     for i = 1:length(GKa)
%         GKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).Velocity));
%         GKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).Velocity_pka));
%     end
%     for i = 1:length(PKa)
%         PKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).Velocity));
%         PKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).Velocity_pka));
%     end
%     for i = 1:length(Control)
%         Control_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).Velocity));
%         Control_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).Velocity_pka))
%     end
% 
%         CorrelationTable = array2table([Control_corrbefore; ...
%         Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
%     writetable(CorrelationTable, [savename, 'Velocity', 'Depol_Control.csv'] )
%         CorrelationTable = array2table([GKa_corrbefore; ...
%         GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
%     writetable(CorrelationTable, [savename, 'Velocity', 'Depol_GKa.csv'] )
%         CorrelationTable = array2table([PKa_corrbefore; ...
%         PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
%     writetable(CorrelationTable, [savename, 'Velocity', 'Depol_PKa.csv'] )
% [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();
% 
% 
% %Phase spread
%     for i = 1:length(GKa)
%         GKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).phasespread_10));
%         GKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).phasespread_10_pka));
%     end
%     for i = 1:length(PKa)
%         PKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).phasespread_10));
%         PKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).phasespread_10_pka));
%     end
%     for i = 1:length(Control)
%         Control_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).phasespread_10));
%         Control_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).phasespread_10_pka))
%     end
% 
%         CorrelationTable = array2table([Control_corrbefore; ...
%         Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
%     writetable(CorrelationTable, [savename, 'PhaseSpread', 'Depol_Control.csv'] )
%         CorrelationTable = array2table([GKa_corrbefore; ...
%         GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
%     writetable(CorrelationTable, [savename, 'PhaseSpread', 'Depol_GKa.csv'] )
%         CorrelationTable = array2table([PKa_corrbefore; ...
%         PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
%     writetable(CorrelationTable, [savename, 'PhaseSpread', 'Depol_PKa.csv'] )
% [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();
%% is the top10 of the wave the bottom 10 of the wave depolarization?
% 
%  for i = 1:length(GKa)
%         GKa_corrbefore(i) = length(intersect(depolarize_all.(timeseries{GKa(i)}).bottom10, repolarize_all.(timeseries{GKa(i)}).bottom10))./length(repolarize_all.(timeseries{GKa(i)}).bottom10);
%         GKa_corrafter(i) = length(intersect(depolarize_all.(timeseries{GKa(i)}).top10av, repolarize_all.(timeseries{GKa(i)}).top10av))./length(repolarize_all.(timeseries{GKa(i)}).top10av);
%  end
% 
%   for i = 1:length(PKa)
%         PKa_corrbefore(i) = length(intersect(depolarize_all.(timeseries{PKa(i)}).bottom10, repolarize_all.(timeseries{PKa(i)}).bottom10))./length(repolarize_all.(timeseries{PKa(i)}).bottom10);
%         PKa_corrafter(i) = length(intersect(depolarize_all.(timeseries{PKa(i)}).top10av, repolarize_all.(timeseries{PKa(i)}).top10av))./length(repolarize_all.(timeseries{PKa(i)}).top10av);
%   end
% 
%    for i = 1:length(Control)
%         Control_corrbefore(i) = length(intersect(depolarize_all.(timeseries{Control(i)}).bottom10, repolarize_all.(timeseries{Control(i)}).bottom10))./length(repolarize_all.(timeseries{Control(i)}).bottom10);
%         Control_corrafter(i) = length(intersect(depolarize_all.(timeseries{Control(i)}).top10av, repolarize_all.(timeseries{Control(i)}).top10av))./length(repolarize_all.(timeseries{Control(i)}).top10av);
%    end
% 
%   CorrelationTable = array2table([Control_corrbefore, GKa_corrbefore, PKa_corrbefore; ...
%         Control_corrafter, GKa_corrafter, PKa_corrafter]', 'VariableNames',{'WRONG', 'Top 10'})
%     writetable(CorrelationTable, [savename, 'DepolRepolConsistency.csv'] )
% % 
% %KL End!
% Variables = {'COG_KL_wave', 'Degree_KL_wave'}
% for j = 1:length(Variables)
%     % Get regional consistency:
%     for i = 1:length(GKa)
%         GKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).(Variables{j})(1).trial));
%         GKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{GKa(i)}).(Variables{j})(2).trial));
%     end
%     for i = 1:length(PKa)
%         PKa_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).(Variables{j})(1).trial));
%         PKa_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{PKa(i)}).(Variables{j})(2).trial));
%     end
%     for i = 1:length(Control)
%         Control_corrbefore(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).(Variables{j})(1).trial));
%         Control_corrafter(i) = mean(nonzeros(wave_all_end.(timeseries{Control(i)}).(Variables{j})(2).trial))
%     end
% 
%     CorrelationTable = array2table([Control_corrbefore',  GKa_corrbefore', PKa_corrbefore'; ...
%     Control_corrafter', GKa_corrafter', PKa_corrafter']')
% 
%     writetable(CorrelationTable, [savename, Variables{j}, '.csv'] )
% 
%         CorrelationTable = array2table([Control_corrbefore; ...
%         Control_corrafter]', 'VariableNames',{'Pre-', 'Vehicle'})
%     writetable(CorrelationTable, [savename, Variables{j}, 'Wave_endControl.csv'] )
%         CorrelationTable = array2table([GKa_corrbefore; ...
%         GKa_corrafter]', 'VariableNames',{'Pre-', 'GKa'})
%     writetable(CorrelationTable, [savename, Variables{j}, 'Wave_endGKa.csv'] )
%         CorrelationTable = array2table([PKa_corrbefore; ...
%         PKa_corrafter]', 'VariableNames',{'Pre-', 'PKa'})
%     writetable(CorrelationTable, [savename, Variables{j}, 'Wave_endPKa.csv'] )
% end
% [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();
% 

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
