%extract data
timeseries = fieldnames(network_all)
GKa = find(contains(timeseries, 'GKa'));
PKa = find(contains(timeseries, 'PKa'));
Control = setxor([1:length(timeseries)], [GKa; PKa]);
Control(5) = [];
GKa(1) = [];
% PKa(3) = [];

savename  = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange/2Dvs3D/'

[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkcells();

    for i = 1:length(GKa)
        GKa_corrbefore(i) = ({timeseries{GKa(i)}});
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = ({timeseries{PKa(i)}});
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = ({timeseries{Control(i)}});
    end
    IsletIndices = [Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore'];
    
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

%Location ------
%find average distance from center:

%threeD
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).threeD.dist_from_center));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).threeD.dist_from_center));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).threeD.dist_from_center));
    end
    average = [Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore'];
    
    
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).threeD.distcenter_top10));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threeD.distcenter_top10));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).threeD.distcenter_top10));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threeD.distcenter_top10));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).threeD.distcenter_top10));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threeD.distcenter_top10));
    end
    
        Network_loc = [Control_corrbefore, GKa_corrbefore, PKa_corrbefore]';
    Wave_loc = [Control_corrafter,  GKa_corrafter, PKa_corrafter]';

    CorrelationTable = array2table([average, Network_loc], 'VariableNames', {'Average Cell','High Degree'})
        writetable(CorrelationTable, [savename, 'ThreeDNetworkvsAverageLoc.csv'] )

    CorrelationTable = array2table([average, Wave_loc], 'VariableNames', {'Average Cell','High Phase'})
        writetable(CorrelationTable, [savename, 'ThreeDWavevsAverageLoc.csv'] )


%quarter
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).threequarter2D.dist_from_center));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).threequarter2D.dist_from_center));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).threequarter2D.dist_from_center));
    end
    average = [Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore'];
    
    
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).threequarter2D.distcenter_top10));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threequarter2D.distcenter_top10));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).threequarter2D.distcenter_top10));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threequarter2D.distcenter_top10));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).threequarter2D.distcenter_top10));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threequarter2D.distcenter_top10));
    end
    
    Network_loc = [Control_corrbefore, GKa_corrbefore, PKa_corrbefore]';
    Wave_loc = [Control_corrafter,  GKa_corrafter, PKa_corrafter]';

    CorrelationTable = array2table([average, Network_loc], 'VariableNames', {'Average Cell','High Degree'})
        writetable(CorrelationTable, [savename, 'threequarter2DNetworkvsAverageLoc.csv'] )

    CorrelationTable = array2table([average, Wave_loc], 'VariableNames', {'Average Cell','High Phase'})
        writetable(CorrelationTable, [savename, 'threequarter2DWavevsAverageLoc.csv'] )


%half2D
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).half2D.dist_from_center));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).half2D.dist_from_center));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).half2D.dist_from_center));
    end
    average = [Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore'];
    
    
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).half2D.distcenter_top10));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).half2D.distcenter_top10));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).half2D.distcenter_top10));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).half2D.distcenter_top10));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(network_all.(timeseries{Control(i)}).half2D.distcenter_top10));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).half2D.distcenter_top10));
    end
    
        Network_loc = [Control_corrbefore, GKa_corrbefore, PKa_corrbefore]';
    Wave_loc = [Control_corrafter,  GKa_corrafter, PKa_corrafter]';

    CorrelationTable = array2table([average, Network_loc], 'VariableNames', {'Average Cell','High Degree'})
        writetable(CorrelationTable, [savename, 'half2DNetworkvsAverageLoc.csv'] )

    CorrelationTable = array2table([average, Wave_loc], 'VariableNames', {'Average Cell','High Phase'})
        writetable(CorrelationTable, [savename, 'half2DWavevsAverageLoc.csv'] )

[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();


%% ---- KL Cell vs Region ---%% 
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

% 3D ---
    for j = 1:length(Variables)
        % Get regional consistency:
        for i = 1:length(GKa)
            GKa_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threeD.(Variables{j})(1).trial));
            GKa_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).threeD.(Variables2{j})(1).trial));
        end
        for i = 1:length(PKa)
            PKa_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threeD.(Variables{j})(1).trial));
            PKa_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).threeD.(Variables2{j})(1).trial));
        end
        for i = 1:length(Control)
            Control_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threeD.(Variables{j})(1).trial));
            Control_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{Control(i)}).threeD.(Variables2{j})(1).trial));
        end
        
         
    end
       CorrelationTable = array2table(fliplr([Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore']), 'VariableNames',{'Cellular', 'Regional'})
        writetable(CorrelationTable, [savename, 'threeDCellvsRegionPhase.csv'] )
    
       CorrelationTable = array2table(fliplr([Control_corrafter'; GKa_corrafter'; PKa_corrafter']), 'VariableNames',{'Cellular', 'Regional'})
        writetable(CorrelationTable, [savename, 'threeDCellvsRegionNetwork.csv'] )   
        

% 2D half ---
    for j = 1:length(Variables)
        % Get regional consistency:
        for i = 1:length(GKa)
            GKa_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).half2D.(Variables{j})(1).trial));
            GKa_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).half2D.(Variables2{j})(1).trial));
        end
        for i = 1:length(PKa)
            PKa_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).half2D.(Variables{j})(1).trial));
            PKa_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).half2D.(Variables2{j})(1).trial));
        end
        for i = 1:length(Control)
            Control_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).half2D.(Variables{j})(1).trial));
            Control_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{Control(i)}).half2D.(Variables2{j})(1).trial));
        end
        
         
    end
       CorrelationTable = array2table(fliplr([Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore']), 'VariableNames',{'Cellular', 'Regional'})
        writetable(CorrelationTable, [savename, 'half2DCellvsRegionPhase.csv'] )
    
       CorrelationTable = array2table(fliplr([Control_corrafter'; GKa_corrafter'; PKa_corrafter']), 'VariableNames',{'Cellular', 'Regional'})
        writetable(CorrelationTable, [savename, 'half2DCellvsRegionNetwork.csv'] )   
        

% three quarter ---
    for j = 1:length(Variables)
        % Get regional consistency:
        for i = 1:length(GKa)
            GKa_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threequarter2D.(Variables{j})(1).trial));
            GKa_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{GKa(i)}).threequarter2D.(Variables2{j})(1).trial));
        end
        for i = 1:length(PKa)
            PKa_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threequarter2D.(Variables{j})(1).trial));
            PKa_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{PKa(i)}).threequarter2D.(Variables2{j})(1).trial));
        end
        for i = 1:length(Control)
            Control_corrbefore(j,i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threequarter2D.(Variables{j})(1).trial));
            Control_corrafter(j,i) = mean(nonzeros(network_all.(timeseries{Control(i)}).threequarter2D.(Variables2{j})(1).trial));
        end
        
         
    end
       CorrelationTable = array2table(fliplr([Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore']), 'VariableNames',{'Cellular', 'Regional'})
        writetable(CorrelationTable, [savename, 'threequarter2DCellvsRegionPhase.csv'] )
    
       CorrelationTable = array2table(fliplr([Control_corrafter'; GKa_corrafter'; PKa_corrafter']), 'VariableNames',{'Cellular', 'Regional'})
        writetable(CorrelationTable, [savename, 'threequarter2DCellvsRegionNetwork.csv'] )   
        

%% -- Rotation -- %%
[Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkarrays();

    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threeD.rot));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threequarter2D.rot));
        GKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).half2D.rot));
    end

    for i = 1:length(GKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threeD.rot));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threequarter2D.rot));
        PKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).half2D.rot));
    end

    for i = 1:length(GKa)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threeD.rot));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threequarter2D.rot));
        Control_corrboth(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).half2D.rot));
    end

       CorrelationTable = array2table(fliplr([[Control_corrbefore'; GKa_corrbefore'; PKa_corrbefore'], ...
           [Control_corrafter'; GKa_corrafter'; PKa_corrafter'], ...
           [Control_corrboth'; GKa_corrboth'; PKa_corrboth']]), 'VariableNames',{'3D', '1/4 Z-stack', '1/2 Z-stack'})
        writetable(CorrelationTable, [savename, 'Waveaxis.csv'] )
    


%% 3D wave axis after proper normalization ---
savename  = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Analysis_After_ChamberChange/'

% Intravaraibility Network
VarNames = {'Pre- intravariability', 'Vehicle intravariability', 'Intervariability'}
    % Get regional consistency:
    for i = 1:length(GKa)
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threeD.rot));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threeD.rot_pka));
        GKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threeD.rot_intra));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threeD.rot));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threeD.rot_pka));
        PKa_corrboth(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threeD.rot_intra));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threeD.rot));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threeD.rot_pka));
        Control_corrboth(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threeD.rot_intra));
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
        GKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threeD.rot));
        GKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{GKa(i)}).threeD.rot_pka));
    end
    for i = 1:length(PKa)
        PKa_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threeD.rot));
        PKa_corrafter(i) = mean(nonzeros(wave_all.(timeseries{PKa(i)}).threeD.rot_pka));
    end
    for i = 1:length(Control)
        Control_corrbefore(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threeD.rot));
        Control_corrafter(i) = mean(nonzeros(wave_all.(timeseries{Control(i)}).threeD.rot_pka));
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



   %%% ------- Making empty array ----- %%%% 
    
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

    
function [Control_corrafter Control_corrbefore GKa_corrbefore GKa_corrafter PKa_corrbefore PKa_corrafter Control_corrboth PKa_corrboth GKa_corrboth]= mkcells()
PKa_corrbefore = cell(1,14);
PKa_corrafter = cell(1,14);
PKa_corrboth = cell(1,14);
GKa_corrbefore = cell(1,14);
GKa_corrafter = cell(1,14);
GKa_corrboth = cell(1,14);
Control_corrbefore = cell(1,14);
Control_corrafter = cell(1,14);
Control_corrboth = cell(1,14);
Control_corrboth = cell(1,14);
PKa_corrboth = cell(1,14);
GKa_corrboth = cell(1,14);
end

