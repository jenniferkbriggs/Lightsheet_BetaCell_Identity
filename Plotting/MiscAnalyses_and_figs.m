
% Normalize the correlation data to range [0, 1]
normalized_data = (correlation_data - min(correlation_data)) / (max(correlation_data) - min(correlation_data));

% Define a colormap (e.g., 'jet')
cmap = (jet(length(normalized_data)))
% Map the normalized data to colormap indices
colormap_indices = round(normalized_data * (colormap_length - 1)) + 1;

% Get the RGB values from the colormap
rgb_values = colormap(rescale(correlation_data, 0, 1));

% Display the RGB values
disp(rgb_values);


figure;
scatter3(Locations(:,1), Locations(:,2),Locations(:,3), 50, cmap, 'filled');
title('Correlation Data Visualization');
xlabel('X-axis');
ylabel('Y-axis');
colorbar;  % Display colorbar to show the colormap
colormap(colormap_name);  % Use the same colormap for colorbar


%% Make table of number of cells: 
%extract data
timeseries = fieldnames(wave_all)
GKa = find(contains(timeseries, 'GKa'));
PKa = find(contains(timeseries, 'PKa'));
Control = setxor([1:length(timeseries)], [GKa; PKa]);
Control(5) = [];
GKa(1) = [];
% PKa(3) = [];



T = table('Size',[length(timeseries), 2], 'VariableTypes', ["string","double"], 'VariableNames',{'Name', 'Number of Cells'})

%Control
ct = 1
for i = 1:length(Control)
   name = timeseries{Control(i)}
   T(ct,1) = cellstr(name(1:end-11));
   T(ct,2) = num2cell(length(wave_all.(timeseries{Control(i)}).phase_pka(1,:)));
   ct = ct+1
end

for i = 1:length(GKa)
   name = timeseries{GKa(i)}
   T(ct,1) = cellstr(name(1:end-11));
   T(ct,2) = num2cell(length(wave_all.(timeseries{GKa(i)}).phase_pka(1,:)));
   ct = ct+1
end

for i = 1:length(PKa)
   name = timeseries{PKa(i)}
   T(ct,1) = cellstr(name(1:end-11));
   T(ct,2) = num2cell(length(wave_all.(timeseries{PKa(i)}).phase_pka(1,:)));
   ct = ct+1
end

writetable(T, 'NumberOfCells.csv')