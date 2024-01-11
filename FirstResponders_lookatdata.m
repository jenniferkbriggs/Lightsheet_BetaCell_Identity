%% This code is the script file to run first responder analysis
datapath = '/Volumes/Briggs_10TB/Merrin/2023Data/EJ130 2G10G/'
figure
for j = ["Good/", "Medium/", "Bad/"]
files = dir(strcat(datapath, j))
nexttile, text
       set(gca, 'box','off')
        set(gca, 'color','none')
for i = 1:length(files)
    %check if it is a folder:
    if isfolder(strcat(datapath, j, files(i).name))
        try %some folders are meta data with permissions so it will throw an error
        %if it is a folder, open it and load the data
       calcium = readmatrix(strcat(datapath, j, files(i).name, '/', files(i).name, '_Plot.csv'));
       locations = readmatrix(strcat(datapath, j, files(i).name, '/', files(i).name, '_Pos.csv'));
       time = calcium(:,1);
       calcium(:,1) = [];
       calcium = rescale(calcium, 0,1);

       
       nexttile, plot(time, mean(calcium')), xlabel('Time'), ylabel('Average Calcium'), 
       set(gca, 'box','off')
       set(gca, 'color','none')
       title(files(i).name)


        end
    end
end

end