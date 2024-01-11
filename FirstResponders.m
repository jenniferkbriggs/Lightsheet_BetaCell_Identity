%% This code is the script file to run first responder analysis
clear all
close all
clc


datapath = '/Volumes/Briggs_10TB/Merrin/2023Data/EJ130 2G10G/'

%indexing for later:
        ctpka = 1;
        ct = 1;



for j = ["Good/", "Medium/"]%, "Bad/"]
files = dir(strcat(datapath, j))
for i = 1:length(files)
    %check if it is a folder:
    if isfolder(strcat(datapath, j, files(i).name))
        try %some folders are meta data with permissions so it will throw an error
        %if it is a folder, open it and load the data
           calcium = readmatrix(strcat(datapath, j, files(i).name, '/', files(i).name, '_Plot.csv'));
           locations = readmatrix(strcat(datapath, j, files(i).name, '/', files(i).name, '_Pos.csv'));
           time = calcium(:,1);
           calcium(:,1) = [];
           %these data are so spiky that maybe we try to take the slow
           %frequencies....
           % fs = mean(diff(time));
           % y = lowpass(mean(calcium'),0.0001, fs);
           %calcium = rescale(calcium, 0,1);
           calcium = movmean(calcium,60,1);

           try
               load(strcat(datapath, j, files(i).name, '/', files(i).name, 'FirstResponderData.mat'), 'top10_2','top10_1','st','t_second_sort','t_first_sort','locations')
           catch
               figure, plot(mean(calcium')), title('Mark the before beginning and after end of each first phase')
               st = ginput(4);
           end

           % -- first depolarization -- $
           %rescale
           startindex = round(st(1:2,1));
           for ce = 1:size(calcium,2)
                calcium2(:,ce) = rescale(calcium(startindex(1):startindex(2),ce), 0, 1);
           end

           %interpolate: 
           calcium2 = interp1(calcium2, [1:0.1:size(calcium2,1)]);

           %find time past 50% for each cell
           for ce = 1:size(calcium,2)
               tcell = find(calcium2(:,ce)>0.5);
               t_first(ce) = tcell(1);
           end

           %sort
           [~, t_first_sort] = sort(t_first);
           top10_1 = t_first_sort(1:ceil(0.1*size(calcium,2)));
           clear calcium2 t_first

           % -- second depolarization -- $
           %rescale
           startindex = round(st(3:4,1));
           for ce = 1:size(calcium,2)
                calcium2(:,ce) = rescale(calcium(startindex(1):startindex(2),ce), 0, 1);
           end

           %find time past 50% for each cell
           for ce = 1:size(calcium,2)
               tcell = find(calcium2(:,ce)>0.5);
               t_second(ce) = tcell(1);
           end

           %sort
           [tt, t_second_sort] = sort(t_second);
           top10_2 = t_second_sort(1:ceil(0.1*size(calcium,2)));
           clear calcium2 t_second 

           %save data for later:
            save(strcat(datapath, j, files(i).name, '/', files(i).name, 'FirstResponderData.mat'), 'top10_2','top10_1','st','t_second_sort','t_first_sort','locations')
        %calculate percent maintained and split between pka and control:

        if contains(files(i).name, "PKa")
            percent_maintained(ctpka,2) = length(intersect(top10_1,top10_2))/length(top10_1);
            ctpka = ctpka + 1;
        else
            percent_maintained(ct,1) = length(intersect(top10_1,top10_2))/length(top10_1);
            ct = ct + 1;
        end
        end
    end
end

end