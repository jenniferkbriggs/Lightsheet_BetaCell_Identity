%% This code is the script file to run first responder analysis
clear all
close all
clc


datapath = '/Volumes/Briggs_10TB/Merrin/2023Data/EJ130 2G10G/'

%indexing for later:
        ctpka = 1;
        ct = 1;

figure

for k = ["FirstPhase/"]%, "Bad/"]
files = dir(strcat(datapath, k))
for i = 1:length(files)
    %check if it is a folder:
    if isfolder(strcat(datapath, k, files(i).name))
        try %some folders are meta data with permissions so it will throw an error
        %if it is a folder, open it and load the data
           calcium = readmatrix(strcat(datapath, k, files(i).name, '/', files(i).name, '_Plot.csv'));
           Locations = readmatrix(strcat(datapath, k, files(i).name, '/', files(i).name, '_Pos.csv'));
           time = calcium(:,1);
           calcium(:,1) = [];
           %calcium = movmean(calcium,60,1);

           try
               load(strcat(datapath, k, files(i).name, '/', files(i).name, 'FirstResponderData.mat'), 'top10_2','top10_1','st','t_second_sort','t_first_sort','Locations')
           catch
               figure, plot(mean(calcium')), title('Mark the before beginning and after end of each first phase')
               st = ginput(2);
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
           clear t_first calcium2


           %Look at locations: 
           figure(1)
           j = 1
           top10 = ceil(0.1*size(calcium,2));
            nexttile
            scatter3(Locations(:,1), Locations(:,2), Locations(:,3), 75, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor',[0.7, 0.7, 0.7] , 'MarkerFaceAlpha', 0.5)
            hold on
            scatter3(Locations(t_first_sort(j,1:top10),1), Locations(t_first_sort(j,1:top10),2), Locations(t_first_sort(j,1:top10),3), 100, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor','blue')
            scatter3(Locations(t_first_sort(j,end-top10:end),1), Locations(t_first_sort(j,end-top10:end),2), Locations(t_first_sort(j,end-top10:end),3), 100, 'MarkerFaceColor', 'red', 'MarkerEdgeColor','red' )
    
            legend('Normal Cell','First Responder','Last Responder','Location', 'Northeast')
            set(gca, 'color','none')
            title(['File Name Number ' files(i).name])
        
            xlabel('X (\mum)')
            ylabel('Y (\mum)')
            zlabel('Z (\mum)')

            %calcium waveforms
            nexttile
             plot(time, calcium, 'color',[0.9,0.9,0.9])
             hold on, line1 = plot(time, calcium(:,t_first_sort(j,1:10)), 'linewidth',1, 'color', 'blue')
             hold on, line2 = plot(time, calcium(:,t_first_sort(j,end-10:end)), 'linewidth',1, 'color', 'red')
            xlim([(time(startindex(1))), time(startindex(2))])
            title(['File Name Number ' files(i).name])
            axg.data(j) = gca
            legend([line1(1), line2(1)], {'First Responder','Last Responder'})
            set(gca, 'color','none')
            xlabel('Time (s)')
            ylabel('Ca^{2+} Fluoresence')
            set(gca, 'box','off')


           %save data for later:
            save(strcat(datapath, k, files(i).name, '/', files(i).name, 'FirstResponderData.mat'),'top10_1','st','t_first_sort','Locations')
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