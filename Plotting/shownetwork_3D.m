function out = shownetwork_3D(adj, pos, tit, savename,hubby)
%function to view 3D network (e.g. Figure 4)
%Jennifer Briggs 2020
%Inputs:
%adj -- adjacency matrix (square)
%pos -- 3xCell Number XYZ position matrix
%both -- graph types
%Hubby -- indices corresponding to hubs
%tit -- title of image (string)
%savename -- name the image will be saved as
savetime =  datestr(datetime('today'),'yyyymmdd');
if 0

importdata([pwd '\' savename '_' savetime '.m'])
names = fieldnames(ans)
indx = strfind(names, 'ans')
for i = 1:length(names)
    if indx{i} ~= 1
    eval([names{i} '=ans.' names{i}]);
    end
end
end

AdjacencyGraph = graph(adj, 'upper');
x = pos(:,1);
y = pos(:,2);
z = pos(:,3);
numcell = length(adj)
repcellh = randi([1 length(hubby)],1,1);
repcellh = hubby(repcellh)

%zero degree cells
s = sum(adj, 'omitnan');
repcelll = find(s == 0);

%1 degree cells
repcell_one = find(sum(adj, 'omitnan') == 1);
repcell_two = find(sum(adj, 'omitnan') == 2);



cellmat = zeros(numcell,numcell);
cellmat([repcellh],:) = 1;
cellmat(:,[repcellh]) = 1;
adj2 = adj.*cellmat;

g1 = graph(adj2,'upper');

% 
% %graph with one or two edges
% cellmat = zeros(numcell,numcell);
% cellmat([repcellh],:) = 1;
% cellmat(:,[repcellh]) = 1;
% adj2 = adj.*cellmat;
% g2 = graph(ones(length([repcell_one'; repcell_two']), length([repcell_one'; repcell_two'])), 'upper');

%color nodes blue and red
g1c = repmat([0 0 0],numcell,1); %BLUE
g1c(repcellh,1) = 0.01%[0.6350];
g1c(repcellh,2) = 0.01%0.0780;
g1c(repcellh,3)= 1%0.1840;

%find connections
[allconnections, ~] = find(adj(:, repcellh) == 1);
allconnections = unique(allconnections);
g1c(allconnections,1) = 0.01%[0.6350];
g1c(allconnections,2) = 0.01%0.0780;
g1c(allconnections,3)= 1%0.1840;

%RED - low degree cells
repcelll = unique(repcelll);
g1c([(repcelll)] ,1) =  1;%0.3650;
g1c([repcelll],2) = .01;%0.9220;
g1c([repcelll],3)= .01;%0.8160;

%color all other nodes black
g1c(g1c == 0) = .5;

Edgec1 = repmat([0.01 0.01 1],size(table2array(g1.Edges),1),1);


Msize = ones(numcell, 1).*7;
Msize(repcellh) = 15;



figure, plot(g1,'Xdata',x,'YData',y,'ZData', z,'EdgeColor', Edgec1, ...
    'NodeColor',g1c, 'LineWidth', 2, 'MarkerSize', Msize)
%plot just low degree (not zero, that was plotted with g1)
% hold on, plot(g2,'Xdata',x([repcell_one'; repcell_two']),'YData',...
%     y([repcell_one'; repcell_two']),'ZData', z([repcell_one'; repcell_two']),...
%     'EdgeColor', [1, 0.1, 0.1], ...
%     'NodeColor',[1, 0.1, 0.1], 'LineWidth', 2, 'MarkerSize', 15)
% 
% figure, plot(g, 'Xdata',x,'YData',y,'ZData', z,'EdgeColor', Edgec, ...
%     'NodeColor',g1c, 'LineWidth', 2, 'MarkerSize', Msize, 'NodeLabel', Nnames)
% hold on
% plot(g2, 'Xdata',x,'YData',y,'ZData', z)
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'ztick',[])
set(gcf, 'color','white')
end