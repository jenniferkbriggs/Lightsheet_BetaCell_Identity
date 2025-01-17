function out = shownetwork_sim(adj, pos, Hubby, tit, savename)
savetime =  datestr(datetime('today'),'yyyymmdd');
AdjacencyGraph = sparse(adj);
%AdjacencyGraph = graph(A);
x = pos(:,1);
y = pos(:,2);
z = pos(:,3);

    

    maxz = max(z);
    minz = min(z);
    meanz = mean(z);
    
    cellplane = find( abs(z - meanz) < .05);
    
    xnew = x(cellplane);
    ynew = y(cellplane);
 
numcells = length(x);
for i=1:numcells
    N (i,1) = nnz(adj(:,i));  % N is matrix containing # of links for each cell (nnz returns number of nonzero elements)
end
histArray=zeros(size(N))'; % prealocate
% a forloop to count how many times you encounter in a particular value:
for n=1:length(N)
    histArray(1,N(n)+1)=histArray(1,N(n)+1)+1; % every time you meet the particular value, you add 1 into to corresponding bin
end
histArray = histArray(2:end); % removing first number corresponding to 0 links
histArrayPerc=histArray.*100/sum(histArray); % converting into % via dividing # of cell with specific # of links by total # of links 
m=find(histArray);    % returns all indexes of non-zero elements of the array
maxNonZeroLinks=max(m);   % number of links to consider when normalizing probabilty
k=1:1:maxNonZeroLinks;            % index of # of links (starting with 0 links)
kpercent=k.*100/(maxNonZeroLinks);  % convert # of links into % of limks

a = find(kpercent > 60)
a = a(1);
hubthreshold = k(a);
adjplane = adj(cellplane,cellplane);

connect = sum(adj);
[hubs] = find(connect>hubthreshold)%Number of cells linked to more than 50% of Islet
Nodec = repmat([.175 .54 .60],numcells,1);
Nodes = ones(numcells,1).*5;

%Hubby = intersect(hubs, cellplane);

for lll = 1:length(Hubby)
% Nodec(Hubby(lll),:)=	[0.6350, 0.0780, 0.1840];
Nodec(Hubby(lll),:)=	[0.6350, 0.0780, 0.1840];
Nodes(Hubby(lll)) = 8;
end
Nodec = Nodec(cellplane,:);
Nodes = Nodes(cellplane);

Hubsplane = find(Nodec(:,1) == 0.6350)
AdjacencyGraph = graph(adj(cellplane,cellplane));
Edgec = repmat([.175 .54 .60],length(table2array(AdjacencyGraph.Edges)),1);
Edges = ones(length(table2array(AdjacencyGraph.Edges)),1).*.75;
e = table2array(AdjacencyGraph.Edges);
for lll = 1:length(Hubsplane)
    indx = find(e(:,1) == Hubsplane(lll) | e(:,2) == Hubsplane(lll));
for tt = 1:length(indx)
    Edgec(indx(tt) ,:)=	[0.6350, 0.0780, 0.1840];
    Edges(indx(tt) ,:) = 1;
end
end


fig = figure
plot(AdjacencyGraph, 'Xdata',xnew,'YData',ynew, 'EdgeColor', Edgec, 'NodeColor',Nodec,'MarkerSize',Nodes, 'LineWidth',Edges )
title(tit, 'FontSize', 20)
set(gca,'ytick',[])
set(gca,'xtick',[])
C = repmat([.175 .54 .60], numcells,1);
for cc = 1:length(hubs)
C(hubs(cc),:) = [0.6350, 0.0780, 0.1840];
end
axes('Position',[.15 .12 .25 .25])
box on
scatter3(x,y,z,12,C,'Filled')
hold on
[X,Y]  = meshgrid(x,y);
Z = (ones(size(X))).*mean(z);
Cc = ones(numcells);
Cc(:,:,1) = .50;
Cc(:,:,2) = .50;
Cc(:,:,3) = .50;

mesh(X, Y, Z, Cc,'FaceAlpha', 0.5)

set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'ztick',[])

try
saveas(gcf, [pwd '\' savename '_' savetime '.fig'])
saveas(gcf, [pwd '\' savename '_' savetime '.png'])
end
end