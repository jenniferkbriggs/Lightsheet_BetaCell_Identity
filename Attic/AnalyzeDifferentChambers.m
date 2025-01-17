fields = fieldnames(wave_all)
for i = 1:length(fields)
    j = fields{i}
    for k = 1:size(wave_all.(j).top10, 1)
        for l = 1:size(wave_all.(j).top10, 1)
            if k <= l
                wave_all.(j).intravarall(k,l) = NaN;
            else
                wave_all.(j).intravarall(k,l) = length(intersect(wave_all.(j).top10(k,:), wave_all.(j).top10(l,:)))/length(wave_all.(j).top10(l,:));
            end
        end
    end
end


var = COG_KL_net.trial
ylab = 'Consistency'




figure, nexttile

COG_KLnet = [mean(mean(network_all.EJ13110GF121noflow_Statistics.intravarall, 'omitnan'), 'omitnan'), ...
    mean(mean(network_all.EJ13110GF15115p_Statistics.intravarall, 'omitnan'), 'omitnan'),...
    mean(mean(network_all.EJ12910GF06.intravarall, 'omitnan'), 'omitnan'),...
    mean(mean(network_all.EJ12910GF09.intravarall, 'omitnan'), 'omitnan'),...
      mean(mean(network_all.EJ12910GF11.intravarall, 'omitnan'), 'omitnan')]
X = categorical({'Control No Flow', 'Control Flow','Pre-Chamber Change 1','Pre-Chamber Change 2','Pre-Chamber Change 3'})
bar(X, COG_KLnet)
ylabel('Average percent of high degree cells maintained')

set(gca, 'FontSize',15)
set(gcf, 'color','white')
set(gca, 'box','off')

nexttile


COG_KLnet = [mean(mean(wave_all.EJ13110GF121noflow_Statistics.intravarall, 'omitnan'), 'omitnan'), ...
    mean(mean(wave_all.EJ13110GF15115p_Statistics.intravarall, 'omitnan'), 'omitnan'),...
    mean(mean(wave_all.EJ12910GF06.intravarall, 'omitnan'), 'omitnan'),...
    mean(mean(wave_all.EJ12910GF09.intravarall, 'omitnan'), 'omitnan'),...
      mean(mean(wave_all.EJ12910GF11.intravarall, 'omitnan'), 'omitnan')]
X = categorical({'Control No Flow', 'Control Flow','Pre-Chamber Change 1','Pre-Chamber Change 2','Pre-Chamber Change 3'})
bar(X, COG_KLnet)


set(gca, 'FontSize',15)
set(gcf, 'color','white')
set(gca, 'box','off')

ylabel('Average percent of high phase cells maintained')


%linear regression

COG_KLnet = [[1:length(network_all.EJ13110GF121noflow_Statistics.correlation)]'\network_all.EJ13110GF121noflow_Statistics.correlation', ...
[1:length(network_all.EJ13110GF15115p_Statistics.correlation)]'\network_all.EJ13110GF15115p_Statistics.correlation', ...
[1:length(network_all.EJ12910GF06.correlation)]'\network_all.EJ12910GF06.correlation', ...
[1:length(network_all.EJ12910GF09.correlation)]'\network_all.EJ12910GF09.correlation', ...
[1:length(network_all.EJ12910GF11.correlation)]'\network_all.EJ12910GF11.correlation']
