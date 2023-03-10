clc;
clear;
close all;

filepath_hic = input('入力Hi-Cファイルパスを入力してください：', "s");
heatmapMatrix = importdata(filepath_hic);
hic_g = graph(heatmapMatrix);
p = plot(hic_g,'MarkerSize',1 ,'EdgeAlpha', 0.0000001);
title('graph-cluster-test');
ucc = centrality(hic_g,'closeness');
p.NodeCData = ucc;
colormap jet;
colorbar
caxis([0.8*10^-4 1.2 *10^-4])
