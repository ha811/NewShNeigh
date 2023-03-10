clc;
clear;
close all;

filepath_hic = input('入力Hi-Cファイルパスを入力してください：', "s");
output_path = input('出力ファイルパスを入力してください：', "s");
heatmapMatrix = importdata(filepath_hic);
[posMatrix,retVals shMat]= ShNeigh(heatmapMatrix,1);
datasize = size(posMatrix,1);
dataX = posMatrix(:,2);
dataY = posMatrix(:,3);
dataZ = posMatrix(:,4);
hic_num = size(heatmapMatrix,1);



hold on
plot3(dataX,dataY,dataZ, "blue");
saveas(gcf,output_path);

dist = pdist(posMatrix(:, 2:4));
dist_mat = squareform(dist);
writematrix(dist_mat, 'ShNeigh_reconmat_40000.txt');
writematrix(shMat, 'ShNeigh_distmat_40000.txt')
writematrix(posMatrix(:, 2:4), 'ShNeigh_XYZ_40000.txt')
