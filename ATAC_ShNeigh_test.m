clc;
clear;
close all;

filepath_hic = input('入力Hi-Cファイルパスを入力してください：', "s");
output_path = input('出力ファイルパスを入力してください：', "s");
heatmapMatrix = importdata(filepath_hic);
[posMatrix,retVals]= ShNeigh(heatmapMatrix,1);
datasize = size(posMatrix,1);
dataX = posMatrix(:,2);
dataY = posMatrix(:,3);
dataZ = posMatrix(:,4);
hic_num = size(heatmapMatrix,1);
chrsize = hic_num * 40000;
bin = chrsize / datasize;


plot3(dataX,dataY,dataZ, "black");
saveas(gcf,output_path);

