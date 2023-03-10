clc;
clear;
close all;

filepath = input('入力ファイルパスを入力してください：', "s");
output_path = input('出力ファイルパスを入力してください：', "s");
heatmapMatrix = importdata(filepath);
[posMatrix,retVals]= ShNeigh(heatmapMatrix,1);
dataX = posMatrix(:,2);
dataY = posMatrix(:,3);
dataZ = posMatrix(:,4);
plot3(dataX,dataY,dataZ);
saveas(gcf,output_path);



