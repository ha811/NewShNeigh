clc;
clear;
close all;

filepath_xyz = input('入力xyzファイルパスを入力してください：', "s");
output_path = input('出力ファイルパスを入力してください：', "s");
posMatrix = importdata(filepath_xyz);

dataX = posMatrix(:,1);
dataY = posMatrix(:,2);
dataZ = posMatrix(:,3);

scatter3(dataX,dataY,dataZ, "blue");
saveas(gcf,output_path);


