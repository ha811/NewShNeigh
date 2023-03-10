clc;
clear;
close all;

filepath_hic = input('入力Hi-Cファイルパスを入力してください：', "s");
filepath_atac = input('入力ATACファイルパスを入力してください：', "s");
hic_resolution = input('Hi-C解像度を入力してください：');
output_path = input('出力パスと名前を入力してください：', "s");
heatmapMatrix = importdata(filepath_hic);
atac_pos = importdata(filepath_atac);
hic_size=size(heatmapMatrix,1);
hic_num = size(heatmapMatrix, 1);
atac_score_mat = zeros(hic_num,1);
for i = 1:hic_num
    atac_score_mat(i,1) = nnz((atac_pos>(i-1)*hic_resolution) & (atac_pos<i*hic_resolution));
end
delete=find(atac_score_mat<hic_resolution/100);
delete_num = size(delete,1);
heatmapMatrix(delete,:) = zeros(delete_num,hic_num);
heatmapMatrix(:,delete) = zeros(hic_num,delete_num);
writematrix(heatmapMatrix,output_path);