clc;
clear;
close all;

filepath_hic = input('入力Hi-Cファイルパスを入力してください：', "s");
filepath_atac = input('入力ATACファイルパスを入力してください：', "s");
hic_resolution = input('Hi-C解像度を入力してください：');
heatmapMatrix = importdata(filepath_hic);
[row,col,v] = find(heatmapMatrix);
sparse_pos = [row col v];
atac_mat = importdata(filepath_atac);
atac_num = size(atac_mat,1);
hic_num = size(heatmapMatrix, 1);
atac_score_mat = zeros(hic_num,1);
for i = 1:hic_num
    for r = 1:atac_num
        if atac_mat(r,1) > (i-1) * hic_resolution && atac_mat(r,1) < i * hic_resolution
            atac_score_mat(i, 1) = atac_score_mat(i,1) + 1;
        end
    end
end
atac_con_map = zeros(hic_num);
for i = 1:hic_num
    for r = 1:hic_num
        if (atac_score_mat(i,1) + atac_score_mat(r,1)) ~= 0
            atac_con_map(i,r) = (atac_score_mat(i,1) * atac_score_mat(r,1))/(atac_score_mat(i,1) + atac_score_mat(r,1)) * heatmapMatrix(i,r);
        end
    end
end
[IDX, isnoise]=DBSCAN_graph(atac_con_map,25,500);
IDX_num = max(IDX);
atac_I_mat = zeros(IDX_num,1);
for i = 1:IDX_num
    cluster_pos = find(IDX==i);
    atac_cluster = atac_score_mat(cluster_pos);
    cluster_score = sum(atac_cluster);
    atac_I_mat(i,1) = cluster_score;
end
[atac_I_mat, I] = sort(atac_I_mat,'descend');
atac_I_mat = [I atac_I_mat]
top_cluster_pos = find(IDX==atac_I_mat(1,1));

