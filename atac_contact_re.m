clc;
clear;
close all;

filepath_hic = input('入力Hi-Cファイルパスを入力してください：', "s");
filepath_atac = input('入力ATACファイルパスを入力してください：', "s");
output_path_matrix = input('出力行列ファイルパスと名前を入力してください：', "s");
output_path_heatmap = input('出力ヒートマップファイルパスと名前を入力してください：', "s");
hic_resolution = input('Hi-C解像度を入力してください：');
heatmapMatrix = importdata(filepath_hic);
[row,col,v] = find(heatmapMatrix);
sparse_pos = [row col v];
atac_mat = importdata(filepath_atac);
atac_pos = atac_mat(:, 1);
atac_num = size(atac_pos,1);
hic_num = size(heatmapMatrix, 1);
atac_con_pos = zeros(atac_num,2);
for i = 1:hic_num
    for r = 1:atac_num
        if atac_mat(r,1) > (i-1) * hic_resolution && atac_mat(r,1) < i * hic_resolution
            atac_con_pos(r, 1) = i;
            atac_con_pos(r, 2) = atac_mat(r,2);
        end
    end
end
atac_con_num = size(atac_con_pos, 1);
atac_con_map_pre = zeros(atac_con_num^2, 3);
for i = 1:atac_con_num
    for r = 1:atac_con_num
        atac_con_map_pre((i-1)*atac_con_num + r, :) = [atac_con_pos(i,1) atac_con_pos(r,1) atac_con_pos(i,2)*atac_con_pos(r,2)/(atac_con_pos(i,2) + atac_con_pos(r,2))];
    end
end
Lia = ismember(atac_con_map_pre(:,1:2),sparse_pos(:,1:2),'rows');
atac_contact_map = atac_con_map_pre(Lia,:);
atac_contact_map_num = size(atac_contact_map,1);
sparse_pos_num = size(sparse_pos,1);
atac_contact_map_pos = atac_contact_map(:, 1:2);

for i = 1:atac_contact_map_num
    atac_contact_map(i,3) = atac_contact_map(i, 3) * heatmapMatrix(atac_contact_map(i,1),atac_contact_map(i, 2));
end
contact_region = max(atac_contact_map(:,1));
data_heatmap = zeros(contact_region);
for i = 1:atac_contact_map_num
    data_heatmap(atac_contact_map(i,1),atac_contact_map(i,2)) = data_heatmap(atac_contact_map(i,1),atac_contact_map(i,2))+atac_contact_map(i,3);
end
fig = heatmap(data_heatmap,'GridVisible','off','ColorLimits',[0 15]);
saveas(fig,output_path_heatmap);



writematrix(atac_contact_map,output_path_matrix);

