clc;
clear;
close all;

filepath_hic = input('入力Hi-Cファイルパスを入力してください：', "s");
filepath_atac = input('入力ATACファイルパスを入力してください：', "s");
output_path_scat = input('出力散布図パスと名前を入力してください：', "s");
hic_resolution = input('Hi-C解像度を入力してください：');
heatmapMatrix = importdata(filepath_hic);
[row,col,v] = find(heatmapMatrix);
sparse_pos = [row col v];
atac_pos = importdata(filepath_atac);
hic_num = size(heatmapMatrix, 1);
atac_score_mat = zeros(hic_num,1);
for i = 1:hic_num
    atac_score_mat(i,1) = nnz((atac_pos>(i-1)*hic_resolution) & (atac_pos<i*hic_resolution));
end
atac_con_map = zeros(hic_num);
for i = 1:hic_num
    for r = 1:hic_num
        if (atac_score_mat(i,1) + atac_score_mat(r,1) ~= 0) && (i < r) && ((r-i)>100)
            atac_con_map(i,r) = (atac_score_mat(i,1) * atac_score_mat(r,1))/(atac_score_mat(i,1) + atac_score_mat(r,1)) * heatmapMatrix(i,r);
        end
    end
end
[row,col,v] = find(atac_con_map);
atac_scat=cat(2,row,col,v);
atac_scat_num = size(atac_scat,1);
for i = 1:atac_scat_num
    if atac_scat(i,3) < 400
        atac_scat(i,3) = 0;
    end
end
[row, col] = find(~atac_scat(:,3));
atac_scat(row,:) = [];
[IDX,isnoise] = DBSCAN(atac_scat(:,1:2),50,20);
k = max(IDX);
Colors = hsv(k);
cluster_score_mat = zeros(k,1);
cluster_pos_mat = zeros(hic_num,k);
for i = 1:k
    cluster_atac=atac_scat(IDX==i, :);
    for r = 1:hic_num
        cluster_pos_mat(r,i) = nnz(cluster_atac(:,1:2)==r);
    end
    cluster_score_mat(i,1) = sum(atac_score_mat(cluster_pos_mat(:,i)>0));
    color = Colors(i,:);    
    scatter(cluster_atac(:, 1), cluster_atac(:, 2), 3, color);
    hold on
end
[cluster_score_mat, I] = sort(cluster_score_mat,'descend');
cluster_info = find(cluster_pos_mat(:,I(1,1)))
cluster_score_mat



saveas(gcf,output_path_scat);


