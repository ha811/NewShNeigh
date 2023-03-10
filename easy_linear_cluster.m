clc;
clear;
close all;

filepath_atac = input('入力ATACファイルパスを入力してください：', "s");
atac_mat = importdata(filepath_atac);
chrsize = 156340000;
bin = 40000;
bin_num = fix(chrsize / bin);
atac_score_mat = zeros(bin_num,1);
for i = 1:bin_num
    atac_score_mat(i,1) = nnz((atac_mat >= (i-1)*bin) &(atac_mat<=i*bin));
end
[row, ~, ~]=find(atac_score_mat>=400);
row_num = size(row,1);
col = zeros(row_num,1);
clus_pts = [row col];
[IDX,isnoise] = DBSCAN(clus_pts,10,8);
k = max(IDX);
cluster_score_mat = zeros(k,1);
cluster_pos_mat = zeros(bin_num,k);
Colors = hsv(7);
for i = 1:k
    cluster_score_mat(i,1) = sum(atac_score_mat(clus_pts(IDX==i),1));
    for r = 1:bin_num
        if nnz(clus_pts((IDX==i),1)==r) > 0
            cluster_pos_mat(r,i) = 1;
        end
    end
    clus_num = nnz(IDX==i);
end
[cluster_score_mat, I] = sort(cluster_score_mat,'descend');
Legends={};
for i=1:7
    cluster_info = find(cluster_pos_mat(:,I(i,1))) * 40000;
    y = zeros(size(cluster_info,1),1);
    scatter(cluster_info,y,5,Colors(i,:));
    Legends{end+1} = num2str(i);
    hold on
end
legend(Legends)
xlim([0 chrsize]);


