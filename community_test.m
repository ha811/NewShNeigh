clc;
clear;
close all;
filepath_com = "/Users/hanayaren/Desktop/chr4_input/chr4_infomap_trimmed.txt";
com_mat = importdata(filepath_com);
filepath_atac = "/Users/hanayaren/Desktop/chr4_input/mm10_chr4_trimmed.sam";
atac_mat = importdata(filepath_atac);
hic_resolution = 40000;
hic_num = numel(com_mat);
chrsize=hic_num * hic_resolution;
atac_num = size(atac_mat,1);
atac_score_mat = zeros(hic_num,1);
for i = 1:hic_num
    atac_score_mat(i,1) = nnz((atac_mat>(i-1)*hic_resolution) & (atac_mat<i*hic_resolution));
end
k = max(max(com_mat));
com_score_mat = zeros(k,1);
for i = 1:k
    com_pos = find(com_mat==i);
    com_score_mat(i,1) = sum(atac_score_mat(com_pos,1));
end
[com_score_mat, I] = sort(com_score_mat,'descend');
com_score_mat = [I com_score_mat];
Colors=hsv(7);
Legends={};
for i = 1:7
    cluster_pos = find(com_mat==I(i,1)) * hic_resolution;
    y=zeros(1,size(cluster_pos,2));
    scatter(cluster_pos,y,5,Colors(i,:));
    Legends{end+1} = num2str(i);
    hold on
end
legend(Legends)
xlim([0 chrsize]);

first_clus = find(com_mat == I(1,1));
first_range = [min(first_clus), max(first_clus)] * hic_resolution;

second_clus = find(com_mat == I(2, 1));
second_range = [min(second_clus), max(second_clus)] * hic_resolution;

third_clus = find(com_mat == I(3, 1));
third_range = [min(third_clus), max(third_clus)] * hic_resolution;

ClusterRanges = [first_range; second_range; third_range];
writematrix(ClusterRanges, "/Users/hanayaren/Desktop/community_cluster/cluster_range.txt")


