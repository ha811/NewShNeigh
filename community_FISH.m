clc;
clear;
close all;
chr = 3;
method = 'info';
filepath_com_form = "/Users/hanayaren/Desktop/2022_hirota_ATAC/community/chr%d_40000_%s_trimmed.txt";
filepath_com = sprintf(filepath_com_form, chr, method);
com_mat = importdata(filepath_com);
filepath_atac_form = "/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC/NoTre_chr%d_trimmed.sam";
filepath_atac = sprintf(filepath_atac_form, chr);
atac_mat = importdata(filepath_atac);
hic_resolution = 40000;
hic_num = numel(com_mat);
chrsize=hic_num * hic_resolution;
atac_num = size(atac_mat,1);
atac_score_mat = zeros(hic_num,1);
for i = 1:hic_num
    atac_score_mat(i,1) = nnz((atac_mat>(i-1)*hic_resolution) & (atac_mat<i*hic_resolution));
end
atac_pos = find(atac_score_mat>1000);
k = max(max(com_mat));
com_score_mat = zeros(k,1);
for i = 1:k
    com_pos = find(com_mat==i);
    com_score_mat(i,1) = sum(atac_score_mat(com_pos,1));
end
[com_score_mat, I] = sort(com_score_mat,'descend');
com_score_mat = [I com_score_mat];
Colors=hsv(3);
Legends = {};
f = figure();
f.Position(3:4) = [560 210];
figure(f)
for i = 1:3
    
    cluster_pos = find(com_mat==I(i,1));
    atac_cluster_pos = intersect(cluster_pos, atac_pos) * hic_resolution;
    y=zeros(1,size(atac_cluster_pos,1));
    scatter(atac_cluster_pos,y,5,Colors(i,:));
    Legends{end+1} = num2str(i);
    hold on
end
yline(0);

lgd = legend(Legends);
title(lgd, 'ATAC cluster ranked')
xlim([0 chrsize]);
ylim([-0.2, 1]);
xlabel('ゲノム座標[bp]')
ax = gca;
ax.YTickLabel = cell(size(ax.YTickLabel));
ax.YTick = [];
title(sprintf('ゲノム座標上でのクラスター形成予測位置(chr%d)', chr))


first_clus = find(com_mat == I(1,1)) * hic_resolution;
first_range = ['1_min', mat2str(mink(first_clus,1)), '1_max',mat2str(maxk(first_clus,1))];

second_clus = find(com_mat == I(2, 1)) * hic_resolution;
second_range = ['2_min', mat2str(mink(second_clus, 1)), '2_max',mat2str(maxk(second_clus, 1))];

third_clus = find(com_mat == I(3, 1)) * hic_resolution;
third_range = ['3_min', mat2str(mink(third_clus, 1)), '3_max', mat2str(maxk(third_clus, 1))];

com_path_form = "/Users/hanayaren/Desktop/2022_hirota_ATAC/clusrer_range/chr%d_cluster_range_%s.txt";
com_path = sprintf(com_path_form, chr, method);


writelines(first_range, com_path)
writelines(second_range, com_path, WriteMode="append")
writelines(third_range, com_path, WriteMode="append")
fig_path_form = '/Users/hanayaren/Desktop/2022_hirota_ATAC/clusrer_range/chr%d_%s_com.png';
fig_path = sprintf(fig_path_form, chr, method);
saveas(gcf, fig_path)


