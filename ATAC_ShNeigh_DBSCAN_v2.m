clc;
clear;
close all;

filepath_hic = input('入力Hi-Cファイルパスを入力してください：', "s");
filepath_atac = input('入力ATACファイルパスを入力してください：', "s");
output_path = input('出力ファイルパスを入力してください：', "s");
heatmapMatrix = importdata(filepath_hic);
atac_pos = importdata(filepath_atac);
[posMatrix,retVals]= ShNeigh(heatmapMatrix,1);
datasize = size(posMatrix,1);
dataX = posMatrix(:,2);
dataY = posMatrix(:,3);
dataZ = posMatrix(:,4);
hic_num = size(heatmapMatrix,1);
chrsize = hic_num * 40000;
bin = chrsize / datasize;
atac_num = size(atac_pos,1);
atac_score_mat = zeros(datasize,1);
for i = 1:datasize
    atac_score_mat(i,1) = nnz((bin*(i-1)<=atac_pos) & (atac_pos<=i*bin));
end
atac_idx = find(atac_score_mat>=300);
atac_size = nnz(atac_score_mat >=300);
data_atac = zeros(atac_size, 3);
for i = 1:atac_size-1
    data_atac(i,1) = dataX(atac_idx(i,1),1);
    data_atac(i,2) = dataY(atac_idx(i,1),1);
    data_atac(i,3) = dataZ(atac_idx(i,1),1);
end
[IDX, isnoise] = DBSCAN(data_atac,0.016,25);
k = max(IDX);
Colors = hsv(k);
cluster_score_mat=zeros(k,1);
cluster_pos_mat = zeros(atac_size,k);
for i = 1:k
    cluster_atac=data_atac(IDX==i, :);
    for r = 1:atac_size
        if IDX(r,1) == i
            cluster_pos_mat(r,i) = 1;
        end
    end
    cluster_score_mat(i,1) = sum(atac_score_mat(IDX==i,1));
    color = Colors(i,:);    
    scatter3(cluster_atac(:, 1), cluster_atac(:, 2), cluster_atac(:, 3), 10, color);
    hold on
end
hold on
plot3(dataX,dataY,dataZ, "black");
saveas(gcf,output_path);
delete(gcf)
[cluster_score_mat, I] = sort(cluster_score_mat,'descend');
Legends = {};
for i = 1:7
    cluster_info = find(cluster_pos_mat(:,I(i,1)));
    cluster_real_pos = atac_idx(cluster_info) * 40000;
    y = zeros(size(cluster_info,1),1);
    scatter(cluster_real_pos,y,5,Colors(i,:));
    Legends{end+1} = num2str(i);
    hold on
end
legend(Legends)
xlim([0 chrsize]);

