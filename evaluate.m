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
atac_con_map_pre = zeros(atac_num^2, 2);
for i = 1:atac_num
    for r = 1:atac_num
        atac_con_map_pre((i-1)*atac_num + r, :) = [atac_con_pos(i,1) atac_con_pos(r,1)];
    end
end
atac_nocon_map_pre =zeros(atac_num^2, 2);
for i = 1:atac_num
    for r = 1:atac_num
        atac_nocon_map_pre((i-1)*atac_num + r, :) = [atac_mat(i,1) atac_mat(r,1)];
    end
end
Lia = ismember(atac_con_map_pre(:,1:2),sparse_pos(:,1:2),'rows');
atac_contact_map = atac_nocon_map_pre(Lia,:);
atac_contact_map_num = size(atac_contact_map,1);
sparse_pos_num = size(sparse_pos,1);
atac_contact_map_pos = atac_contact_map(:, 1:2);
[posMatrix,retVals]= ShNeigh(heatmapMatrix,1);
data = posMatrix(:,2:4);
datasize = size(data,1);
dataX_s = zeros(1,datasize * 100);
dataY_s = zeros(1,datasize * 100);
dataZ_s = zeros(1,datasize * 100);
dataX = posMatrix(:,2);
dataY = posMatrix(:,3);
dataZ = posMatrix(:,4);
for i = 1:datasize-1
    dataX_spliced = linspace(dataX(i),dataX(i+1));
    dataY_spliced = linspace(dataY(i),dataY(i+1));
    dataZ_spliced = linspace(dataZ(i),dataZ(i+1));
    dataX_s(1+100*(i-1):100*i) = dataX_spliced;
    dataY_s(1+100*(i-1):100*i) = dataY_spliced;
    dataZ_s(1+100*(i-1):100*i) = dataZ_spliced;
end
chrsize = 156340000;
split = datasize * 100;
bin = chrsize / split;
atac_size = size(atac_mat,1);
atac_idx = zeros(1, atac_size-1);
for i = 1:atac_size
    for s = 1:split
        if atac_mat(i,1) > bin * (s-1)&& atac_mat(i,1) < bin * s
            atac_idx(i) = s;
        end
    end
end
data_atac = zeros(atac_size, 3);
for i = 1:atac_size-1
    data_atac(i,1) = dataX_s(1,atac_idx(1,i));
    data_atac(i,2) = dataY_s(1,atac_idx(1,i));
    data_atac(i,3) = dataZ_s(1,atac_idx(1,i));
end
[IDX, isnoise] = DBSCAN(data_atac,0.025,25);
k = max(IDX);
Colors = hsv(k);
InClusterPair_all = [0 0];
for i = 1:k
    cluster_index = find(IDX==i);
    cluster_atac_pos = atac_mat(cluster_index,1);
    cluster_num = size(cluster_atac_pos,1);
    InClusterPair = zeros(cluster_num^2,2);
    for r = 1:cluster_num
        for s = 1:cluster_num
            if r <= s
                InClusterPair(r*(s-1)+s,:) = [cluster_atac_pos(r,1) cluster_atac_pos(s,1)];
            end
        end
    end
    InClusterPair_all = cat(1,InClusterPair_all,InClusterPair);
end
all_pair_num = size(InClusterPair_all,1);
InClusterPair_all(1,:) = [];

sia = ismember(InClusterPair_all,atac_contact_map,"rows");
InClusterPair_valid = InClusterPair_all(sia,:);
valid_pair_num = size(InClusterPair_valid,1);
per = valid_pair_num/all_pair_num * 100
