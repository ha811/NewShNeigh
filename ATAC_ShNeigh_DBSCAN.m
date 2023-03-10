clc;
clear;
close all;

filepath_hic = input('入力Hi-Cファイルパスを入力してください：', "s");
filepath_atac = input('入力ATACファイルパスを入力してください：', "s");
output_path = input('出力ファイルパスを入力してください：', "s");
heatmapMatrix = importdata(filepath_hic);
atac_pos = importdata(filepath_atac);
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
atac_size = size(atac_pos,1);
atac_idx = zeros(1, atac_size-1);
for i = 1:atac_size
    for s = 1:split
        if atac_pos(i) > bin * (s-1)&& atac_pos(i) < bin * s
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
[IDX, isnoise] = DBSCAN(data_atac,0.015,10);
k = max(IDX);
Colors = hsv(k);
for i = 1:k
    cluster_atac=data_atac(IDX==i, :);
    color = Colors(i,:);    
    scatter3(cluster_atac(:, 1), cluster_atac(:, 2), cluster_atac(:, 3), 10, color);
    hold on
end
hold on
plot3(dataX,dataY,dataZ, "black");
saveas(gcf,output_path);
