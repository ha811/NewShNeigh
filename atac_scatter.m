clc;
clear;
close all;

filepath_hic = ("/Users/hanayaren/Desktop/hg19_40000_dense/hg19_40000_iced_chr3_dense.matrix");
output_path_scat = ("/Users/hanayaren/Desktop/hg19_juicer_fig/hg19_chr3_40000_contact.png");
hic_resolution = 40000;
heatmapMatrix = importdata(filepath_hic);
hic_num = size(heatmapMatrix, 1);
for i = 1:hic_num
    for r = 1:hic_num
        if i > r
            heatmapMatrix(i, r) = 0;
        end
    end
end
[row, col, v] = find(heatmapMatrix>1.3);
row = row * hic_resolution;
col = col * hic_resolution;

fig = scatter(row,col,[],v);
fig.Marker = '.';


saveas(fig,output_path_scat);


