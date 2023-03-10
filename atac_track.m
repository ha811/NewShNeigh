clc;
clear;
close all;
chr = 1;
filepath_minus = sprintf("/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC/Dox_minus_chr%d_trimmed.sam", chr);
dox_minus = importdata(filepath_minus);
filepath_plus = sprintf("/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC/Dox_plus_chr%d_trimmed.sam", chr);
dox_plus = importdata(filepath_plus);
filepath_notre = sprintf("/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC/NoTre_chr%d_trimmed.sam", chr);
notre = importdata(filepath_notre);
filepath_TSA = sprintf("/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC/TSA_chr%d_trimmed.sam", chr);
TSA = importdata(filepath_TSA);
filepath_triptolide = sprintf("/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC/triptolide_chr%d_trimmed.sam", chr);
triptolide = importdata(filepath_triptolide);
bin_resolution = 40000;
bin_num = 6232;
genome_length = 6232 * 40000;

dox_minus_len = size(dox_minus, 1);
dox_plus_len = size(dox_plus, 1);
notre_len = size(notre, 1);
TSA_len = size(TSA, 1);
triptolide_len = size(triptolide, 1);




minus_score_mat = zeros(bin_num,1);
plus_score_mat = zeros(bin_num,1);
notre_score_mat = zeros(bin_num, 1);
TSA_score_mat = zeros(bin_num, 1);
triptolide_score_mat = zeros(bin_num, 1);
for i = 1:bin_num
    minus_score_mat(i,1) = nnz((dox_minus>(i-1)*bin_resolution) & (dox_minus<i*bin_resolution));
    plus_score_mat(i,1) = nnz((dox_plus>(i-1)*bin_resolution) & (dox_plus<i*bin_resolution));
    notre_score_mat(i,1) = nnz((notre>(i-1)*bin_resolution) & (notre<i*bin_resolution));
    TSA_score_mat(i,1) = nnz((TSA>(i-1)*bin_resolution) & (TSA<i*bin_resolution));
    triptolide_score_mat(i,1) = nnz((triptolide>(i-1)*bin_resolution) & (triptolide<i*bin_resolution));
end


x = linspace(1,genome_length, bin_num);
y = ones(bin_num, 1);
plot(x, y)
hold on

plusminus_score_mat = zeros(bin_num, 1);
TSA_notre_score_mat = zeros(bin_num, 1);
triptolide_notre_score_mat = zeros(bin_num, 1);
triptolide_TSA_score_mat = zeros(bin_num, 1);
for i = 1:bin_num
    plusminus_score_mat(i, 1) = minus_score_mat(i, 1) / plus_score_mat(i, 1);
    TSA_notre_score_mat(i, 1) = TSA_score_mat(i, 1) / notre_score_mat(i, 1);
    triptolide_notre_score_mat(i, 1) = triptolide_score_mat(i, 1) / notre_score_mat(i, 1);
    TSA_triptolide_score_mat(i, 1) = TSA_score_mat(i, 1) / triptolide_score_mat(i, 1);
    triptolide_TSA_score_mat(i, 1) = triptolide_score_mat(i, 1) / TSA_score_mat(i, 1);
end

plot(x, plusminus_score_mat)
hold off
legend('controll', 'minus/plus')




saveas(gcf, sprintf('/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC_track/chr%d_dox.png', chr))
plot(x, y)
hold on
plot(x, TSA_notre_score_mat)
legend('notre', 'TSA')
hold off
saveas(gcf, sprintf('/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC_track/chr%d_TSA.png', chr))
plot(x,y)
hold on
plot(x, triptolide_notre_score_mat)
legend('notre', 'triptolide')
saveas(gcf, sprintf('/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC_track/chr%d_triptolide.png', chr))
hold off
plot(x, y)
hold on
plot(x, TSA_triptolide_score_mat)
ylim([0, 20])
legend('controll', 'TSA/triptolide')
saveas(gcf, sprintf('/Users/hanayaren/Desktop/2022_hirota_ATAC/ATAC_track/chr%d_TSA_triptolide.png', chr))


