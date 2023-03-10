clc;
clear;
close all;
filepath_mm10 = "/Users/hanayaren/Desktop/hg19_40000_iced_chr1_dense.matrix";
filepath_hg19 = "/Users/hanayaren/Desktop/hg19_40000_dense/hg19_40000_iced_chr1_dense.matrix";
mm10_mat = importdata(filepath_mm10);
hg19_mat = importdata(filepath_hg19);
mm10_size = size(mm10_mat, 1);
hg19_size = size(hg19_mat, 1);
mm10_plus_num = nnz(mm10_mat);
hg19_plus_num = nnz(hg19_mat);
mm10_ratio = mm10_plus_num / (mm10_size * mm10_size);
hg19_ratio = hg19_plus_num / (hg19_size * hg19_size);

mm10_avg = sum(mm10_mat, 'all') / mm10_plus_num;
hg19_avg = sum(hg19_mat, 'all') / hg19_plus_num;
