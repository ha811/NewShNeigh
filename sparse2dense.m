clc;
clear;
close all;
filepath_sparse = "/Users/hanayaren/Desktop/chr21_5kb_trimmed_KRnormed_trimmed.RAWobserved";
sparse_mat = importdata(filepath_sparse);
sparse_num = size(sparse_mat, 1);
mat_size = sparse_mat(end, 2); 
dense_mat = zeros(mat_size);
for i = 1:sparse_num
    dense_mat(sparse_mat(i, 1), sparse_mat(i, 2)) = sparse_mat(i, 3);
    dense_mat(sparse_mat(i, 2), sparse_mat(i, 1)) = sparse_mat(i, 3);
end
dense_mat = dense_mat(4000:9620, 4000:9620);
writematrix(dense_mat, "/Users/hanayaren/Desktop/chr21_5kb_KRnormed_dense_FLA.txt")