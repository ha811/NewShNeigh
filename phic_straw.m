clc;
clear;
close all;

filepath_mat = input('入力hicファイルパスを入力してください：', "s");
hic_mat = importdata(filepath_mat);
start_pos = fix(input('input start position:') / 40000);
end_pos = fix(input('input end position:') / 40000);

straw_mat = hic_mat(start_pos:end_pos, start_pos:end_pos);
writematrix(straw_mat, '/Users/hanayaren/Desktop/chr4_input/phic_straw.txt', 'Delimiter', 'space');
M = max(straw_mat, [], 'all')
