clc;
clear;
close all;

%filepath_rough = input('低解像度Hi-Cファイルパスを入力してください：', "s");
%filepath_detailed = input('高解像度Hi-Cファイルパスを入力して下さい：', "s");
%output_path = input('出力ファイルパスを入力してください：', "s");
filepath_rough = "/Users/hanayaren/Desktop/chr4_input/mm10_1000000_iced_chr4_dense.matrix";
filepath_detailed = "/Users/hanayaren/Desktop/chr4_input/B2_40000_iced_chr4_dense.matrix";
output_path = "/Users/hanayaren/Desktop/GEMDSh_out/GEMDSh_40000.fig";
rough_Matrix = importdata(filepath_rough);
detailed_Matrix = importdata(filepath_detailed);


rough_reso = 1000000;
de_reso = 40000;
RD = rough_reso / de_reso;
[roughMatrix,retVals, ~]= ShNeigh(rough_Matrix,1);
rough_num = size(roughMatrix, 1);

block_pos_list = zeros(1, 3);
block_sepos_list = zeros(2, 3, rough_num);
scaling_const = 100;
rough_model = roughMatrix(:, 2:4) * scaling_const;
[k, rough_v] = convhull(rough_model);
[detMatrix, retVals, shMat] = ShNeigh(detailed_Matrix, 1);
detailed_num = size(shMat, 1);
block_num = fix(detailed_num / RD) + 1;
rough_model = rough_model(rough_num - block_num + 1:end,:);
start_num = rem(detailed_num, RD);
count = 0;
flat_count=0;
error_count=0;

%粗いモデルの中にブロックごとのモデルをそのまま配置
for i = 1:block_num
    error = 0;
    if i == 1
        block_mat = shMat(1:start_num, 1:start_num);
    else
        block_mat = shMat(start_num + (i-2)*RD + 1:start_num + (i-1)*RD, start_num + (i-2)*RD + 1:start_num + (i-1)*RD);
    end
    if nnz(block_mat) < 100
        continue
    else
        XYZ = mdscale(block_mat, 3);
        count = count + 1;
        genome_length = size(block_mat, 1) * 10000;
        block_model = XYZ;
    end
    
    if (nnz(block_model(:, 1)) == 0) || (nnz(block_model(:, 2)) == 0) || (nnz(block_model(:, 3)) == 0)
        
        block_model(:, 1) = rough_model(count, 1);
        block_model(:, 2) = rough_model(count, 2);
        block_model(:, 3) = rough_model(count, 3);
        flat_count=flat_count+1;
        error = 1;
    elseif anynan(block_model)
        block_model(:, 1) = rough_model(count, 1);
        block_model(:, 2) = rough_model(count, 2);
        block_model(:, 3) = rough_model(count, 3);
        error_count=error_count+1;
        error = 1;

    else
        [k, block_v] = convhull(block_model);
        if (std(block_model(:, 1)) < 0.05) || (std(block_model(:, 2)) < 0.05) || (std(block_model(:, 3)) < 0.05)
            block_model(:, 1) = rough_model(count, 1);
            block_model(:, 2) = rough_model(count, 2);
            block_model(:, 3) = rough_model(count, 3);
            error = 1;
            flat_count=flat_count+1;
        else    
            block_scaling_const = nthroot((rough_v*genome_length)/(rough_num*1000000*block_v), 3);
            block_model(:, 1) = block_model(:, 1)*block_scaling_const + rough_model(count, 1);
            block_model(:, 2) = block_model(:, 2)*block_scaling_const + rough_model(count, 2);
            block_model(:, 3) = block_model(:, 3)*block_scaling_const + rough_model(count, 3);
        end
    end
    
    if count == 1
        block_pos_list = block_model;
        block_sepos_list(1, :, count) = block_model(1, :);
        block_sepos_list(2, :, count) = block_model(end, :);
    elseif error
        block_pos_list = cat(1, block_pos_list, block_model(:, 1:3));

    else
        center_block = rough_model(count, :);
        v1 = block_model(1, :) - center_block;
        v2 = block_pos_list(end, :) - center_block;
        [angle, axis_rotation] = find_angle_axis(v1, v2);
        block_normed = zeros(size(block_model, 1), 4);
        block_normed(:, 1) = block_model(:, 1) - center_block(1, 1);
        block_normed(:, 2) = block_model(:, 2) - center_block(1, 2);
        block_normed(:, 3) = block_model(:, 3) - center_block(1, 3);
        block_normed(:, 4) = 1;
        Rotation_mat = build_Rotation_Matrix(angle, axis_rotation);
        block_normed_trans = block_normed.';
        block_normed_rotated_trans = Rotation_mat * block_normed_trans;
        block_normed_rotated = block_normed_rotated_trans.';
        block_rotated = block_normed_rotated(:, 1:3) + center_block;
        block_pos_list = cat(1, block_pos_list, block_rotated(:, 1:3));
        
    end

    
end






dataX = block_pos_list(:,1);
dataY = block_pos_list(:,2);
dataZ = block_pos_list(:,3);

hold on
plot3(dataX,dataY,dataZ, "blue");
saveas(gcf,output_path);

logmat = [flat_count error_count];
writematrix(logmat, 'GEMSh_log.txt');

dist = pdist(block_pos_list);
dist_mat = squareform(dist);
writematrix(dist_mat, '/Users/hanayaren/Desktop/GEMDSh_out/GEMDSh_reconmat_40000.txt')
writematrix(shMat, '/Users/hanayaren/Desktop/GEMDSh_out/GEMDSh_distmat_40000.txt')
writematrix(block_pos_list, '/Users/hanayaren/Desktop/GEMDSh_out/GEMDSh_XYZ_40000.txt')
