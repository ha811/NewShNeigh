clc;
clear;
close all;

FLA_path = "/Users/hanayaren/Desktop/chr21_trimmed_XYZ.txt";
GEMDSh_path = "/Users/hanayaren/Desktop/IMR90_GEMDSh_reconmat_5000.txt";
Shdist_path = "/Users/hanayaren/Desktop/ShNeigh_distmat_IMR90_5000.txt";
ShNeigh_recon_path = "/Users/hanayaren/Desktop/ShNeigh_reconmat_IMR90_5000.txt";
contact_path = "/Users/hanayaren/Desktop/chr21_5kb_KRnormed_dense_FLA.txt";
GEMDSh_mat = importdata(GEMDSh_path);
contact_mat = importdata(contact_path);
Shdist_mat = importdata(Shdist_path);
ShNeigh_recon_mat = importdata(ShNeigh_recon_path);
FLA_mat = importdata(FLA_path);
FLA_size = FLA_mat(end, 1);
FLA_full_mat = zeros(FLA_size, 3);
for i = 1:FLA_size
    if nnz(FLA_mat(:, 1) == i) == 1
        I = find(FLA_mat(:, 1) == i);
        FLA_full_mat(i, 1:3) = FLA_mat(I(1,1), 2:4);
    else
        FLA_full_mat(i, 1:3) = [nan, nan, nan];
    end
end
FLA_dist = pdist(FLA_full_mat);
FLA_dist_mat = squareform(FLA_dist);
I = isnan(FLA_dist_mat);
FLA_dist_mat(I) = 0;

Shdist_size = size(Shdist_mat, 1);



GEMDSh_size = size(GEMDSh_mat, 1);
contact_size = size(contact_mat, 1);
alpha = 1;
GEMDSh_recon_contact = (GEMDSh_mat+0.01).^-alpha;
FLA_recon_contact = (FLA_dist_mat+0.01).^-0.25;
ShNeigh_recon_contact = (ShNeigh_recon_mat+0.01).^-alpha;
R_FLA_contact = corr2(FLA_recon_contact, contact_mat);
if contact_size > GEMDSh_size
    start_contact = contact_size - GEMDSh_size + 1;
    contact_mat = contact_mat(start_contact:end, start_contact:end);
end
R_GEMDSh_contact = corr2(GEMDSh_recon_contact, contact_mat);
R_ShNeigh_contact = corr2(ShNeigh_recon_contact, contact_mat);

if FLA_size > Shdist_size
    start_FLA = FLA_size - Shdist_size +1;
    FLA_dist_mat = FLA_dist_mat(start_FLA:end, start_FLA:end);
end
R_FLA_dist = corr2(FLA_dist_mat, Shdist_mat);
R_GEMDSh_dist = corr2(GEMDSh_mat, Shdist_mat);
R_ShNeigh_dist = corr2(ShNeigh_recon_mat, Shdist_mat);

y = [R_FLA_contact, R_ShNeigh_contact, R_GEMDSh_contact];
x = [1,2,3];
bar(x, y)
xticks([1, 2, 3]);
xticklabels({'FLAMINGO','ShNeigh', 'GEMDSh'})
title("correlation of reconstructed contact mat")
saveas(gcf, "/Users/hanayaren/Desktop/GEMDSh_out/corr_FLA_GEMDSh_5000.png")
hold off

y = [R_FLA_dist, R_ShNeigh_dist, R_GEMDSh_dist];
x = [1, 2, 3];
bar(x, y)
xticks([1, 2, 3])
xticklabels({'FLAMINGO','ShNeigh' 'GEMDSh'})
title("correlation of distance mat")
saveas(gcf, "/Users/hanayaren/Desktop/GEMDSh_out/corr_FLA_GEMDSh_5000_dist.png")
hold off




