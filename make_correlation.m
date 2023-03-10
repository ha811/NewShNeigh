clc;
clear;
close all;

ShNeigh_path = "/Users/hanayaren/Desktop/ShNeigh/ShNeigh_reconmat_40000.txt";
GEMDSh_path = "/Users/hanayaren/Desktop/GEMDSh_out/GEMDSh_reconmat_40000.txt";
dist_path = "/Users/hanayaren/Desktop/ShNeigh/ShNeigh_distmat_40000.txt";
contact_path = "/Users/hanayaren/Desktop/chr4_input/B2_40000_iced_chr4_dense.matrix";
GEMDSh_mat = importdata(GEMDSh_path);
dist_mat = importdata(dist_path);
contact_mat = importdata(contact_path);
ShNeigh_mat = importdata(ShNeigh_path);
ShNeigh_size = size(ShNeigh_mat, 1);
GEMDSh_size = size(GEMDSh_mat, 1);
dist_size = size(dist_mat, 1);
contact_size = size(contact_mat, 1);
alpha = 1;
if dist_size > GEMDSh_size
    start_pos = dist_size - GEMDSh_size +1;
    dist_mat = dist_mat(start_pos:end, start_pos:end);
else
    start_pos = GEMDSh_size - dist_size +1;
    GEMDSh_mat = GEMDSh_mat(start_pos:end, start_pos:end);
end
GEMDSh_recon_contact = (GEMDSh_mat+0.01).^-alpha;
ShNeigh_recon_contact = (ShNeigh_mat+0.01).^-alpha;
R_GEMDSh = corr2(GEMDSh_mat, dist_mat);
R_Sh = corr2(ShNeigh_mat, dist_mat);
if contact_size > GEMDSh_size
    start_contact = contact_size - GEMDSh_size + 1;
    contact_mat = contact_mat(start_contact:end, start_contact:end);
end
R_GEMDSh_contact = corr2(GEMDSh_recon_contact, contact_mat);
R_Sh_contact = corr2(ShNeigh_recon_contact, contact_mat);
y = [R_Sh, R_GEMDSh];
x = [1,2];
bar(x, y)
xticks([1, 2]);
xticklabels({'ShNeigh', 'GEMDSh'})
title("correlation of distance mat")
saveas(gcf, "/Users/hanayaren/Desktop/GEMDSh_out/mm10_chr4_corr_40000_dist.png")

y = [R_Sh_contact, R_GEMDSh_contact];
x = [1, 2];

bar(x, y)
xticks([1, 2]);
xticklabels({'ShNeigh', 'GEMDSh'})
title("correlation of reconstructed contact mat")
saveas(gcf, "/Users/hanayaren/Desktop/GEMDSh_out/mm10_chr4_corr_40000_recon_contact.png")








