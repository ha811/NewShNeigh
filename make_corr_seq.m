clc;
clear;
close all;
R_GEMDSh_dist_list = zeros(19, 1);
R_Sh_dist_list = zeros(19, 1);
R_GEMDSh_contact_list = zeros(19, 1);
R_Sh_contact_list = zeros(19, 1);
de_reso = 40000;
for i = 1:19
    ShNeigh_path_form = "/data2/hanakeren/ShNeigh_out/mm10/mm10_chr%d_ShNeigh_reconmat_%d.txt";
    GEMDSh_path_form = "/data2/hanakeren/ShNeigh_out/mm10/mm10_chr%d_GEMDSh_reconmat_%d.txt";
    dist_path_form = "/data2/hanakeren/ShNeigh_out/mm10/mm10_chr%d_GEMDSh_distmat_%d.txt";
    contact_path_form = "/data2/hanakeren/HiC_prac/out_data/hic_results/matrix/mm10/iced/%d/mm10_%d_iced_chr%d_dense.matrix";
    ShNeigh_path = sprintf(ShNeigh_path_form, i, de_reso);
    GEMDSh_path = sprintf(GEMDSh_path_form, i, de_reso);
    dist_path = sprintf(dist_path_form, i, de_reso);
    contact_path = sprintf(contact_path_form, de_reso, i, de_reso);
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
    R_GEMDSh_dist = corr2(GEMDSh_mat, dist_mat);
    R_Sh_dist = corr2(ShNeigh_mat, dist_mat);
    R_GEMDSh_dist_list(i) = R_GEMDSh_dist;
    R_Sh_dist_list(i) = R_Sh_dist;
    if contact_size > GEMDSh_size
        start_contact = contact_size - GEMDSh_size + 1;
        contact_mat = contact_mat(start_contact:end, start_contact:end);
    end
    R_GEMDSh_contact = corr2(GEMDSh_recon_contact, contact_mat);
    R_Sh_contact = corr2(ShNeigh_recon_contact, contact_mat);
    R_GEMDSh_contact_list(i) = R_GEMDSh_contact;
    R_Sh_contact_list(i) = R_Sh_contact;
end

boxplot([R_Sh_dist_list, R_GEMDSh_dist_list], 'Labels', {'ShNeigh', 'New method'})
title("correlation of distance mat")
ylabel('correlation')
xlabel('modeling method')
distboxform = "/data2/hanakeren/ShNeigh_out/mm10/mm10_corr_%d_dist_box.png";
distbox_path = sprintf(distboxform, de_reso);
saveas(gcf, distbox_path)


boxplot([R_Sh_contact_list, R_GEMDSh_contact_list], 'Labels', {'ShNeigh', 'New method'})
title("correlation of reconstructed contact mat")
ylabel('correlation')
xlabel('modeling method')
reconboxform = "/data2/hanakeren/ShNeigh_out/mm10/mm10_corr_%d_recon_contact_box.png";
reconbox_path = sprintf(reconboxform, de_reso);
saveas(gcf, reconbox_path)










