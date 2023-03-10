h = openfig('/Users/hanayaren/Desktop/GEMDSh_out/GEMDSh_40000.fig');
% 動画ファイルの保存
v = VideoWriter('/Users/hanayaren/Desktop/ShNeigh_video/GEMDSh_40000.mp4','MPEG-4');
open(v);
for n = 1:3:360
    view([n, 30]);
    drawnow; % 描画更新(※)
    frame = getframe(h);
    writeVideo(v,frame);
end
close(v);