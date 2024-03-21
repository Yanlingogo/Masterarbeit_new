clear;
clc;

opts = mpoption;
opts.opf.violation   = 1e-12;
opts.mips.costtol    = 1e-12;
opts.mips.gradtol    = 1e-12;
opts.mips.comptol    = 1e-12;
opts.opf.ignore_angle_lim = true;
opts.out.all = 0; 

mpc = ext2int(loadcase('case33_modified'));
mpc = runopf(mpc,opts);
mpc = ext2int(mpc);

idx_slack = find(mpc.bus(:,2)==3);
idx_sgen = find(mpc.gen(:,1)==idx_slack);
idx_ngen = find(mpc.gen(:,1)~=idx_slack);

mpc.bus(:,9) = deg2rad(mpc.bus(:,9)); % 角度转弧度

idx_slack = find(mpc.bus(:,2)==3);
plot_rng = [-2 2 -2.2 2.2]; %画图范围
resolution = 30; % 画图分辨率
contour_plot(mpc,idx_slack,plot_rng,resolution);


