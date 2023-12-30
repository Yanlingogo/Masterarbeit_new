clc;
clear;
mpc = loadcase('case9');

mpc = opf_initialization(mpc,0);
mpc.uncertainty.gamma0 = 0.2;

idx_slack = find(mpc.bus(:,2)==3);
idx_sgen = find(mpc.gen(:,1)==idx_slack);
idx_ngen = find(mpc.gen(:,1)~=idx_slack);

alpha_p = mpc.gen(idx_ngen,9)'/sum(mpc.gen(idx_ngen,9));
alpha_q = mpc.gen(idx_ngen,4)'/sum(mpc.gen(idx_ngen,4));
mpc.gen([idx_sgen,idx_ngen'],end+1) = [0 alpha_p]; % 分配slack bus 不平衡量给PV节点的比例
mpc.gen([idx_sgen,idx_ngen'],end+1) = [0 alpha_q];
mpc.gen(:,10) = -mpc.gen(:,9)*0.2; % PV节点PMIN

idx_slack = find(mpc.bus(:,2)==3);
plot_rng = [-2 3 -1 1]; %画图范围
resolution = 20; % 画图分辨率
contour_plot(mpc,idx_slack,plot_rng,resolution);