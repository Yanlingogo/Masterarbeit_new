%% Only consider 6 time periods
clc
clear
close all
%%Index setting
% bus idx
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
% branch idx
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% gen idx
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
% cost idx
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

mpc = ext2int(loadcase('case33_modified'));
mpc = ext2int(mpc);
load('Pd2_test.mat');
load("Qd2_test.mat");
%% parameters
baseMVA     = mpc.baseMVA;
baseKV      = mpc.bus(1,BASE_KV);
Umax        = mpc.bus(:,VMAX).^2;                                            
Umin        = mpc.bus(:,VMIN).^2;
Pgmin       = mpc.gen(:,PMIN)/baseMVA;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; 
Pgmax       = mpc.gen(:,PMAX)/baseMVA; 
Qgmax       = mpc.gen(:,QMAX)/baseMVA; 
Pramp       = [1;0.1;0.1];
Pd          = mpc.bus(:,PD)/baseMVA; Pd=[Pd Pd2];%Pd = repmat(Pd,1,6);
Qd          = mpc.bus(:,QD)/baseMVA; Qd=[Qd Qd2];%Qd = repmat(Qd,1,6);

id_gen      = mpc.gen(:,GEN_BUS);
id_slack    =  find(mpc.bus(:,BUS_TYPE) == REF);
id_gen_slack  = find(id_gen == id_slack);
id_gen_nslack = find(id_gen ~= id_slack);
Ngen_nslack = numel(id_gen_nslack);
Nbus        = size(mpc.bus,1);
Ngen        = numel(id_gen);
Nbranch     = size(mpc.branch,1);

branch_r   = mpc.branch(:,BR_R);
branch_x   = mpc.branch(:,BR_X);

from_bus       = mpc.branch(:, F_BUS);                         
to_bus         = mpc.branch(:, T_BUS);  
Cf             = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
Ct             = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
C              = Cf - Ct;
Cg             = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);
Cg_nslack      = Cg(:,id_gen_nslack);
Cg_slack       = Cg(:,id_gen_slack);
%% Problem formulation
% Time period
T = 6;
% beta_V
A_v = [eye(Nbus*T) zeros(Nbus*T,2*T*Ngen+2*T*Nbranch);
       -eye(Nbus*T) zeros(Nbus*T,2*T*Ngen+2*T*Nbranch)];
% beta_Generator at non-slack bus 
A_gen_ns = [zeros(Ngen*T,Nbus*T) eye(Ngen*T) zeros(Ngen*T,Ngen*T+2*T*Nbranch);
            zeros(Ngen*T,Nbus*T) -eye(Ngen*T) zeros(Ngen*T,Ngen*T+2*T*Nbranch);
            zeros(Ngen*T,Nbus*T+Ngen*T) eye(Ngen*T) zeros(Ngen*T,2*T*Nbranch);
            zeros(Ngen*T,Nbus*T+Ngen*T) -eye(Ngen*T) zeros(Ngen*T,2*T*Nbranch);];% bounds for P Q
% temporal coupling constraints
C_diff = zeros(1,Ngen+1);
C_diff(1,1) = -1;
C_diff(1,end) = 1;
C_tc = zeros((T-1)*Ngen,Ngen*T);% temporal coupling
for i = 1:(T-1)*Ngen
    C_tc(i,i:i+Ngen) = C_diff;
end
C_tc = vertcat(C_tc,-C_tc);
A_tc = [zeros(2*(T-1)*Ngen,Nbus*T),C_tc,zeros(size(C_tc)),zeros(2*(T-1)*Ngen,2*T*Nbranch)];% no ramp for Q
% beta_voltage_constraints for injection power
C_comb = C;
BR_comb = diag(branch_r);
BX_comb = diag(branch_x);
for i = 2:T
    C_comb = blkdiag(C_comb,C);
    BR_comb = blkdiag(BR_comb,diag(branch_r));
    BX_comb = blkdiag(BX_comb,diag(branch_x));
end
A_inj = [C_comb zeros(Nbranch*T,2*T*Ngen) -2*BR_comb -2*BX_comb;
         -C_comb zeros(Nbranch*T,2*T*Ngen) 2*BR_comb  2*BX_comb];
% beta_power flow balance
Cg_comb = Cg;
C_comb2 = C';
for i = 2:T
    Cg_comb = blkdiag(Cg_comb,Cg);
    C_comb2 = blkdiag(C_comb2,C');
end
A_eq = [zeros(Nbus*T) Cg_comb zeros(Nbus*T, Ngen*T) -C_comb2 zeros(Nbus*T, Nbranch*T);
        zeros(Nbus*T) -Cg_comb zeros(Nbus*T, Ngen*T) C_comb2 zeros(Nbus*T, Nbranch*T);
        zeros(Nbus*T) zeros(Nbus*T, Ngen*T) Cg_comb zeros(T*Nbus, Nbranch*T) -C_comb2;
        zeros(Nbus*T) zeros(Nbus*T, Ngen*T) -Cg_comb zeros(T*Nbus, Nbranch*T) C_comb2;];

% Integrate all matrices
A = vertcat(A_v,A_gen_ns,A_tc,A_inj,A_eq);

% lower and upper bounds on U
b0_U_u = repmat(Umax,T,1);
b0_U_l = -repmat(Umin,T,1);
% lower and upper bounds on P/Q
b0_P_u = repmat(Pgmax,T,1);
b0_P_l = -repmat(Pgmin,T,1);
b0_Q_u = repmat(Qgmax,T,1);
b0_Q_l = -repmat(Qgmin,T,1);
% ramp rate
b0_ramp = repmat(Pramp,2*(T-1),1);
% power flow equation
tol_cons = 1e-6;
b0_inj = tol_cons*ones(2*Nbranch*T,1);
Pd = Pd(:);
Qd = Qd(:);
b0_pf = vertcat(Pd+tol_cons,-Pd+tol_cons,Qd+tol_cons,-Qd+tol_cons);

b0 = vertcat(b0_U_u,b0_U_l,b0_P_u,b0_P_l,b0_Q_u,b0_Q_l,b0_ramp,b0_inj,b0_pf);
%% condense model with umbrella constraints
%[A_u, b_u] = E_UCI(A, b0);
A_u = A;
b_u = b0;
%% find the feasible region
D_all = cell(1,T);
v_all = cell(1,T);
for i = 1:T
    idx_pqs = [1+Ngen*(i-1) 1+Ngen*(i-1)+Ngen*T];
    B = A_u(:,Nbus*T+idx_pqs);
    remainingIndices = setdiff(1:size(A_u,2), Nbus*T+idx_pqs);
    A_de = A_u(:, remainingIndices);
    % generate a large enough space
    D0 = [ 1 0;
     -1 0;
      0 1;
      0 -1;];
    v = [Pgmax(id_gen_slack);
         -Pgmin(id_gen_slack);
         Qgmax(id_gen_slack);
         -Qgmin(id_gen_slack)];
    % parameter for M-method
    M = 80;
    [D0,v] = feas_cut(B,A_de,b_u,D0,v,M);
    D_all{i} = D0;
    v_all{i} = v;
end
%% get the intersection of all linear equations
intersections = cell(1,T);
for i = 1:T
    [D_all{i}, v_all{i}] = E_UCI(D_all{i},v_all{i});
    intersections{i} = intersection(D_all{i}, v_all{i});
end

figure;

% 设置颜色地图，以便每个多边形有不同的颜色
colors = jet(length(intersections)); % 'jet' 是 MATLAB 的一个颜色地图

% 循环通过每个 cell 绘制多边形
for i = 1:length(intersections)
    % 获取当前 cell 的坐标点
    xy = (intersections{i})';
    centroid = mean(xy, 1);
    angles = atan2(xy(:,2) - centroid(2), xy(:,1) - centroid(1));
    [~, order] = sort(angles);
    sortedPoints = xy(order, :);
    % 假设 xy 是一个 n x 2 的矩阵，其中第一列是 x 坐标，第二列是 y 坐标
    x = sortedPoints(:, 1);
    y = sortedPoints(:, 2);
    
    % 在三维空间中添加多边形。Z 坐标由 i 控制，以将它们堆叠起来
    z = ones(size(x)) * i;
    patch(z, x, y, colors(i, :), 'EdgeColor', 'none');
    
    % 可以选择添加一些透明度
    alpha(0.5);

    hold on; % 保持当前图形，以便在上面绘制
    plot3(z, x, y, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i, :));
end

% 调整视角
ylim([-0.2 0.7]);
zlim([-1.2 1.2]);

set(gca, 'YDir','reverse')
view(3);

% 添加轴标签
xlabel('T(time period)');
ylabel('P(Mvar)');
zlabel('Q(Mvar)');
title('Feasible region with temporal coupling');

% 优化图形显示
grid on;
%axis tight;
rotate3d on; % 允许使用鼠标旋转视图