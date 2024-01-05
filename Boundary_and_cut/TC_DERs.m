%% Only consider 6 time periods

% Data_preprocessing;
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
% load('Pd2_test.mat');
% load("Qd2_test.mat");
%% parameters
% Time period
T = 6;
% network parameters
baseMVA     = mpc.baseMVA;
baseKV      = mpc.bus(1,BASE_KV);
Umax        = mpc.bus(:,VMAX).^2;                                            
Umin        = mpc.bus(:,VMIN).^2;
Pgmin       = mpc.gen(:,PMIN)/baseMVA;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; 
Pgmax       = mpc.gen(:,PMAX)/baseMVA; 
Qgmax       = mpc.gen(:,QMAX)/baseMVA; 

Pgmin(1)    = -5; Pgmax(1)= 5; Qgmin(1)=-5;Qgmax(1) = 5;

Pramp       = [1;0.1;0.1];
% load data from case file
% Pd          = mpc.bus(:,PD)/baseMVA; Pd=[Pd Pd2];
% Qd          = mpc.bus(:,QD)/baseMVA; Qd=[Qd Qd2]; 
% load data generated by rand function
Pd          = Lp;
Qd          = Lq;

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

% parameters for EES
id_ess = [16;26];
Ness = numel(id_ess);
P_ess_max = [0.05;0.03];
P_ess_min = [-0.05;-0.03];
E_ess_max = [0.5;0.5];
E_ess_min = [0;0];
E_ess_0   = [0.05;0.05];
C_ess = sparse(id_ess,1:Ness,ones(Ness,1),Nbus,Ness);

% parameters for uncertainty/DER(Wind, PV),Loads
gamma = 0; % radius γ ∈ R determine the shape and size of confidence ellipsoid
cov_wind = scaled_cov(1,:);
cov_pv = scaled_cov(2,:);
cov_lp = scaled_cov(3,:);
cov_lq = scaled_cov(4,:);

id_wind = [25 30];
Nwind = numel(id_wind);
C_wind = sparse(id_wind,1:Nwind,ones(Nwind,1),Nbus,Nwind);
id_pv = [11 19];
Npv   = numel(id_pv);
C_pv = sparse(id_pv,1:Npv,ones(Npv,1),Nbus,Npv);

Pmax_wind =[];
Pmax_pv = [];
Pd_max = [];
Qd_max = [];

for i = 1:T
    Pmax_wind = [Pmax_wind; WT_power(:,i)-gamma*norm(cov_wind{i},2)];
    Pmax_pv = [Pmax_pv; PV_power(:,i)-gamma*norm(cov_pv{i},2)];
    Pd_max(:,i) = Pd(:,i)+gamma*norm(cov_lp{i},2);
    Qd_max(:,i) = Qd(:,i)+gamma*norm(cov_lq{i},2);
end
% ensure the output from Wind/PV are non-negative
Pmax_wind(Pmax_wind<0) = 0;
Pmax_pv(Pmax_pv<0) = 0;
Pd_max = [zeros(1,T);Pd_max];
Qd_max = [zeros(1,T);Qd_max];
%% Problem formulation
% beta_V
A_v = [eye(Nbus*T) zeros(Nbus*T,2*T*(Ngen+Nbranch+Ness)+T*(Nwind+Npv));
       -eye(Nbus*T) zeros(Nbus*T,2*T*(Ngen+Nbranch+Ness)+T*(Nwind+Npv))];
% beta_Generator at non-slack bus 
A_gen_ns = [zeros(Ngen*T,Nbus*T) eye(Ngen*T) zeros(Ngen*T,(Nwind+Npv+Ngen)*T+2*T*(Nbranch+Ness));
            zeros(Ngen*T,Nbus*T) -eye(Ngen*T) zeros(Ngen*T,(Nwind+Npv+Ngen)*T+2*T*(Nbranch+Ness));
            zeros(Ngen*T,(Nbus+Ngen)*T) eye(Ngen*T) zeros(Ngen*T,(Nwind+Npv)*T+2*T*(Nbranch+Ness));
            zeros(Ngen*T,(Nbus+Ngen)*T) -eye(Ngen*T) zeros(Ngen*T,(Nwind+Npv)*T+2*T*(Nbranch+Ness));];
            % bounds for P Q
% temporal coupling constraints
C_diff = zeros(1,Ngen+1);
C_diff(1,1) = -1;
C_diff(1,end) = 1;
C_tc = zeros((T-1)*Ngen,Ngen*T);% temporal coupling
for i = 1:(T-1)*Ngen
    C_tc(i,i:i+Ngen) = C_diff;
end
C_tc = vertcat(C_tc,-C_tc);
% no ramp for Q
A_ramp = [zeros(2*(T-1)*Ngen,Nbus*T),C_tc,...
    zeros(2*(T-1)*Ngen,(Nwind+Npv+Ngen)*T+2*T*(Nbranch+Ness))];% [U,rampP,Q,flow,DER,ESS)
% beta_voltage_constraints for injection power
C_comb = C;
BR_comb = diag(branch_r);
BX_comb = diag(branch_x);
for i = 2:T
    C_comb = blkdiag(C_comb,C);
    BR_comb = blkdiag(BR_comb,diag(branch_r));
    BX_comb = blkdiag(BX_comb,diag(branch_x));
end
A_inj = [C_comb zeros(Nbranch*T,2*T*Ngen) -2*BR_comb -2*BX_comb zeros(Nbranch*T,2*T*Ness+T*(Nwind+Npv));
         -C_comb zeros(Nbranch*T,2*T*Ngen) 2*BR_comb  2*BX_comb zeros(Nbranch*T,2*T*Ness+T*(Nwind+Npv))];

%% beta_ESS
A_aux   = zeros(2*Ness, (Nbus+2*Ngen+2*Nbranch+Nwind+Npv+Ness)*T); % Auxiliary matrix
% Matrix for end
A_esse  = zeros(2*Ness,Ness*T); 
A_esse(1:Ness,Ness*(T-1)+1:Ness*T) = eye(Ness);
A_esse(Ness+1:2*Ness,Ness*(T-1)+1:Ness*T) = -eye(Ness);

A_essbe = [A_aux A_esse];
% Matrix for charging/discharging
E_comb = zeros(T*Ness,T*Ness);
E_diff_1 = eye(Ness); % note the difference at the E(0)
E_diff = [-eye(Ness) eye(Ness)];
for t = 1:T
    if t == 1
        E_comb(1:Ness,1:Ness) = E_diff_1;
    else
        E_comb((t-1)*Ness+1:t*Ness,(t-2)*Ness+1:t*Ness) = E_diff;
    end
end
A_E_comb = [eye(Ness*T) E_comb];
A_E_comb = [A_E_comb;-A_E_comb];
A_ess = [zeros(Ness*T,2*T*(Ngen+Nbranch)+(Nbus+Nwind+Npv)*T) eye(Ness*T) zeros(Ness*T);
         zeros(Ness*T,2*T*(Ngen+Nbranch)+(Nbus+Nwind+Npv)*T) -eye(Ness*T) zeros(Ness*T);
         zeros(Ness*T,2*T*(Ngen+Nbranch)+(Nbus+Nwind+Npv+Ness)*T) eye(Ness*T);
         zeros(Ness*T,2*T*(Ngen+Nbranch)+(Nbus+Nwind+Npv+Ness)*T) -eye(Ness*T);
         A_essbe;
         zeros(2*T*Ness,2*T*(Ngen+Nbranch)+(Nbus+Nwind+Npv)*T) A_E_comb];
%% DERs/ wind and pv
A_wind = [zeros(Nwind*T,(Nbus+2*Ngen+2*Nbranch)*T) eye(Nwind*T) zeros(Nwind*T,(Npv+2*Ness)*T);
          zeros(Nwind*T,(Nbus+2*Ngen+2*Nbranch)*T) -eye(Nwind*T) zeros(Nwind*T,(Npv+2*Ness)*T);];
A_pv = [zeros(Npv*T,(Nbus+2*Ngen+2*Nbranch+Nwind)*T) eye(Npv*T) zeros(Npv*T,2*Ness*T);
        zeros(Npv*T,(Nbus+2*Ngen+2*Nbranch+Nwind)*T) -eye(Npv*T) zeros(Npv*T,2*Ness*T);];
A_DERs = vertcat(A_wind,A_pv);
%% beta_power flow balance
Cg_comb = Cg;
C_comb2 = C'; % constant matrix for flow variables 
Cess_comb = C_ess;
Cw_comb = C_wind;
Cp_comb = C_pv;
for i = 2:T
    Cg_comb = blkdiag(Cg_comb,Cg);
    C_comb2 = blkdiag(C_comb2,C');
    Cess_comb = blkdiag(Cess_comb,C_ess);
    Cw_comb = blkdiag(Cw_comb,C_wind);
    Cp_comb = blkdiag(Cp_comb,C_pv);
end
A_eq = [zeros(Nbus*T) Cg_comb zeros(Nbus*T, Ngen*T) -C_comb2 zeros(Nbus*T, Nbranch*T) Cw_comb Cp_comb Cess_comb zeros(Nbus*T,Ness*T);
        zeros(Nbus*T) -Cg_comb zeros(Nbus*T, Ngen*T) C_comb2 zeros(Nbus*T, Nbranch*T) -Cw_comb -Cp_comb -Cess_comb zeros(Nbus*T,Ness*T);
        zeros(Nbus*T, (Nbus+Ngen)*T) Cg_comb zeros(T*Nbus, Nbranch*T) -C_comb2 zeros(T*Nbus, (Nwind+Npv+2*Ness)*T);
        zeros(Nbus*T, (Nbus+Ngen)*T) -Cg_comb zeros(T*Nbus, Nbranch*T) C_comb2 zeros(T*Nbus, (Nwind+Npv+2*Ness)*T);];

% Integrate all matrices
A = vertcat(A_v,A_gen_ns,A_ramp,A_inj,A_ess,A_DERs,A_eq);
%% bounds on constraints
tol_cons = 1e-6;
% lower and upper bounds on U
b0_U = vertcat(repmat(Umax,T,1), -repmat(Umin,T,1));

% lower and upper bounds on P/Q
b0_P = vertcat(repmat(Pgmax,T,1),-repmat(Pgmin,T,1));
b0_Q = vertcat(repmat(Qgmax,T,1),-repmat(Qgmin,T,1));

% lower and upper bounds on Wind Turbine
b0_wind = vertcat(Pmax_wind+tol_cons,zeros(Nwind*T,1)+tol_cons);
% lower and upper bounds on PV
b0_pv = vertcat(Pmax_pv+tol_cons,zeros(Npv*T,1)+tol_cons);

% ramp rate
b0_ramp = repmat(Pramp,2*(T-1),1);

% lower and upper bounds on ESS
b0_Pess_u = repmat(P_ess_max,T,1);
b0_Pess_l = -repmat(P_ess_min,T,1);
b0_E_u    = repmat(E_ess_max,T,1);
b0_E_l    = -repmat(E_ess_min,T,1);
b0_end_u  = E_ess_0+tol_cons;
b0_end_l  = -E_ess_0+tol_cons;
b0_tc_u   = [E_ess_0+tol_cons; zeros((T-1)*Ness,1)+tol_cons];
b0_tc_l   = [-E_ess_0+tol_cons; zeros((T-1)*Ness,1)+tol_cons];
b0_tc     = [b0_tc_u;b0_tc_l]; 

b0_ESS    = vertcat(b0_Pess_u,b0_Pess_l,b0_E_u,b0_E_l,b0_end_u,b0_end_l,b0_tc);

% power flow equation
b0_inj = tol_cons*ones(2*Nbranch*T,1);
Pd = Pd_max(:);
Qd = Qd_max(:);
% Pd = Pd(:);
% Qd = Qd(:);
b0_pf = vertcat(Pd+tol_cons,-Pd+tol_cons,Qd+tol_cons,-Qd+tol_cons);

% Integrate all vectors
b0 = vertcat(b0_U,b0_P,b0_Q,b0_ramp,b0_inj,b0_ESS,b0_wind,b0_pv,b0_pf);
%% condense model with umbrella constraints
%[A_u, b_u] = E_UCI(A, b0);
A_u = A;
b_u = b0;
%% find the feasible region
D_all = cell(1,T);
v_all = cell(1,T);
for i = 1:T
    disp(['searching boundaries at T = ',num2str(i)])
    idx_pqs = [1+Ngen*(i-1) 1+Ngen*(T+i-1)];
    B = A_u(:,Nbus*T+idx_pqs);
    remainingIndices = setdiff(1:size(A_u,2), Nbus*T+idx_pqs);
    A_de = A_u(:, remainingIndices);
    % generate a large enough space
    D0 = [ 1 0;
     -1 0;
      0 1;
      0 -1;];
    Ps_max = sum(Pd(Nbus*(T-1)+1:Nbus*T));
    Qs_max = sum(Qd(Nbus*(T-1)+1:Nbus*T));
    Ps_min = sum(Pd(Nbus*(T-1)+1:Nbus*T))-sum(Pgmax(id_gen_nslack))-sum(WT_power(:,i))-sum(PV_power(:,i));
    Qs_min = sum(Qd(Nbus*(T-1)+1:Nbus*T))-sum(Qgmax(id_gen_nslack));
    v = [Ps_max*2;
         -Ps_min*2;
         Qs_max*2;
         -Qs_min*2];
    % parameter for M-method
    M = 80;
    [D0,v] = feas_cut(B,A_de,b_u,D0,v,M);

    D_all{i} = D0;
    v_all{i} = v;
end
%% get the intersection of all linear equations
intersections = cell(1,T);
for i = 1:T
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
ylim([-2 2]);
zlim([-2 2]);

set(gca, 'YDir','reverse')
view(3);

% 添加轴标签
xlabel('T(time period)');
ylabel('P(p.u.)');
zlabel('Q(p.u.)');
title(['Feasible region with temporal coupling (\gamma=', num2str(gamma), ')']);

% 优化图形显示
grid on;
axis normal;
rotate3d on; % 允许使用鼠标旋转视图