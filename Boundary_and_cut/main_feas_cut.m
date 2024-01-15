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
%% parameters
baseMVA     = mpc.baseMVA;
baseKV      = mpc.bus(1,BASE_KV);
Umax        = mpc.bus(:,VMAX).^2;                                            
Umin        = mpc.bus(:,VMIN).^2;
Pgmin       = mpc.gen(:,PMIN)/baseMVA;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; 
Pgmax       = mpc.gen(:,PMAX)/baseMVA; 
Qgmax       = mpc.gen(:,QMAX)/baseMVA; 
Pd          = mpc.bus(:,PD)/baseMVA;
Qd          = mpc.bus(:,QD)/baseMVA;

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
beta = [eye(Nbus) zeros(Nbus,2*Ngen_nslack+2*Nbranch);
        -eye(Nbus) zeros(Nbus,2*Ngen_nslack+2*Nbranch);
        zeros(Ngen_nslack,Nbus) eye(Ngen_nslack) zeros(Ngen_nslack,Ngen_nslack+2*Nbranch);
        zeros(Ngen_nslack,Nbus) -eye(Ngen_nslack) zeros(Ngen_nslack,Ngen_nslack+2*Nbranch);
        zeros(Ngen_nslack,Nbus+Ngen_nslack) eye(Ngen_nslack) zeros(Ngen_nslack,2*Nbranch);
        zeros(Ngen_nslack,Nbus+Ngen_nslack) -eye(Ngen_nslack) zeros(Ngen_nslack,2*Nbranch);

        C zeros(Nbranch,2*Ngen_nslack) diag(-2*branch_r) diag(-2*branch_x);
        -C zeros(Nbranch,2*Ngen_nslack) diag(2*branch_r) diag(2*branch_x);

        zeros(Nbus) Cg_nslack zeros(Nbus, Ngen_nslack) -C' zeros(Nbus, Nbranch);
        zeros(Nbus) -Cg_nslack zeros(Nbus, Ngen_nslack) C' zeros(Nbus, Nbranch);

        zeros(Nbus) zeros(Nbus, Ngen_nslack) Cg_nslack zeros(Nbus, Nbranch) -C';
        zeros(Nbus) zeros(Nbus, Ngen_nslack) -Cg_nslack zeros(Nbus, Nbranch) C';

        zeros(4,Nbus+2*Ngen_nslack+2*Nbranch)];

gamma = [zeros(2*Nbus+4*Ngen_nslack+2*Nbranch,2);
         Cg_slack zeros(size(Cg_slack));
         -Cg_slack zeros(size(Cg_slack));
         zeros(size(Cg_slack)) Cg_slack;
         zeros(size(Cg_slack)) -Cg_slack;
         1 0; -1 0; 0 1; 0 -1];
tol_cons = 1e-6;
b0 = [Umax;-Umin;Pgmax(id_gen_nslack);-Pgmin(id_gen_nslack);
      Qgmax(id_gen_nslack);-Qgmin(id_gen_nslack);
      tol_cons*ones(Nbranch,1);tol_cons*ones(Nbranch,1);
      Pd+tol_cons;-Pd+tol_cons;Qd+tol_cons;-Qd+tol_cons;
      Pgmax(id_gen_slack);-Pgmin(id_gen_slack);
      Qgmax(id_gen_slack);-Qgmin(id_gen_slack)];
% initilization 
T = [ 1 0;
     -1 0;
      0 1;
      0 -1;];
v = [Pgmax(id_gen_slack);
     -Pgmin(id_gen_slack);
     Qgmax(id_gen_slack);
     -Qgmin(id_gen_slack)];

%M = 10*sum(Pd);
M = 80;
P_half = 0.5*(2*sum(Pd)-sum(Pgmin(id_gen_nslack))-sum(Pgmax(id_gen_nslack)));
Q_half = 0.5*(2*sum(Qd)-sum(Qgmin(id_gen_nslack))-sum(Qgmax(id_gen_nslack)));
%z_0 = [P_half;Q_half];
z_0 = [0.2;0];
% tolerance
tol = 1e-3;
[T,v] = feas_cut(gamma, beta, b0, T, v, M);
% res = 20;
% z_p = linspace(-0.2,0.8,res);
% z_q = linspace(-2,1.5,res);
% mesh = zeros(res);
% for i = 1:res
%     for j =1:res
%         z_B = [z_p(i);z_q(j)];
%         [mesh(i,j),~] = Boundart_check(b0, beta, gamma, z_B);
%     end
% end


% K = 1; K_2 = 0;
% while K ~= 0
%     [K, z_s] = Boundary_search(b0, T, v, gamma, beta, M);
% 
%     lambda_l = 0;
%     lambda_u = 1;
%     lambda = 0.5*(lambda_l + lambda_u);
%     if K <= 1e-16
%         break;
%     else
%         while ~(lambda_u - lambda_l <= tol && K_2 > 0)
%             lambda = 0.5*(lambda_l + lambda_u);
%             z_B = lambda*z_s + (1-lambda)*z_0;
%             [K_2, h_s] = Boundart_check(b0, beta, gamma, z_B);
%             if K_2 == 0
%                 lambda_l = lambda;
%             else 
%                 lambda_u = lambda;
%             end
%         end
%         T = [T;-(h_s'*gamma)];
%         v = [v;-(h_s'*b0)];
%     end
% end

% 设置 x 轴和 y 轴的范围
x = linspace(-0.2, 0.7, 400);
y = linspace(-1.2, 1.2, 400);
[X, Y] = meshgrid(x, y);

% 初始化绘图
fig=figure; box on; hold all; set(fig, 'Position', [100, 100, 850, 650]);

% 设置坐标轴范围和网格线间隔
xlim([-0.2, 0.7]);
ylim([-1.2, 1.2]);
xticks(-0.2:0.1:0.7);
yticks(-1.2:0.5:1.2);
grid on;

% 绘制每个不等式定义的线
for i = 1:size(T,1)
    if T(i,1) == 0
        line(xlim, [v(i)/T(i,2) v(i)/T(i,2)], 'Color', 'r');
    elseif T(i,2) == 0
        line([v(i)/T(i,1) v(i)/T(i,1)], ylim, 'Color', 'r');
    else
        plot(x, (v(i) - T(i,1)*x)/T(i,2), 'r');
    end
end

% 检查每个点是否满足所有不等式
% inside = all((A * [X(:), Y(:)]' <= b)');

% 绘制可行区域
% fill(X(inside), Y(inside), 'g', 'FaceAlpha', 0.3);

% 设置坐标轴标签和标题
xlabel('x');
ylabel('y');
title('Linear Inequality Constraints');
hold off;