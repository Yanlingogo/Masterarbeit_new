clear;
clc;
mpc = loadcase('case33_modified.m');
mpc = runopf(mpc);
%% Index 
% bus idx
tic
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
%% constant parameters
id_bus      = mpc.bus(:,BUS_I);
id_gen      = mpc.gen(:,GEN_BUS);
Nbus        = size(mpc.bus,1);
Ngen        = numel(id_gen);
Nbranch     = size(mpc.branch,1);

id_sline = 1;

R   = diag(mpc.branch(:,BR_R));
X   = diag(mpc.branch(:,BR_X));
Z2  = R.^2 + X.^2; 

bus = mpc.bus;
gen = mpc.gen;
branch =mpc.branch;
V_max = bus(:,VMAX).^2;
V_min = bus(:,VMIN).^2;
pg_max = (gen(:,GEN_STATUS).*gen(:,PMAX))/mpc.baseMVA;
pg_min = (gen(:,GEN_STATUS).*gen(:,PMIN))/mpc.baseMVA;
qg_max = (gen(:,GEN_STATUS).*gen(:,QMAX))/mpc.baseMVA;
qg_min = (gen(:,GEN_STATUS).*gen(:,QMIN))/mpc.baseMVA;

Cg          = sparse(id_gen,1:Ngen,1,Nbus,Ngen);
Cg          = Cg(2:end,2:end); % without the slack generator
Csline      = sparse(1, id_sline, 1,1,Nbranch);
Pd          = mpc.bus(2:end,PD)/mpc.baseMVA;   
Qd          = mpc.bus(2:end,QD)/mpc.baseMVA;

from_bus       = mpc.branch(:, F_BUS);                         
to_bus         = mpc.branch(:, T_BUS); 
Cf             = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
Ct             = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
B              = (Cf + Ct)';
A              = [zeros(Nbranch,1) eye(Nbranch)]*B - eye(Nbranch);
C              = (eye(Nbranch) - A)\eye(size(A));
DR             = (eye(Nbranch) - A)\(A*R);
DX             = (eye(Nbranch) - A)\(A*X);
DX_p = max(DX,zeros(size(DX))).*(DX>0);
DX_m = min(DX,zeros(size(DX))).*(DX<0);

Mp = C'*R*C;
Mq = C'*X*C;

H   = C'*(2*(R*DR+X*DX) + Z2);
H_p = max(H,zeros(size(H))).*(H>0);
H_m = min(H,zeros(size(H))).*(H<0);
% parameters of current position
V_0 = bus(:,VM).^2;
pg_0 = gen(:,PG)/mpc.baseMVA;
qg_0 = gen(:,QG)/mpc.baseMVA;
P_0 = branch(:,PT)/mpc.baseMVA;
Q_0 = branch(:,QT)/mpc.baseMVA;
l_0 = (P_0.^2+Q_0.^2)./V_0(to_bus);

p_0 = Cg*pg_0(2:end) - Pd;
q_0 = Cg*qg_0(2:end) - Qd;
V_nslack0 = 1 + Mp*p_0+Mq*q_0 - H*l_0;
P_test = C*p_0 - DR*l_0; % to_bus
% Jacobian and Hessian
J = [diag(2*P_0/V_0(from_bus)); diag(2*Q_0/V_0(from_bus)); diag(- (P_0.^2+Q_0.^2)/V_0(from_bus).^2)];
J_p = max(J,zeros(size(J))).*(J>0);
J_m = min(J,zeros(size(J))).*(J<0);
Hes = [diag(2./V_0(from_bus)) zeros(Nbus-1) diag((-2*P_0)./(V_0(from_bus).^2));
       zeros(Nbus-1) diag(2./V_0(from_bus)) diag((-2*Q_0)./(V_0(from_bus).^2));
       diag((-2*P_0)./(V_0(from_bus).^2)) diag((-2*Q_0)./(V_0(from_bus).^2)) diag(2*(P_0.^2+Q_0.^2)./(V_0(from_bus).^3))];

%% Define the variables
pg = sdpvar(Ngen-1,1);
qg = sdpvar(Ngen-1,1);
l_u = sdpvar(Nbranch,1);
l_l = sdpvar(Nbranch,1);

p_inj = Cg*pg - Pd;
q_inj = Cg*qg - Qd;

V_u = (V_0(1)*ones(Nbus-1,1) + Mp*p_inj+Mq*q_inj-H_p*l_l-H_m*l_u);
V_l = (V_0(1)*ones(Nbus-1,1) + Mp*p_inj+Mq*q_inj-H_p*l_u-H_m*l_l);
P_u = (C*p_inj-DR*l_l);
P_l = (C*p_inj-DR*l_u);
Q_u = (C*q_inj -DX_p*l_l-DX_m*l_u);
Q_l = (C*q_inj -DX_p*l_u-DX_m*l_l);

theta_1 = [P_u - P_0; Q_u - Q_0; V_u(from_bus) - V_0(from_bus)];
theta_2 = [P_u - P_0; Q_u - Q_0; V_l(from_bus) - V_0(from_bus)];
theta_3 = [P_u - P_0; Q_l - Q_0; V_u(from_bus) - V_0(from_bus)];
theta_4 = [P_u - P_0; Q_l - Q_0; V_l(from_bus) - V_0(from_bus)];
theta_5 = [P_l - P_0; Q_u - Q_0; V_u(from_bus) - V_0(from_bus)];
theta_6 = [P_l - P_0; Q_u - Q_0; V_l(from_bus) - V_0(from_bus)];
theta_7 = [P_l - P_0; Q_l - Q_0; V_u(from_bus) - V_0(from_bus)];
theta_8 = [P_l - P_0; Q_l - Q_0; V_l(from_bus) - V_0(from_bus)];

cons = -l_u + abs(2*(J_p'*theta_1 + J_m'*theta_8))+l_0 <= 0;
cons = [cons, -l_u + theta_1'*Hes*theta_1+l_0 <= 0;
              -l_u + theta_2'*Hes*theta_2+l_0 <= 0;
              -l_u + theta_3'*Hes*theta_3+l_0 <= 0;
              -l_u + theta_4'*Hes*theta_4+l_0 <= 0;
              -l_u + theta_5'*Hes*theta_5+l_0 <= 0;
              -l_u + theta_6'*Hes*theta_6+l_0 <= 0;
              -l_u + theta_7'*Hes*theta_7+l_0 <= 0;
              -l_u + theta_8'*Hes*theta_8+l_0 <= 0;];
cons = [cons, l_l - (J_p'*theta_8 + J_m'*theta_1+l_0) <= 0];
cons = [cons,V_u <= V_max(2:end);
             V_l >= V_min(2:end);
             l_l >= 0];
cons = [cons,pg <= pg_max(2:end), pg >= pg_min(2:end),...
    qg <= qg_max(2:end), qg >= qg_min(2:end)];

delta_angle = deg2rad(20);

Obj_0 = sum(pg);
options = sdpsettings('solver', 'gurobi');
sol_0 = optimize(cons, Obj_0, options);

Obj_1 = P_u(1)+Q_u(1);
options = sdpsettings('solver', 'gurobi');
sol_1 = optimize(cons, Obj_1, options);

result_Piju = value(P_u);
result_Qiju = value(Q_u);

P1 = [-result_Piju(1),-result_Qiju(1)];


Obj_2 = P_l(1)+Q_l(1);
options = sdpsettings('solver', 'gurobi');
sol_2 = optimize(cons, -Obj_2, options);

result_Pijl = value(P_l);
result_Qijl = value(Q_l);

P2 = [-result_Pijl(1),-result_Qijl(1)];

% 计算矩形的左下角坐标
leftBottomX = min(P1(1), P2(1));
leftBottomY = min(P1(2), P2(2));

% 计算矩形的宽度和高度
width = abs(P1(1) - P2(1));
height = abs(P1(2) - P2(2));

% 画出矩形
figure; % 打开一个新的图形窗口
rectangle('Position', [leftBottomX, leftBottomY, width, height], 'EdgeColor', 'r');
axis equal; % 保持横纵轴的比例相同
grid on; % 显示网格
xlabel('X');
ylabel('Y');
title('Rectangle from Two Points');