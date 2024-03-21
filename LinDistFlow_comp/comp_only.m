clc
close all
clear;
%% index
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
% cost idx
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
%% load the result
load('vert_org_1.mat')
lightorange = [250, 188, 113]/255;
orange = [255, 99, 0]/255;
h1 = fill(vert(:,1), vert(:,2), lightorange, 'FaceAlpha',0.7);
set(h1, 'facealpha', 0.6, 'EdgeColor', orange,'LineWidth', 2);
hold all

load('AC_org_1.mat')
lightzs = [181, 180, 214]/255;
violett = [158, 49, 252]/255;
h2 = fill(sortedPoints(:,1), sortedPoints(:,2), lightzs);
set(h2, 'facealpha', 0.9, 'EdgeColor', violett,'LineWidth',2);
%% compensation for sampled points: first method
% e_st   = sparse(id_slack,1,1,Nbus,1);
% 
% A = [e_st' zeros(1,2*Nbus);
%      C -2*R -2*X zeros(Nbranch,2);
%      zeros(Nbus) -C' zeros(Nbus,Nbranch) e_st zeros(Nbus,1);
%      zeros(Nbus,Nbus+Nbranch) -C' zeros(Nbus,1) e_st];
% B = [zeros(Nbus,2*(Ngen-1));
%      Cg_ns zeros(Nbus,Ngen-1);
%      zeros(Nbus,Ngen-1), Cg_ns];
% 
% b = [1;zeros(Nbranch,1);Pd;Qd];
% 
% P_p = Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
% Q_p = Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
% H =blkdiag(P_p,Q_p); % participation factor
% 
% u = (-sortedPoints + [sum(Pd) sum(Qd)])*H';
% x = A\(b - B*u');
% U_comp = x(1:Nbus,:);
% Pij_comp = x(Nbus+1:Nbus+Nbranch,:);
% Qij_comp = x(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
% l_comp = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:);
% PQ_loss = [(branch_r'*l_comp)' (branch_x'*l_comp)'];
% 
% comp_pcc = sortedPoints + PQ_loss;
% plot(comp_pcc(:,1),comp_pcc(:,2),'gx');
%% Compensation with z_0, second
% grid the feasible region
pcc_grid = grid_region(vert);

mpc = loadcase('case33_org');
% parameters
id_bus      = mpc.bus(:,BUS_I);
id_gen      = mpc.gen(:,GEN_BUS);
Nbus        = size(mpc.bus,1);
Ngen        = numel(id_gen);
Nbranch     = size(mpc.branch,1);

% idx
id_slack =  find(mpc.bus(:,BUS_TYPE) == REF);
Nslack = numel(id_slack);
gen_nslack = find(id_gen ~= id_slack);
id_gen_nslack = id_gen(gen_nslack);
id_gen_slack = find(id_gen == id_slack);

baseMVA     = mpc.baseMVA;              % baseMVA
cost_param  = mpc.gencost(:,5:end);     % objective coefficients
Umax        = mpc.bus(:,VMAX).^2;                                            
Umin        = mpc.bus(:,VMIN).^2;
Phimax      = mpc.branch(:,ANGMAX)/180*pi;
Phimin      = mpc.branch(:,ANGMIN)/180*pi;
Pgmin       = mpc.gen(:,PMIN)/baseMVA; Pgmin(id_gen_slack) = -10;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; Qgmin(id_gen_slack) = -10;
Pgmax       = mpc.gen(:,PMAX)/baseMVA; Pgmax(id_gen_slack) = 10;
Qgmax       = mpc.gen(:,QMAX)/baseMVA; Qgmax(id_gen_slack) = 10;
Fmax        = mpc.branch(:,RATE_A)/baseMVA;

Pd          = mpc.bus(:,PD)/baseMVA;   
Qd          = mpc.bus(:,QD)/baseMVA;
Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);
Cg_ns       = sparse(id_gen_nslack,1:Ngen-1,ones(Ngen-1,1),Nbus,Ngen-1);

% branch info
branch_r   = mpc.branch(:,BR_R);
branch_x   = mpc.branch(:,BR_X);
R   = diag(mpc.branch(:,BR_R));
X   = diag(mpc.branch(:,BR_X));
Z2  = R.^2 + X.^2; 
% 
[Ybus, Yf, Yt] = makeYbus(mpc);
Gbus           = real(Ybus);
Bbus           = imag(Ybus);
Gf             = real(Yf);
Bf             = imag(Yf);
Gt             = real(Yt);
Bt             = imag(Yt);
from_bus       = mpc.branch(:, F_BUS);                         
to_bus         = mpc.branch(:, T_BUS);  
Cf             = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
Ct             = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
C              = Cf - Ct;

B              = (Cf + Ct)';
A              = [zeros(Nbranch,1) eye(Nbranch)]*B - eye(Nbranch);
IN              = (eye(Nbranch) - A)\eye(size(A));
DR             = (eye(Nbranch) - A)\(A*R);
DX             = (eye(Nbranch) - A)\(A*X);
DX_p = max(DX,zeros(size(DX))).*(DX>0);
DX_m = min(DX,zeros(size(DX))).*(DX<0);

Mp = IN'*R*IN;
Mq = IN'*X*IN;

H_l   = IN'*(2*(R*DR+X*DX) + Z2);

mpc = runpf(mpc); %pg = 0, qg = 0


U_0 = mpc.bus(:,8).^2;
p_pcc0 = mpc.gen(id_gen_slack,PG)/mpc.baseMVA;
q_pcc0 = mpc.gen(id_gen_slack,QG)/mpc.baseMVA;
s_0 = [p_pcc0,q_pcc0];
P = mpc.branch(:,14)/mpc.baseMVA;
Q = mpc.branch(:,15)/mpc.baseMVA;
L = (P.^2+Q.^2)./U_0(from_bus);
z_0 = [U_0; P; Q; p_pcc0; q_pcc0; L];

e_st   = sparse(id_slack,1,1,Nbus,1);
P_p = Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
Q_p = Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
H =blkdiag(P_p,Q_p); 

D    = full([1, zeros(1,Nbus*2+Nbranch);
        C, -2*diag(branch_r),-2*diag(branch_x), zeros(Nbranch,2);
        zeros(Nbus,Nbus),C',zeros(Nbus,Nbranch),-e_st,zeros(Nbus,1);
        zeros(Nbus,Nbus+Nbranch),C',zeros(Nbus,1),-e_st]);
RR  = sparse(diag(branch_r));
XX  = sparse(diag(branch_x));
PP  = sparse(diag(P));
QQ  = sparse(diag(Q));
Ga  = sparse(diag(L.^(-1)));

E   = [sparse(1,Nbranch);RR^2+XX^2;Ct' *RR;Ct'*XX];
F   = [Cf, - 2*PP*Ga, -2*QQ*Ga,sparse(Nbranch,2)];
G   = (PP^2+QQ^2) * Ga^2;

jac_z = [D E;F G]; % jacobian 
jac_u = [zeros(Nbus,2*(Ngen-1)); 
         -Cg_ns zeros(Nbus,Ngen-1);
         zeros(Nbus,Ngen-1) -Cg_ns;
         zeros(Nbranch,2*(Ngen-1))];

delta_s = pcc_grid - s_0;
z_comp = z_0 + jac_z\jac_u*H*delta_s';
U_comp = z_comp(1:Nbus,:);
Pij_comp = z_comp(Nbus+1:Nbus+Nbranch,:);
Qij_comp = z_comp(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
L_comp = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:);
PQ_loss = ([branch_r';branch_x']*L_comp)';
pcc_grid = pcc_grid + PQ_loss;
scatter(pcc_grid(:,1),pcc_grid(:,2),5,'filled');

% validColumns = all(U_comp >= 0.81 & U_comp <= 1.21);
% pcc_grid_1 = pcc_grid(validColumns,:);
% scatter(pcc_grid_1(:,1),pcc_grid_1(:,2),5,'filled');

p_inj = Cg_ns*P_p*delta_s(:,1)'-Pd; p_inj = p_inj(2:end,:);
q_inj = Cg_ns*Q_p*delta_s(:,2)'-Qd; q_inj = q_inj(2:end,:);
U_comp = U_comp(1,:) + Mp*p_inj + Mq*q_inj + H_l*L_comp;  

validColumns = all(U_comp >= 0.81 & U_comp <= 1.21);
pcc_grid_2 = pcc_grid(validColumns,:);

scatter(pcc_grid_2(:,1),pcc_grid_2(:,2),5,'filled');
%% The third method
% mpc = runpf(mpc);
% U_0 = mpc.bus(:,8).^2;
% p_pcc0 = mpc.gen(id_gen_slack,PG)/mpc.baseMVA;
% q_pcc0 = mpc.gen(id_gen_slack,QG)/mpc.baseMVA;
% s_0 = [p_pcc0,q_pcc0];
% P = mpc.branch(:,14)/mpc.baseMVA;
% Q = mpc.branch(:,15)/mpc.baseMVA;
% L0 = (P.^2+Q.^2)./U_0(from_bus);
% z_0 = [U_0; P; Q; p_pcc0; q_pcc0; L0];
% 
% e_st   = sparse(id_slack,1,1,Nbus,1);
% 
% A = [e_st' zeros(1,2*Nbus);
%      C -2*R -2*X zeros(Nbranch,2);
%      zeros(Nbus) -C' zeros(Nbus,Nbranch) e_st zeros(Nbus,1);
%      zeros(Nbus,Nbus+Nbranch) -C' zeros(Nbus,1) e_st];
% B = [zeros(Nbus,2*(Ngen-1));
%      Cg_ns zeros(Nbus,Ngen-1);
%      zeros(Nbus,Ngen-1), Cg_ns];
% 
% b = [1;-Z2*L0;Pd+Ct'*R*L0;Qd+Ct'*X*L0];
% 
% P_p = Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
% Q_p = Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
% H =blkdiag(P_p,Q_p);
% 
% u = (-sortedPoints + [sum(Pd) sum(Qd)])*H';
% 
% x = A\(b - B*u');
% U_comp = x(1:Nbus,:);
% Pij_comp = x(Nbus+1:Nbus+Nbranch,:);
% Qij_comp = x(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
% l_comp = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:);
% PQ_loss = [(branch_r'*l_comp)' (branch_x'*l_comp)'];
% 
% comp_pcc = sortedPoints + PQ_loss;
% plot(comp_pcc(:,1),comp_pcc(:,2),'gx');
%% compensation using hessian
% mpc = runpf(mpc);
% U_0 = mpc.bus(:,8).^2;
% p_pcc0 = mpc.gen(id_gen_slack,PG)/mpc.baseMVA;
% q_pcc0 = mpc.gen(id_gen_slack,QG)/mpc.baseMVA;
% s_0 = [p_pcc0,q_pcc0];
% Pij_0 = mpc.branch(:,14)/mpc.baseMVA;
% Qij_0 = mpc.branch(:,15)/mpc.baseMVA;
% L_0 = (Pij_0.^2+Qij_0.^2)./U_0(from_bus);
% theta_0 = [diag(Pij_0); diag(Qij_0); diag(U_0(2:end))];
% 
% e_st   = sparse(id_slack,1,1,Nbus,1);
% 
% A = [e_st' zeros(1,2*Nbus);
%      C -2*R -2*X zeros(Nbranch,2);
%      zeros(Nbus) -C' zeros(Nbus,Nbranch) e_st zeros(Nbus,1);
%      zeros(Nbus,Nbus+Nbranch) -C' zeros(Nbus,1) e_st];
% B = [zeros(Nbus,2*(Ngen-1));
%      Cg_ns zeros(Nbus,Ngen-1);
%      zeros(Nbus,Ngen-1), Cg_ns];
% 
% b = [1;zeros(Nbranch,1);Pd;Qd];
% 
% P_p = Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
% Q_p = Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
% H =blkdiag(P_p,Q_p); % participation factor
% 
% u = (-sortedPoints + [sum(Pd) sum(Qd)])*H';
% x = A\(b - B*u');
% 
% J = [diag(2*Pij_0./U_0(from_bus)); diag(2*Qij_0./U_0(from_bus)); diag(- (Pij_0.^2+Qij_0.^2)./U_0(from_bus).^2)];
% Hes = [diag(2./U_0(from_bus)) zeros(Nbus-1) diag((-2*Pij_0)./(U_0(from_bus).^2));
%        zeros(Nbus-1) diag(2./U_0(from_bus)) diag((-2*Qij_0)./(U_0(from_bus).^2));
%        diag((-2*Pij_0)./(U_0(from_bus).^2)) diag((-2*Qij_0)./(U_0(from_bus).^2)) diag(2*(Pij_0.^2+Qij_0.^2)./(U_0(from_bus).^3))];
% 
% l_comp = zeros(Nbranch,size(sortedPoints,1));
% for i = 1: size(sortedPoints,1)
%     x_index = x(:,i);
%     U_comp = x_index(1:Nbus,:);
%     Pij_comp = x_index(Nbus+1:Nbus+Nbranch,:);
%     Qij_comp = x_index(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
%     theta_comp = [diag(Pij_comp);diag(Qij_comp);diag(U_comp(2:end,:))];
% 
%     l_comp(:,i) = diag(diag(L_0) + J'*(theta_comp-theta_0) + 0.5*(theta_comp-theta_0)'*Hes*(theta_comp-theta_0));
%     % l_comp(:,i) = diag(diag(L_0) + J'*(theta_comp-theta_0));
% end
% PQ_loss = [(branch_r'*l_comp)' (branch_x'*l_comp)'];
% 
% comp_pcc = sortedPoints + PQ_loss;
% plot(comp_pcc(:,1),comp_pcc(:,2),'gx');
%% Newton-Raphson Method
% e_st   = sparse(id_slack,1,1,Nbus,1);
% 
% A = [e_st' zeros(1,2*Nbus);
%      C -2*R -2*X zeros(Nbranch,2);
%      zeros(Nbus) -C' zeros(Nbus,Nbranch) e_st zeros(Nbus,1);
%      zeros(Nbus,Nbus+Nbranch) -C' zeros(Nbus,1) e_st];
% B = [zeros(Nbus,2*(Ngen-1));
%      Cg_ns zeros(Nbus,Ngen-1);
%      zeros(Nbus,Ngen-1), Cg_ns];
% 
% b = [1;zeros(Nbranch,1);Pd;Qd];
% 
% P_p = Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
% Q_p = Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
% H =blkdiag(P_p,Q_p); % participation factor
% 
% u = (-sortedPoints + [sum(Pd) sum(Qd)])*H';
% x = A\(b - B*u');
% U_comp = x(1:Nbus,:);
% comp_pcc = zeros(size(sortedPoints));
% for i = 1:size(sortedPoints,1)
%     mpc.gen(2:3,PG) = u(i,1:2)*10;
%     mpc.gen(2:3,QG) = u(i,3:4)*10;
%     mpc.gen(:,VG) = sqrt(U_comp(id_gen,i));
% 
%     mpc = runpf(mpc);
%     if mpc.success
%         comp_pcc(i,:) = mpc.gen(1,[PG QG]);
%     else
%         comp_pcc(i,:) = [inf inf];
%     end
% end 
% infRows = any(isinf(comp_pcc), 2);
% 
% comp_pcc(infRows, :) = [];
% comp_pcc = comp_pcc/mpc.baseMVA;
% plot(comp_pcc(:,1),comp_pcc(:,2),'gx');