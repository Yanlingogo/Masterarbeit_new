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

Mp = 2*C'*R*C;
Mq = 2*C'*X*C;

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
P_test = C*p_0 - DR*l_0; % = P_ij(to_bus)
% Jacobian and Hessian
J = [diag(2*P_0/V_0(from_bus)); diag(2*Q_0/V_0(from_bus)); diag(- (P_0.^2+Q_0.^2)/V_0(from_bus).^2)];
J_p = max(J,zeros(size(J))).*(J>0);
J_m = min(J,zeros(size(J))).*(J<0);
Hes = [diag(2./V_0(from_bus)) zeros(Nbus-1) diag((-2*P_0)./(V_0(from_bus).^2));
       zeros(Nbus-1) diag(2./V_0(from_bus)) diag((-2*Q_0)./(V_0(from_bus).^2));
       diag((-2*P_0)./(V_0(from_bus).^2)) diag((-2*Q_0)./(V_0(from_bus).^2)) diag(2*(P_0.^2+Q_0.^2)./(V_0(from_bus).^3))];
%% 
import casadi.*

pg = SX.sym('pg',Ngen,1);
U = SX.sym('U', Ngen,1);


