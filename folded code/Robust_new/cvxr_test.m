clc;
clear;
%% Index 
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

%% initilize the problem
mpc = runopf(case9);
% number of all elements
Nbus = size(mpc.bus,1);
Ngen = size(mpc.gen,1);
idx_load = find((mpc.bus(:,PD).^2+mpc.bus(:,QD).^2) ~= 0);
Nload = numel(idx_load);
Npq = sum(mpc.bus(:,BUS_TYPE)==1);
Nbranch = size(mpc.branch,1);

% parameters 
bus = mpc.bus;
gen = mpc.gen;
branch = mpc.branch;

bus_type = bus(:,BUS_TYPE);
idx_slack = find(bus_type == 3);
idx_nslack = find(bus_type ~= 3); Nnslack = numel(idx_nslack);
idx_sgen = find(gen(:,GEN_BUS)==idx_slack);

idx_fr = branch(:,F_BUS);
idx_to = branch(:,T_BUS);
id_gen = gen(:,GEN_BUS); 
id_gen_ns = setdiff(id_gen,idx_slack);
idx_ngen = find(bus(:,BUS_TYPE)~=2); Nngen = numel(idx_ngen);
idx_pq = find(bus_type == 1);
onoff_pq = zeros(Nbus,1); onoff_pq(idx_pq) = 1;
idx_pvs = find(bus_type ~=1);

%gencost = mpc.gencost(:,COST:end);

E_fr = sparse(idx_fr, 1:Nbranch, 1, Nbus, Nbranch);
E_to = sparse(idx_to, 1:Nbranch, 1, Nbus, Nbranch);
E = E_fr - E_to;
Cg = sparse(id_gen,1:Ngen,1,Nbus,Ngen);
Cg_ns = sparse(id_gen_ns,1:Ngen-1,1,Nbus,Ngen-1);
Cl = sparse(idx_load,1:Nload,1,Nbus,Nload);

v0 = bus(:,VM);
theta0 = bus(:,VA);
Phi0 = E'*theta0;
% alpha0 = gen(:,end);
% delta0 = mpc.delta;
pg0 = (gen(:,GEN_STATUS).*gen(:,PG))/mpc.baseMVA;
qg0 = (gen(:,GEN_STATUS).*gen(:,QG))/mpc.baseMVA;
pl0 = bus(idx_load,PD)/mpc.baseMVA;
ql0 = bus(idx_load,QD)/mpc.baseMVA;
p_inj0 = Cg*pg0-Cl*pl0;
q_inj0 = Cg*qg0-Cl*ql0;

% limitation of the network
vmax = bus(:,VMAX);
vmin = bus(:,VMIN);
Phi_max = deg2rad(branch(:,ANGMAX));
Phi_min = deg2rad(branch(:,ANGMIN));
pg_max = (gen(:,GEN_STATUS).*gen(:,PMAX))/mpc.baseMVA;
pg_min = (gen(:,GEN_STATUS).*gen(:,PMIN))/mpc.baseMVA;
qg_max = (gen(:,GEN_STATUS).*gen(:,QMAX))/mpc.baseMVA;
qg_min = (gen(:,GEN_STATUS).*gen(:,QMIN))/mpc.baseMVA;

s_line_max = branch(:,RATE_A)/mpc.baseMVA;

phase_shift = false;
[~,~,~,M,M_line] = makeYbus_cvxr(mpc,phase_shift);
M_eq = M;
J1 = -diag(v0(idx_fr).*v0(idx_to).*sin(E'*theta0))*E';
J2 = diag(v0(idx_to).*cos(E'*theta0))*E_fr'+diag(v0(idx_fr).*cos(E'*theta0))*E_to';
J3 = diag(v0(idx_fr).*v0(idx_to).*cos(E'*theta0))*E';
J4 = diag(v0(idx_to).*sin(E'*theta0))*E_fr'+diag(v0(idx_fr).*sin(E'*theta0))*E_to';
D = diag(v0);
J_psi0 = [zeros(Nbus,2*Nnslack) sum(Cg_ns,2) zeros(Nbus,1);
          zeros(Nbus,2*Nnslack) zeros(Nbus,1) sum(Cg_ns,2);
          J1(:,idx_nslack) J2(:,idx_nslack) zeros(Nbranch,2);
          J3(:,idx_nslack) J4(:,idx_nslack) zeros(Nbranch,2);
          zeros(Nbus,Nnslack) 2*D(:,idx_nslack) zeros(Nbus,2)];
J_psi02 = [zeros(Nbus,2*Nnslack) Cg_ns zeros(Nbus,2);
           zeros(Nbus,2*Nnslack) zeros(Nbus,2) Cg_ns;
           J1(:,idx_nslack) J2(:,idx_nslack) zeros(Nbranch,4);
           J3(:,idx_nslack) J4(:,idx_nslack) zeros(Nbranch,4);
           zeros(Nbus,Nnslack) 2*D(:,idx_nslack) zeros(Nbus,4)];
J_inv_M = (M_eq*J_psi0)\M_eq;
J_inv_M2 = pinv(full(M_eq*J_psi02))*M_eq;
C = blkdiag(E(idx_nslack,:)',eye(Nnslack+4));
A = [-eye(Nbranch+Nnslack+4); eye(Nbranch+Nnslack+4)]*C;
K = A*J_inv_M2;

K_plus  = max(K,zeros(size(K))).*(K>0);
K_minus = min(K,zeros(size(K))).*(K<0);

M_line_plus  = max(M_line, zeros(size(M_line))).*(M_line>0);
M_line_minus = min(M_line, zeros(size(M_line))).*(M_line<0);

% original value of residues and basis function
g_vvcos0 = v0(idx_fr).*v0(idx_to).*cos(Phi0) - onoff_pq(idx_fr).*v0(idx_fr).*v0(idx_to).*cos(Phi0)...
    - v0(idx_fr).*onoff_pq(idx_to).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi0;
g_vvsin0 = v0(idx_fr).*v0(idx_to).*sin(Phi0) - onoff_pq(idx_fr).*v0(idx_fr).*v0(idx_to).*sin(Phi0)...
    - v0(idx_fr).*onoff_pq(idx_to).*v0(idx_to).*sin(Phi0) - v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi0;
g_vv0 = v0.^2-2*onoff_pq.*v0.*v0;
Psi_vvcos0 = v0(idx_fr).*v0(idx_to).*cos(Phi0);
Psi_vvsin0 = v0(idx_fr).*v0(idx_to).*sin(Phi0);
Psi_vv0 = v0.^2;

%% problem formulation
entries{1}  = 1:Nbus;        %v_u
entries{2}  = 1:Nbus;        %v_l
entries{3}  = 1:Nbus;        %Delta_v_u
entries{4}  = 1:Nbus;        %Delta_v_l

entries{5}  = 1:Nbranch;     %Phi_u
entries{6}  = 1:Nbranch;     %Phi_l
entries{7}  = 1:Nbranch;     %Delta_Phi_u
entries{8}  = 1:Nbranch;     %Delta_Phi_l

entries{9}  = 1:Nbus;        %g_u_pinj
entries{10} = 1:Nbus;        %g_l_pinj
entries{11} = 1:Nbus;        %g_u_qinj
entries{12} = 1:Nbus;        %g_l_qinj  
entries{13} = 1:Nbranch;     %g_u_vvcos %定义顺序与 julia不同
entries{14} = 1:Nbranch;     %g_l_vvcos
entries{15} = 1:Nbranch;     %g_u_vvsin
entries{16} = 1:Nbranch;     %g_l_vvsin
entries{17} = 1:Nbus;        %g_u_vv
entries{18} = 1:Nbus;        %g_l_vv

entries{19} = 1:Nbranch;     %Psi_u_vvcos
entries{20} = 1:Nbranch;     %Psi_l_vvcos
entries{21} = 1:Nbranch;     %Psi_u_vvsin
entries{22} = 1:Nbranch;     %Psi_l_vvsin
entries{23} = 1:Nbus;        %Psi_u_vv
entries{24} = 1:Nbus;        %Psi_l_vv

entries{25} = 1:8*Nbranch;   %D_vv_u_uu...
entries{26} = 1:8*Nbranch;   %D_cos_u_u...
entries{27} = 1:4*Nbranch;   %p_line_fr_u...

entries{28} = 1:2*Ngen;      %pg_opt, qg_opt


%% lower & upper bounds for the variables 
lbx = [vmin;vmin;vmin-v0;vmin-v0;Phi_min;Phi_min;Phi_min-Phi0;Phi_min-Phi0;-inf(8*Nbus+28*Nbranch+2*Ngen+2*Nload+5,1)];
ubx = [vmax;vmax;vmax-v0;vmax-v0;Phi_max;Phi_max;Phi_max-Phi0;Phi_max-Phi0; inf(8*Nbus+28*Nbranch+2*Ngen+2*Nload+5,1)];
%% initial state x0
x0 = vertcat(v0,v0,zeros(2*Nbus,1),Phi0,Phi0,zeros(4*Nbus+2*Nbranch,1),...
    g_vvcos0,g_vvcos0,g_vvsin0,g_vvsin0,g_vv0,g_vv0,...
    Psi_vvcos0,Psi_vvcos0,Psi_vvsin0,Psi_vvsin0,Psi_vv0,Psi_vv0,...
    zeros(20*Nbranch,1),pg0,qg0,pl0,ql0,zeros(5,1));
%% equality/inequality constraints
ineq_residues = @(x)create_ineq_residues(x(entries{13}),x(entries{14}),x(entries{15}),x(entries{16}),x(entries{17}),x(entries{18}),...
    K_plus,K_minus,Phi_max,Phi_min,vmax,vmin,idx_nslack);

%% objective function
obj_p = @(x)create_obj_p(entries{1}, vang, Pg, Qg, Gbus, Bbus, Pd, Qd, Cg,id_slack);
obj_q = @(x)create_obj_q(vmag, vang, Pg, Qg, Gbus, Bbus, Pd, Qd, Cg,id_slack);