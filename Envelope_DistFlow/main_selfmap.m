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
P_test = C*p_0 - DR*l_0; % to_bus
% Jacobian and Hessian
J = [diag(2*P_0/V_0(from_bus)); diag(2*Q_0/V_0(from_bus)); diag(- (P_0.^2+Q_0.^2)/V_0(from_bus).^2)];
J_p = max(J,zeros(size(J))).*(J>0);
J_m = min(J,zeros(size(J))).*(J<0);
Hes = [diag(2./V_0(from_bus)) zeros(Nbus-1) diag((-2*P_0)./(V_0(from_bus).^2));
       zeros(Nbus-1) diag(2./V_0(from_bus)) diag((-2*Q_0)./(V_0(from_bus).^2));
       diag((-2*P_0)./(V_0(from_bus).^2)) diag((-2*Q_0)./(V_0(from_bus).^2)) diag(2*(P_0.^2+Q_0.^2)./(V_0(from_bus).^3))];
%% variables


% p_pcc = SX.sym('p_pcc',1,1);
% q_pcc = SX.sym('q_pcc',1,1);
% pg = SX.sym('pg',Ngen-1,1);
% qg = SX.sym('qg',Ngen-1,1);
% V_u = SX.sym('V_u',Nbranch,1);
% V_l = SX.sym('V_l',Nbranch,1);
% P_u = SX.sym('P_u',Nbranch,1);
% P_l = SX.sym('P_l',Nbranch,1);
% Q_u = SX.sym('Q_u',Nbranch,1);
% Q_l = SX.sym('Q_l',Nbranch,1);
% l_u = SX.sym('l_u',Nbranch,1);
% l_l = SX.sym('l_l',Nbranch,1);
% x = vertcat(pg,qg,V_u,V_l,P_u,P_l,Q_u,Q_l,l_u,l_l);
entries{1}  = 1:2;                                                % p_pcc,q_pcc
entries{2}  = 3:2*Ngen;                                           % pg,qg
entries{3}  = 2*Ngen+1:2*Ngen+Nbus;                               % V_u voltages for buses
entries{4}  = 2*Ngen + Nbus+1:2*Ngen+2*Nbus;                      % V_l
entries{5}  = 2*Ngen+2*Nbus+1:2*Ngen+2*Nbus+Nbranch;              % P_u(to_bus) active power for branch
entries{6}  = 2*Ngen+2*Nbus+Nbranch+1:2*Ngen+2*Nbus+2*Nbranch;
entries{7}  = 2*Ngen+2*Nbus+2*Nbranch+1:2*Ngen+2*Nbus+3*Nbranch;  % Q_u(to_bus)
entries{8}  = 2*Ngen+2*Nbus+3*Nbranch+1:2*Ngen+2*Nbus+4*Nbranch;  
entries{9}  = 2*Ngen+2*Nbus+4*Nbranch+1:2*Ngen+2*Nbus+5*Nbranch;  % l_u
entries{10} = 2*Ngen+2*Nbus+5*Nbranch+1:2*Ngen+2*Nbus+6*Nbranch;

x0 = vertcat(pg_0(1),qg_0(1),pg_0(2:end),qg_0(2:end),V_0,V_0,P_0,P_0,Q_0,Q_0,l_0,l_0);
lbx = vertcat(pg_min(1),qg_min(1),pg_min(2:end),qg_min(2:end),V_min,V_min,-inf(4*Nbranch,1),zeros(2*Nbranch,1));
ubx = vertcat(pg_max(1),qg_max(1),pg_min(2:end),qg_min(2:end),V_max,V_max,inf(6*Nbranch,1));
%% constraints
% envelope on the variables

eq_1 = @(x)create_VPQ(x(entries{2}(1:Ngen-1)),x(entries{2}(Ngen:2*Ngen-2)),...
    x(entries{3}),x(entries{4}),x(entries{5}),x(entries{6}),x(entries{7}),x(entries{8}),x(entries{9}),x(entries{10}),...
    Cg,C,DR,DX_p,DX_m,Pd,Qd,V_0,Mp,Mq,H_p,H_m);

% envelope on the current 
ineq_current = @(x)create_current(x(entries{3}),x(entries{4}),x(entries{5}),x(entries{6}),x(entries{7}),x(entries{8}),x(entries{9}),x(entries{10}),...
    V_0,P_0,Q_0,l_0,J_p,J_m,Hes,from_bus);


% bounds for pcc
% ineq_pcc = @(x)create_pcc(x(entries{1}(1)),x(entries{1}(2)),x(entries{5}),x(entries{6}),x(entries{7}),x(entries{8}),x(entries{9}),x(entries{10}),R,X,id_sline);

g = @(x)vertcat(eq_1(x),ineq_current(x));

g_test = g(x0);
num_g = numel(g_test);

lbg = vertcat(zeros(6*Nbranch,1),-inf(num_g-6*Nbranch,1));
ubg = vertcat(zeros(6*Nbranch,1),zeros(num_g-6*Nbranch,1));
%% solver options
import casadi.*
% tolerance
tol        = 1e-6;
options.ipopt.tol             = tol;
options.ipopt.constr_viol_tol = tol;
options.ipopt.compl_inf_tol   = tol;
options.ipopt.acceptable_tol  = tol;
options.ipopt.acceptable_constr_viol_tol = tol;
options.ipopt.print_level = 5;
% options.ipopt.grad_f = fgrad;
options.print_time        = 5;
options.ipopt.max_iter    = 200;

Nx         =  numel(x0);
x_SX       =   SX.sym('x',Nx,1);
constraint = g(x_SX);
%% objective and solution
sum_p = @(x) x(entries{5}(1));
obj = sum_p(x_SX);
nlp = struct('x',x_SX,'f',-obj,'g',constraint);
S   = nlpsol('solver','ipopt', nlp,options);
sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
        'lbx', lbx, 'ubx', ubx);
xopt= full(sol.x);
%% plot