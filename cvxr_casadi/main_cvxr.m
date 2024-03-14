clear;
clc;
mpc = loadcase('case9.m');
mpc = runopf(mpc);
%% index
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

nb = size(mpc.bus,1);
Ngen = size(mpc.gen,1);
nl = size(mpc.branch,1);
% parameters 
bus = mpc.bus;
gen = mpc.gen;
branch = mpc.branch;

bus_type = bus(:,BUS_TYPE);
idx_slack = find(bus_type == 3);
idx_nslack = find(bus_type ~= 3);
Nnslack = numel(idx_nslack);
idx_sgen = find(gen(:,GEN_BUS)==idx_slack);

idx_fr = branch(:,F_BUS);
idx_to = branch(:,T_BUS);
id_gen = gen(:,GEN_BUS);
%idx_load mentioned before
idx_pq = find(bus_type == 1);
onoff_pq = zeros(nb,1); onoff_pq = ones(size(onoff_pq));
idx_pvs = find(bus_type ~=1);

%gencost = mpc.gencost(:,COST:end);

E_fr = sparse(idx_fr, 1:nl, 1, nb, nl);
E_to = sparse(idx_to, 1:nl, 1, nb, nl);
E = E_fr - E_to;
Cg = sparse(id_gen,1:Ngen,1,nb,Ngen);
Cl = sparse(idx_nslack,1:Nnslack,1,nb,Nnslack);

v0 = bus(:,VM);
theta0 = deg2rad(bus(:,VA));
Phi0 = E'*theta0;
% alpha0 = gen(:,end);
% delta0 = mpc.delta;
pg0 = (gen(:,GEN_STATUS).*gen(:,PG))/mpc.baseMVA;
qg0 = (gen(:,GEN_STATUS).*gen(:,QG))/mpc.baseMVA;
pl0 = bus(idx_nslack,PD)/mpc.baseMVA;
ql0 = bus(idx_nslack,QD)/mpc.baseMVA;

% limitation of the network
vmax = bus(:,VMAX);
vmin = bus(:,VMIN);
Phi_max = deg2rad(branch(:,ANGMAX));
Phi_min = deg2rad(branch(:,ANGMIN));
pg_max = (gen(:,GEN_STATUS).*gen(:,PMAX))/mpc.baseMVA;
pg_min = (gen(:,GEN_STATUS).*gen(:,PMIN))/mpc.baseMVA;
qg_max = (gen(:,GEN_STATUS).*gen(:,QMAX))/mpc.baseMVA;
qg_min = (gen(:,GEN_STATUS).*gen(:,QMIN))/mpc.baseMVA;

% s_line_max = branch(:,RATE_A)/mpc.baseMVA;

% Jacobian 
idx_sq = 1:nb;
% idx_sq(idx_sq == 2) =[];
idx_sq(idx_sq == 3) =[]; % Gen3 for reactive power balance

J1 = -diag(v0(idx_fr).*v0(idx_to).*sin(E'*theta0))*E'; 
J2 = diag(v0(idx_to).*cos(E'*theta0))*E_fr'+diag(v0(idx_fr).*cos(E'*theta0))*E_to';
J3 = diag(v0(idx_fr).*v0(idx_to).*cos(E'*theta0))*E';
J4 = diag(v0(idx_to).*sin(E'*theta0))*E_fr'+diag(v0(idx_fr).*sin(E'*theta0))*E_to';
D = diag(v0);

 
% J_psi0 = [zeros(2*nb,2*nb-1);
%           J1, J2(:,idx_sq);
%           J3, J4(:,idx_sq);
%           zeros(nb-1,nb), 2*D];
J_psi0 = [zeros(2*nb,2*nb);
          J1, J2;
          J3, J4;
          zeros(nb), 2*D];
% original value of residues and basis function
g_pinj0 = Cg*pg0 - Cl*pl0;
g_qinj0 = Cg*qg0 - Cl*ql0;

g_vvcos0 = v0(idx_fr).*v0(idx_to).*cos(Phi0) - onoff_pq(idx_fr).*v0(idx_fr).*v0(idx_to).*cos(Phi0)...
    - v0(idx_fr).*onoff_pq(idx_to).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi0;
g_vvsin0 = v0(idx_fr).*v0(idx_to).*sin(Phi0) - onoff_pq(idx_fr).*v0(idx_fr).*v0(idx_to).*sin(Phi0)...
    - v0(idx_fr).*onoff_pq(idx_to).*v0(idx_to).*sin(Phi0) - v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi0;
g_vv0 = v0.^2-2*onoff_pq.*v0.*v0;
g0 = [g_pinj0;g_qinj0;g_vvcos0;g_vvsin0;g_vv0];

Psi_vvcos0 = v0(idx_fr).*v0(idx_to).*cos(Phi0);
Psi_vvsin0 = v0(idx_fr).*v0(idx_to).*sin(Phi0);
Psi_vv0 = v0.^2;

[~,~,~,M,M_line] = makeYbus_cvxr(mpc,false); 
M_eq = M([1:nb nb+idx_sq],:);
M_ineq = M(nb+3,:);

J = M*J_psi0;
% for full J 
rowsToDelete = 12; 
colsToDelete = 1; 

rowsToKeep = true(1, size(J, 1)); 
rowsToKeep(rowsToDelete) = false; 

colsToKeep = true(1, size(J, 2)); 
colsToKeep(colsToDelete) = false; 

J = J(rowsToKeep, colsToKeep); 

J_invM = J\M_eq;

A = [E(idx_nslack,:)' zeros(nl,nb);
     zeros(nb,nb-1) eye(nb);
     -E(idx_nslack,:)' zeros(nl,nb);
     zeros(nb,nb-1) -eye(nb)]; 

K = -A*J_invM;

K_plus  = max(K,zeros(size(K))).*(K>0);
K_minus = min(K,zeros(size(K))).*(K<0);

M_ineq_p = max(M_ineq,zeros(size(M_ineq))).*(M_ineq>0);
M_ineq_m = max(M_ineq,zeros(size(M_ineq))).*(M_ineq<0);


%% optimization problem
entries{1}  = 1:nb;                         %v_u
entries{2}  = nb+1:2*nb;                    %v_l
entries{3}  = 2*nb+1:3*nb;                  %Delta_v_u
entries{4}  = 3*nb+1:4*nb;                  %Delta_v_l

entries{5}  = 4*nb+1:4*nb+nl;               %Phi_u
entries{6}  = 4*nb+nl+1:4*nb+2*nl;          %Phi_l
entries{7}  = 4*nb+2*nl+1:4*nb+3*nl;        %Delta_Phi_u
entries{8}  = 4*nb+3*nl+1:4*nb+4*nl;        %Delta_Phi_l

entries{9}  = 4*nb+4*nl+1:5*nb+4*nl;        %g_u_pinj
entries{10} = 5*nb+4*nl+1:6*nb+4*nl;        %g_l_pinj
entries{11} = 6*nb+4*nl+1:7*nb+4*nl;        %g_u_qinj
entries{12} = 7*nb+4*nl+1:8*nb+4*nl;        %g_l_qinj

entries{13} = 8*nb+4*nl+1:8*nb+5*nl;        %g_u_vvcos
entries{14} = 8*nb+5*nl+1:8*nb+6*nl;        %g_l_vvcos
entries{15} = 8*nb+6*nl+1:8*nb+7*nl;        %g_u_vvsin
entries{16} = 8*nb+7*nl+1:8*nb+8*nl;        %g_l_vvsin

entries{17} = 8*nb+8*nl+1:9*nb+8*nl;        %g_u_vv
entries{18} = 9*nb+8*nl+1:10*nb+8*nl;       %g_l_vv

entries{19} = 10*nb+8*nl+1:10*nb+16*nl;     %D_vv_u_uu...
entries{20} = 10*nb+16*nl+1:10*nb+24*nl;    %D_cos_u_u...

% entries{21} = 10*nb+24*nl+1:10*nb+25*nl;        %Psi_uvvc
% entries{22} = 10*nb+25*nl+1:10*nb+26*nl;        %Psi_lvvc
% entries{23} = 10*nb+26*nl+1:10*nb+27*nl;        %Psi_uvvs
% entries{24} = 10*nb+27*nl+1:10*nb+28*nl;        %Psi_lvvs

entries{25} = 10*nb+28*nl+1:11*nb+28*nl;        %Psi_uvv
entries{26} = 11*nb+28*nl+1:12*nb+28*nl;        %Psi_lvv

entries{27} = 12*nb+28*nl+1:12*nb+28*nl+2*Ngen;      %pg_opt, qg_opt

%% lower & upper bounds for the variables 
lbx = [vmin;vmin;vmin-v0;vmin-v0;Phi_min;Phi_min;Phi_min-Phi0;Phi_min-Phi0;-inf(8*nb+24*nl+2*Ngen,1)];
ubx = [vmax;vmax;vmax-v0;vmax-v0;Phi_max;Phi_max;Phi_max-Phi0;Phi_max-Phi0; inf(8*nb+24*nl+2*Ngen,1)];
%% initial state x0

x0 = vertcat(v0,v0,zeros(2*nb,1),Phi0,Phi0,zeros(8*nb+26*nl,1),pg0,qg0);
%% constraints

cons0 = @(x) create_base(x(entries{1}),x(entries{2}),x(entries{3}),x(entries{4}),x(entries{5}),x(entries{6}),x(entries{7}),x(entries{8}),v0,Phi0,idx_slack);

cons1 = @(x) create_vv(x(entries{19}(1:nl)),x(entries{19}(nl+1:2*nl)),x(entries{19}(2*nl+1:3*nl)),x(entries{19}(3*nl+1:4*nl)),...
        x(entries{19}(4*nl+1:5*nl)),x(entries{19}(5*nl+1:6*nl)),x(entries{19}(6*nl+1:7*nl)),x(entries{19}(7*nl+1:8*nl)),...
        x(entries{3}),x(entries{4}),v0,idx_fr,idx_to);

cons2 = @(x) create_cossin(x(entries{20}(1:nl)),x(entries{20}(nl+1:2*nl)),x(entries{20}(2*nl+1:3*nl)),x(entries{20}(3*nl+1:4*nl)),...
        x(entries{20}(4*nl+1:5*nl)),x(entries{20}(5*nl+1:6*nl)),x(entries{20}(6*nl+1:7*nl)),x(entries{20}(7*nl+1:8*nl)),...
        x(entries{7}),x(entries{8}),Phi0);

cons3 = @(x) create_g_inj(x(entries{9}),x(entries{10}),x(entries{11}),x(entries{12}),x(entries{27}(1:Ngen)),x(entries{27}(Ngen+1:2*Ngen)),Cg,Cl,pl0,ql0,pg_max,pg_min,qg_max,qg_min);
% cons3 = @(x) create_g_inj_slack(x(entries{9}),x(entries{10}),x(entries{11}),x(entries{12}),x(entries{27}(1:Ngen)),x(entries{27}(Ngen+1:2*Ngen)),x(entries{28}(1:Ngen)),x(entries{28}(Ngen+1:2*Ngen)),Cg,Cl,pl0,ql0,pg_max,pg_min,qg_max,qg_min);

cons4 = @(x) create_gvvcs(x(entries{13}),x(entries{14}),x(entries{15}),x(entries{16}),...
        x(entries{19}(1:nl)),x(entries{19}(nl+1:2*nl)),x(entries{19}(2*nl+1:3*nl)),x(entries{19}(3*nl+1:4*nl)),...
        x(entries{19}(4*nl+1:5*nl)),x(entries{19}(5*nl+1:6*nl)),x(entries{19}(6*nl+1:7*nl)),x(entries{19}(7*nl+1:8*nl)),...
        x(entries{20}(1:nl)),x(entries{20}(nl+1:2*nl)),x(entries{20}(2*nl+1:3*nl)),x(entries{20}(3*nl+1:4*nl)),...
        x(entries{20}(4*nl+1:5*nl)),x(entries{20}(5*nl+1:6*nl)),x(entries{20}(6*nl+1:7*nl)),x(entries{20}(7*nl+1:8*nl)),...
        x(entries{1}),x(entries{2}),x(entries{5}),x(entries{6}),Psi_vvcos0,Psi_vvsin0,Phi0,v0,idx_fr,idx_to,onoff_pq);

cons5 = @(x) create_gvv(x(entries{17}),x(entries{18}),...
        x(entries{1}),x(entries{2}),v0,onoff_pq,Psi_vv0);

% cons6 = @(x) create_Psi(x(entries{21}),x(entries{22}),x(entries{23}),x(entries{24}),x(entries{25}),x(entries{26}),...
%     x(entries{19}(1:nl)),x(entries{19}(4*nl+1:5*nl)),x(entries{19}(5*nl+1:6*nl)),...
%     x(entries{20}(1:nl)),x(entries{20}(nl+1:2*nl)),x(entries{20}(2*nl+1:3*nl)),x(entries{20}(3*nl+1:4*nl)),...cos set
%     x(entries{20}(4*nl+1:5*nl)),x(entries{20}(7*nl+1:8*nl)),...sin set
%     x(entries{1}),x(entries{2}),x(entries{3}),x(entries{4}),... v_u/l set
%     Psi_vvcos0,Psi_vvsin0,Psi_vv0,Phi0,v0,idx_fr,idx_to);

cons7 = @(x) create_Kg(x(entries{9}),x(entries{10}),x(entries{11}),x(entries{12}),x(entries{13}),x(entries{14}),x(entries{15}),x(entries{16}),x(entries{17}),x(entries{18}),x(entries{1}),x(entries{2}),x(entries{5}),x(entries{6}),K_plus,K_minus);

% cons8 = @(x) create_Mineq(x(entries{21}),x(entries{22}),x(entries{23}),x(entries{24}),x(entries{25}),x(entries{26}),Cg,Cl,qg_max,qg_min,ql0,M_ineq_p,M_ineq_m);
g = @(x)vertcat(cons0(x),cons1(x),cons2(x),cons3(x),cons4(x),cons5(x),cons7(x));

g_test = g(x0);
num_g = numel(g_test);

lbg = vertcat(zeros(2+2*nb+2*nl,1),-inf(num_g-(2+2*nb+2*nl),1));
ubg = vertcat(zeros(2+2*nb+2*nl,1),zeros(num_g-(2+2*nb+2*nl),1));

obj_p = @(x) x(entries{27}(1));
obj_q = @(x) x(entries{27}(Ngen+1));

%% solution
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
% options.ipopt.max_iter    = 500;
options.ipopt.warm_start_init_point ='yes';

Nx         = numel(x0);
x_SX       = SX.sym('x',Nx,1);
constraint = g(x_SX);

%% sampling 
c_set = [-1 -1;-1 0;-1 1;0 -1;0 1;1 -1;1 0;1 1];
flag = [];
for i = 1:8
    c1 = c_set(i,1);
    c2 = c_set(i,2);
    f_samp   = @(x)c1*obj_p(x)+c2*obj_q(x); % p_{k,l}, q_{k,l} of PCC
    objective = f_samp(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    status = S.stats();
    if strcmp(status.return_status, 'Solve_Succeeded')
        flag = [flag;1];
    end
    xopt= full(sol.x);
    obj_values(i)= full(sol.f);
    obj_p_opt = obj_p(xopt);
    obj_q_opt = obj_q(xopt);
    Points(i,:) = [obj_p_opt,obj_q_opt];
    i = i+1;
end
D_angle = deg2rad(10);
Points_c = mean(Points);
for i = 0:D_angle:2*pi

    c1 = cos(i);
    c2 = sin(i); 

    angle_cons = @(x) [c2 -c1]*([obj_p(x) obj_q(x)]-Points_c)';
    g_ext = @(x) vertcat(g(x),angle_cons(x));
    cons_ext = g_ext(x_SX);
    lbg_ext = vertcat(lbg,0);
    ubg_ext = vertcat(ubg,0);

    Obj_ext = @(x) [c1 c2]*([obj_p(x) obj_q(x)]-Points_c)';
    objective = Obj_ext(x_SX);

    nlp = struct('x',x_SX,'f',-objective,'g',cons_ext);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg_ext,'ubg', ubg_ext,...
            'lbx', lbx, 'ubx', ubx);
    status = S.stats();
    if strcmp(status.return_status, 'Solve_Succeeded')
        flag = [flag;1];
    end
    xopt= full(sol.x);
    obj_p_opt = obj_p(xopt);
    obj_q_opt = obj_q(xopt);
    Points = [Points;[obj_p_opt,obj_q_opt]];
end
Points = Orderpoints(Points);