clc
clear
close all
%%Index setting
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

mpc = loadcase('case33_modified.m');
mpc = runopf(mpc);
% mpc.gen(:,PMAX) = inf;
% mpc.gen(:,PMIN) = -inf;

% mpc,branch(:,RATE_A) = 
% add Generators in mpc
% Info_gen = [6	0	0	100	-100	1	1	1	100	10	0	0	0	0	0	0	0	0	0	0	0;
% 	        16	0	0	100	-100	1	1	1	100	10	0	0	0	0	0	0	0	0	0	0	0;];
% mpc.gen = vertcat(mpc.gen,Info_gen);

%% parameters 
bus = mpc.bus;
gen = mpc.gen;
branch =mpc.branch;

id_bus      = mpc.bus(:,BUS_I);
id_gen      = mpc.gen(:,GEN_BUS);
Nbus        = size(mpc.bus,1);
Ngen        = numel(id_gen);
Nbranch     = size(mpc.branch,1);

% idx
id_slack =  find(mpc.bus(:,BUS_TYPE) == REF);
Nslack = numel(id_slack);
id_gen_nslack = find(id_gen ~= id_slack);
id_gen_slack = find(id_gen == id_slack);

baseMVA     = mpc.baseMVA;              % baseMVA

Umax        = mpc.bus(:,VMAX).^2;                                            
Umin        = mpc.bus(:,VMIN).^2;

Pgmin       = mpc.gen(:,PMIN)/baseMVA; %Pgmin(id_gen_slack) = -inf;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; %Qgmin(id_gen_slack) = -inf;
Pgmax       = mpc.gen(:,PMAX)/baseMVA; %Pgmax(id_gen_slack) = inf;
Qgmax       = mpc.gen(:,QMAX)/baseMVA; %Qgmax(id_gen_slack) = inf;
Fmax        = mpc.branch(:,RATE_A)/baseMVA;

Pd          = mpc.bus(:,PD)/baseMVA;   
Qd          = mpc.bus(:,QD)/baseMVA;
Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);

% branch info
branch_r   = mpc.branch(:,BR_R);
branch_x   = mpc.branch(:,BR_X);
R = diag(branch_r);
X = diag(branch_x);

% 
from_bus       = mpc.branch(:, F_BUS);                         
to_bus         = mpc.branch(:, T_BUS);  
Cf             = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
Ct             = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
C              = Cf - Ct;
% Jacobian and Hessian
U_0 = bus(:,VM).^2;
pg0 = gen(:,PG)/mpc.baseMVA;
qg0 = gen(:,QG)/mpc.baseMVA;
Pij0 = branch(:,PF)/mpc.baseMVA;
Qij0 = branch(:,QF)/mpc.baseMVA;
l_0 = (Pij0.^2+Qij0.^2)./U_0(from_bus);

J = [diag(2*Pij0/U_0(from_bus)); diag(2*Qij0/U_0(from_bus)); diag(- (Pij0.^2+Qij0.^2)/U_0(from_bus).^2)];
J_p = max(J,zeros(size(J))).*(J>0);
J_m = min(J,zeros(size(J))).*(J<0);
Hes = [diag(2./U_0(from_bus)) zeros(Nbus-1) diag((-2*Pij0)./(U_0(from_bus).^2));
       zeros(Nbus-1) diag(2./U_0(from_bus)) diag((-2*Qij0)./(U_0(from_bus).^2));
       diag((-2*Pij0)./(U_0(from_bus).^2)) diag((-2*Qij0)./(U_0(from_bus).^2)) diag(2*(Pij0.^2+Qij0.^2)./(U_0(from_bus).^3))];
%% set the variables 
import casadi.*
U          = SX.sym('U',Nbus,1);
Pij        = SX.sym('Pij',Nbranch,1);
Qij        = SX.sym('Qij',Nbranch,1);
Pg         = SX.sym('Pg',Ngen,1);
Qg         = SX.sym('Qg',Ngen,1);
l          = SX.sym('l',Nbranch,1);
x          = vertcat(U, Pij, Qij, Pg, Qg,l);
%% lower & upper bounds
lbx         = [Umin;-inf(2*Nbranch,1);Pgmin;Qgmin;-inf(Nbranch,1)];       
ubx         = [Umax; inf(2*Nbranch,1);Pgmax;Qgmax;inf(Nbranch,1)];
%% initial state x0


x0      = vertcat(U_0, Pij0, Qij0, pg0, qg0, l_0);
%% DistFlow
% voltage constraints
volt_eq    = C*U -2*(R* Pij + X* Qij) + (R.^2+X.^2)*l;
% volt_eq    = C*U -2*(R* Pij + X* Qij);
% power balance
pf_p_eq    = Cg*Pg - Pd - C'*Pij - Ct'*R*l;
pf_q_eq    = Cg*Qg - Qd - C'*Qij - Ct'*X*l;

% pf_p_eq    = Cg*Pg - Pd - C'*Pij;
% pf_q_eq    = Cg*Qg - Qd - C'*Qij;

% current_eq = Cf*U.*l - (Pij.^2+Qij.^2);
theta = [Pij-Pij0;Qij-Qij0;U(from_bus)-U_0(from_bus)];
current_eq = l - (l_0+J'*theta+ 0.5*(theta'*Hes*theta));
%% Problem formulation
gfun = vertcat(pf_p_eq, pf_q_eq, volt_eq, current_eq);
lbg = zeros(2*Nbus+2*Nbranch,1);
ubg = zeros(2*Nbus+2*Nbranch,1);
% objective
obj_p = Pg(id_gen_slack);
obj_q = Qg(id_gen_slack);

%% solver options

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

constraint = gfun;
%% sampling points
i = 1;
Points = zeros(8,2);
obj_values = zeros(8,1);
for c1 = -1:1
    for c2 = -1:1
        if c1== 0 && c2 ==0
            continue;
        else
            f_samp   = c1*obj_p+c2*obj_q; % p_{k,l}, q_{k,l} of PCC
            objective = f_samp;
            nlp = struct('x',x,'f',objective,'g',gfun);
            S   = nlpsol('solver','ipopt', nlp,options);
            sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                    'lbx', lbx, 'ubx', ubx);
            xopt= full(sol.x);
            obj_values(i)= full(sol.f);
            obj_p_opt = xopt(Nbus+2*Nbranch+id_gen_slack);
            obj_q_opt = xopt(Nbus+2*Nbranch+Ngen+id_gen_slack);
            Points(i,:) = [obj_p_opt,obj_q_opt];
            i = i+1;
        end
    end
end
%% discretization at x-axis
resolution = 20;
s_pmin = min(Points(:,1));
s_pmax = max(Points(:,1));
s_p_step = linspace(s_pmin,s_pmax,resolution+2);
s_p_step = s_p_step(2:end-1);
gridding_p= zeros(20,2); % store the solution under gridding
eq_sp = obj_p;
g_ext = vertcat(gfun,eq_sp);
for c3 = [-1, 1]
    for i = 1:resolution
        lbg_ext = vertcat(lbg,s_p_step(i));
        ubg_ext = vertcat(ubg,s_p_step(i));
        f_grid = c3*obj_q; % q_{k,l} of PCC
        constraint_grid = g_ext;
        objective = f_grid;
        nlp = struct('x',x,'f',objective,'g',constraint_grid);
        S   = nlpsol('solver','ipopt', nlp,options);
        sol = S('x0', x0,'lbg', lbg_ext,'ubg', ubg_ext,...
                'lbx', lbx, 'ubx', ubx);
        xopt2= full(sol.x);
        obj_q_opt = xopt2(Nbus+2*Nbranch+Ngen+id_gen_slack);
        if c3 == -1
            gridding_p(i,1) = obj_q_opt;
        else
            gridding_p(i,2) = obj_q_opt;
        end
    end
end
%% discretization at y-axis
s_qmin = min(Points(:,2));
s_qmax = max(Points(:,2));
s_q_step = linspace(s_qmin,s_qmax,resolution+2);
s_q_step = s_q_step(2:end-1);
gridding_q= zeros(resolution,2); % store the solution under gridding
eq_sq = obj_q;
g_ext = vertcat(gfun,eq_sq);
for c3 = [-1, 1]
    for i = 1:resolution
        lbg_ext = vertcat(lbg,s_q_step(i));
        ubg_ext = vertcat(ubg,s_q_step(i));
        f_grid = c3*obj_p; % p_{k,l} of PCC
        constraint_grid = g_ext;
        objective = f_grid;
        nlp = struct('x',x,'f',objective,'g',constraint_grid);
        S   = nlpsol('solver','ipopt', nlp,options);
        sol = S('x0', x0,'lbg', lbg_ext,'ubg', ubg_ext,...
                'lbx', lbx, 'ubx', ubx);
        xopt2= full(sol.x);
        obj_p_opt = xopt2(Nbus+2*Nbranch+id_gen_slack);
        if c3 == -1
            gridding_q(i,1) = obj_p_opt;
        else
            gridding_q(i,2) = obj_p_opt;
        end
    end
end
%% output 画出的图像中，P/Q > 0 表示向DSO中注入功率
plot(Points(:,1),Points(:,2),'rx');
hold on;
Points_tot = [Points;[s_p_step';s_p_step'],[gridding_p(:,1);gridding_p(:,2)];...
    [gridding_q(:,1);gridding_q(:,2)],[s_q_step';s_q_step']];
centroid = mean(Points_tot, 1);
angles = atan2(Points_tot(:,2) - centroid(2), Points_tot(:,1) - centroid(1));
[~, order] = sort(angles);
sortedPoints = Points_tot(order, :);
plot([sortedPoints(:,1); sortedPoints(1,1)], [sortedPoints(:,2); sortedPoints(1,2)], '-bo');
xlabel('P/p.u.');   % X 轴标签
ylabel('Q/p.u.');   % Y 轴标签
title('Feasible Region of Slack Bus(LinDistFlow)'); % 图像标题
grid on;            % 显示网格

%% cost plot
toc