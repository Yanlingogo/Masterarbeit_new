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

mpc = ext2int(loadcase('case33_modified'));
mpc = ext2int(mpc);

%% parameters 
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

% B              = (Cf + Ct)';
% A              = [zeros(Nbranch,1) eye(Nbranch)]*B - eye(Nbranch);
% D              = (eye(Nbranch) - A)\eye(size(A));
% DR             = (eye(Nbranch) - A)\(A*R);
% DX             = (eye(Nbranch) - A)\(A*X);
% DX_p = max(DX,zeros(size(DX))).*(DX>0);
% DX_m = min(DX,zeros(size(DX))).*(DX<0);
% 
% Mp = D'*R*D;
% Mq = D'*X*D;
%% set the variables 
import casadi.*
U          = SX.sym('U',Nbus,1);
Pij        = SX.sym('Pij',Nbranch,1);
Qij        = SX.sym('Qij',Nbranch,1);
Pg         = SX.sym('Pg',Ngen,1);
Qg         = SX.sym('Qg',Ngen,1);
x          = vertcat(U, Pij, Qij, Pg, Qg);
%% lower & upper bounds
lbx         = [Umin;-inf(2*Nbranch,1);Pgmin;Qgmin];       
ubx         = [Umax; inf(2*Nbranch,1);Pgmax;Qgmax];
%% initial state x0

U0          = mpc.bus(:,VM).^2;
Pg0         = mpc.gen(:,PG)/baseMVA;
Qg0         = mpc.gen(:,QG)/baseMVA;
% Pij0        = max(mpc.branch(:,14), mpc.branch(:,16));
% Qij0        = max(mpc.branch(:,15), mpc.branch(:,17));
Pij0        = zeros(Nbranch,1);
Qij0        = zeros(Nbranch,1);
x0      = vertcat(U0, Pij0, Qij0, Pg0, Qg0);
%% LinDistFlow

% voltage constraints
volt_eq    = C*U -2*(branch_r .* Pij + branch_x .* Qij);
% volt_eq = U(2:end) - (2*Mq*(Cg(2:end,2:end)*Qg(2:end)-Qd(2:end))+2*Mp*(Cg(2:end,2:end)*Pg(2:end)-Pd(2:end)))-1;
% power balance
pf_p_eq    = Cg*Pg - Pd - C'*Pij;
pf_q_eq    = Cg*Qg - Qd - C'*Qij;
% ref bus
ref_eq     = U(id_slack) - mpc.bus(id_slack, VM).^2;
%% Problem formulation
gfun = vertcat(pf_p_eq, pf_q_eq, volt_eq, ref_eq);
lbg = zeros(2*Nbus+Nbranch+1,1);
ubg = zeros(2*Nbus+Nbranch+1,1);
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
options.ipopt.max_iter    = 100;

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
resolution = 30;
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
gridding_q= zeros(20,2); % store the solution under gridding
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
hold all;
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
toc;
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

% mpc = runpf(mpc); %pg = 0, qg = 0
% U_0 = mpc.bus(:,8).^2;
% p_pcc0 = mpc.gen(id_gen_slack,PG)/mpc.baseMVA;
% q_pcc0 = mpc.gen(id_gen_slack,QG)/mpc.baseMVA;
% s_0 = [p_pcc0,q_pcc0];
% P = mpc.branch(:,14)/mpc.baseMVA;
% Q = mpc.branch(:,15)/mpc.baseMVA;
% L = (P.^2+Q.^2)./U_0(from_bus);
% z_0 = [U_0; P; Q; p_pcc0; q_pcc0; L];
% 
% e_st   = sparse(id_slack,1,1,Nbus,1);
% P_p = Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
% Q_p = Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
% H =blkdiag(P_p,Q_p); 
% 
% D    = full([1, zeros(1,Nbus*2+Nbranch);
%         C, -2*diag(branch_r),-2*diag(branch_x), zeros(Nbranch,2);
%         zeros(Nbus,Nbus),C',zeros(Nbus,Nbranch),-e_st,zeros(Nbus,1);
%         zeros(Nbus,Nbus+Nbranch),C',zeros(Nbus,1),-e_st]);
% RR  = sparse(diag(branch_r));
% XX  = sparse(diag(branch_x));
% PP  = sparse(diag(P));
% QQ  = sparse(diag(Q));
% Ga  = sparse(diag(L.^(-1)));
% 
% E   = [sparse(1,Nbranch);RR^2+XX^2;Ct' *RR;Ct'*XX];
% F   = [Cf, - 2*PP*Ga, -2*QQ*Ga,sparse(Nbranch,2)];
% G   = (PP^2+QQ^2) * Ga^2;
% 
% jac_z = [D E;F G]; % jacobian 
% jac_u = [zeros(Nbus,2*(Ngen-1)); 
%          -Cg_ns zeros(Nbus,Ngen-1);
%          zeros(Nbus,Ngen-1) -Cg_ns;
%          zeros(Nbranch,2*(Ngen-1))];
% 
% delta_s = sortedPoints - s_0;
% z_comp = z_0 + jac_z\jac_u*H*delta_s';
% U_comp = z_comp(1:Nbus,:);
% Pij_comp = z_comp(Nbus+1:Nbus+Nbranch,:);
% Qij_comp = z_comp(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
% L_comp = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:);
% PQ_loss = ([branch_r';branch_x']*L_comp)';
% 
% comp_pcc = sortedPoints + PQ_loss;
% plot(comp_pcc(:,1),comp_pcc(:,2),'gx');
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
mpc = runpf(mpc);
U_0 = mpc.bus(:,8).^2;
p_pcc0 = mpc.gen(id_gen_slack,PG)/mpc.baseMVA;
q_pcc0 = mpc.gen(id_gen_slack,QG)/mpc.baseMVA;
s_0 = [p_pcc0,q_pcc0];
Pij_0 = mpc.branch(:,14)/mpc.baseMVA;
Qij_0 = mpc.branch(:,15)/mpc.baseMVA;
L_0 = (Pij_0.^2+Qij_0.^2)./U_0(from_bus);
theta_0 = [diag(Pij_0); diag(Qij_0); diag(U_0(2:end))];

e_st   = sparse(id_slack,1,1,Nbus,1);

A = [e_st' zeros(1,2*Nbus);
     C -2*R -2*X zeros(Nbranch,2);
     zeros(Nbus) -C' zeros(Nbus,Nbranch) e_st zeros(Nbus,1);
     zeros(Nbus,Nbus+Nbranch) -C' zeros(Nbus,1) e_st];
B = [zeros(Nbus,2*(Ngen-1));
     Cg_ns zeros(Nbus,Ngen-1);
     zeros(Nbus,Ngen-1), Cg_ns];

b = [1;zeros(Nbranch,1);Pd;Qd];

P_p = Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
Q_p = Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
H =blkdiag(P_p,Q_p); % participation factor

u = (-sortedPoints + [sum(Pd) sum(Qd)])*H';
x = A\(b - B*u');

J = [diag(2*Pij_0./U_0(from_bus)); diag(2*Qij_0./U_0(from_bus)); diag(- (Pij_0.^2+Qij_0.^2)./U_0(from_bus).^2)];
Hes = [diag(2./U_0(from_bus)) zeros(Nbus-1) diag((-2*Pij_0)./(U_0(from_bus).^2));
       zeros(Nbus-1) diag(2./U_0(from_bus)) diag((-2*Qij_0)./(U_0(from_bus).^2));
       diag((-2*Pij_0)./(U_0(from_bus).^2)) diag((-2*Qij_0)./(U_0(from_bus).^2)) diag(2*(Pij_0.^2+Qij_0.^2)./(U_0(from_bus).^3))];

l_comp = zeros(Nbranch,size(sortedPoints,1));
for i = 1: size(sortedPoints,1)
    x_index = x(:,i);
    U_comp = x_index(1:Nbus,:);
    Pij_comp = x_index(Nbus+1:Nbus+Nbranch,:);
    Qij_comp = x_index(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
    theta_comp = [diag(Pij_comp);diag(Qij_comp);diag(U_comp(2:end,:))];
    
    l_comp(:,i) = diag(diag(L_0) + J'*(theta_comp-theta_0) + 0.5*(theta_comp-theta_0)'*Hes*(theta_comp-theta_0));
    % l_comp(:,i) = diag(diag(L_0) + J'*(theta_comp-theta_0));
end
PQ_loss = [(branch_r'*l_comp)' (branch_x'*l_comp)'];

comp_pcc = sortedPoints + PQ_loss;
plot(comp_pcc(:,1),comp_pcc(:,2),'gx');
