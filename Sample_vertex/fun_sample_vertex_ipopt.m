function [vert] = fun_sample_vertex_ipopt(mpc)
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
id_gen_nslack = find(id_gen ~= id_slack);
id_gen_slack = find(id_gen == id_slack);

baseMVA     = mpc.baseMVA;              % baseMVA
cost_param  = mpc.gencost(:,5:end);     % objective coefficients
Umax        = mpc.bus(:,VMAX).^2;                                            
Umin        = mpc.bus(:,VMIN).^2;
Phimax      = mpc.branch(:,ANGMAX)/180*pi;
Phimin      = mpc.branch(:,ANGMIN)/180*pi;
Pgmin       = mpc.gen(:,PMIN)/baseMVA; Pgmin(id_gen_slack) = -inf;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; Qgmin(id_gen_slack) = -inf;
Pgmax       = mpc.gen(:,PMAX)/baseMVA; Pgmax(id_gen_slack) = inf;
Qgmax       = mpc.gen(:,QMAX)/baseMVA; Qgmax(id_gen_slack) = inf;
Fmax        = mpc.branch(:,RATE_A)/baseMVA;


Pd          = mpc.bus(:,PD)/baseMVA;  
% change the demand to be the supplier
Pd([19 25 28]) = -1*Pd([19 25 28]);
Qd          = mpc.bus(:,QD)/baseMVA;
Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);

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
D              = (eye(Nbranch) - A)\eye(size(A));
DR             = (eye(Nbranch) - A)\(A*R);
DX             = (eye(Nbranch) - A)\(A*X);
DX_p = max(DX,zeros(size(DX))).*(DX>0);
DX_m = min(DX,zeros(size(DX))).*(DX<0);

Mp = D'*R*D;
Mq = D'*X*D;
%% data of the load
% load("Data_load.mat")
% load_error = allData; % in percentage
% load_error = load_error.*Pd(2:33);
% e_l = min(load_error,[],2);
% e_u = max(load_error,[],2);
% mu_load = mean(load_error,2);
% sig_load = var(load_error, 0, 2);
% sig_load = 0.01*Pd(2:end);
% sig_load = sig_load*2;

% cov_load= cov(load_error');

% Pdnd = sqrt(sig_load);
% Pdnd = 0.1*Pd(2:end);
Pdnd = zeros(32,1);

epsl = 0.32;
% epsl = 0.05;
% epsl = 0.003;
Keps = sqrt((1-epsl)/epsl);
%% set the variables 
import casadi.*
U          = SX.sym('U',Nbus,1);
Pij        = SX.sym('Pij',Nbranch,1);
Qij        = SX.sym('Qij',Nbranch,1);
Pg         = SX.sym('Pg',Ngen,1);
Qg         = SX.sym('Qg',Ngen,1);
Pnd        = SX.sym('Pnd',Ngen-1,1);
Und        = SX.sym('Und',Nbranch,1);
x          = vertcat(U, Pij, Qij, Pg, Qg, Pnd, Und);
%% lower & upper bounds
Pgmin_v = Pgmin;
Pgmin_v(2:end) = -inf;
Pgmax_v = Pgmax;
Pgmax_v(2:end) = inf;
lbx         = [Umin(1);-inf(Nbranch,1);-inf(2*Nbranch,1);Pgmin_v;Qgmin;-inf(Ngen-1+Nbranch,1)];       
ubx         = [Umax(1);inf(Nbranch,1); inf(2*Nbranch,1); Pgmax_v;Qgmax;inf(Ngen-1+Nbranch,1)];
%% initial state x0

U0          = mpc.bus(:,VM).^2;
Pg0         = mpc.gen(:,PG)/baseMVA;
Qg0         = mpc.gen(:,QG)/baseMVA;
% Pij0        = max(mpc.branch(:,14), mpc.branch(:,16));
% Qij0        = max(mpc.branch(:,15), mpc.branch(:,17));
Pij0        = zeros(Nbranch,1);
Qij0        = zeros(Nbranch,1);
x0      = vertcat(U0, Pij0, Qij0, Pg0, Qg0, zeros(Ngen+Nbranch-1,1));
%% LinDistFlow

% voltage constraints
volt_eq    = C*U -2*(branch_r .* Pij + branch_x .* Qij);
% volt_eq = U(2:end) - (2*Mq*(Cg(2:end,2:end)*Qg(2:end)-Qd(2:end))+2*Mp*(Cg(2:end,2:end)*Pg(2:end)-Pd(2:end)))-1;
% power balance
pf_p_eq    = Cg*Pg - Pd - C'*Pij;
pf_q_eq    = Cg*Qg - Qd - C'*Qij;
% uncertainty for Pgen, U
un_Gnd      = sum(Pnd) - sum(Pdnd);
un_Und      = Und - (2*Mp*(Cg(2:end,2:end)*Pnd - Pdnd));
% bounds for chance constraints
CC_ineq    = create_PCE(Pg,Pnd,U,Und,Keps,Pgmax,Pgmin,Umax,Umin);
%% Problem formulation
gfun = vertcat(volt_eq,pf_p_eq,pf_q_eq,un_Gnd,un_Und,CC_ineq);
lbg = [zeros(2*Nbranch+2*Nbus+1,1);-inf(4*(Ngen-1)+4*Nbranch,1);];
ubg = [zeros(2*Nbranch+2*Nbus+1,1);zeros(4*(Ngen-1)+4*Nbranch,1);];

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
% options.ipopt.max_iter    = 200;

constraint = gfun;
%% sampling points

c1c2 = [1 0; -1 0; 0 1; 0 -1];
vert = zeros(4,2);
for i = 1:4
    f_samp   = c1c2(i,:)*[obj_p;obj_q]; % p_{k,l}, q_{k,l} of PCC
    objective = f_samp;
    nlp = struct('x',x,'f',objective,'g',gfun);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    xopt= full(sol.x);
    obj_p_opt = xopt(Nbus+2*Nbranch+id_gen_slack);
    obj_q_opt = xopt(Nbus+2*Nbranch+Ngen+id_gen_slack);
    vert(i,:) = [obj_p_opt, obj_q_opt];
end
vert = unique_coor(vert, 1e-4);
% order the vertices
vert = order_points(vert);
% find the normal vectors
Normals = Norm_vector(vert);

% refine the vertices
H = 100;
vert_c = sum(vert,1)/size(vert,1);

while any(H >= 1e-4)
    H = zeros(size(Normals,1),1);
    for i = 1:size(Normals,1)
        f_ext = Normals(i,:)*([obj_p;obj_q]-vert_c');
        objective = f_ext;
        nlp = struct('x',x,'f',-objective,'g',gfun);
        S   = nlpsol('solver','ipopt', nlp,options);
        sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                'lbx', lbx, 'ubx', ubx);
        xopt= full(sol.x);
        obj_p_opt = xopt(Nbus+2*Nbranch+id_gen_slack);
        obj_q_opt = xopt(Nbus+2*Nbranch+Ngen+id_gen_slack);
        
        H(i) = Normals(i,:)*[obj_p_opt;obj_q_opt]-1;
        if  H(i) >= 1e-4 % no movement outwards
            vert(end+1,:) = [obj_p_opt,obj_q_opt];
        end
    end
    % update the central point and normal vector
    vert_c = sum(vert,1)/size(vert,1);
    vert = unique_coor(vert, 1e-4);
    vert = order_points(vert);
    Normals_new = Norm_vector(vert);
    Normals = unique_normal(Normals, Normals_new,1e-4);
end
toc
%% plot the results
vert(end+1,:) = vert(1,:);
plot(vert(:,1), vert(:,2), '-bo');
hold all;
% plot(sortedmu(:,1), sortedmu(:,2), '-ro');
xlabel('P/p.u.');   % X 轴标签
ylabel('Q/p.u.');   % Y 轴标签
title('Feasible Region of Slack Bus(LinDistFlow)'); % 图像标题
grid on;            % 显示网格


end