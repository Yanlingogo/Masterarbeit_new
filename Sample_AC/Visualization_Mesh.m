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
%mpc = loadcase('case9');
tic
mpc = loadcase('case33_mesh');
% mpc.gen(:,PMAX) = inf;
% mpc.gen(:,PMIN) = -inf;

% mpc,branch(:,RATE_A) = 
% add Generators in mpc
% Info_gen = [6	1.63	0.0654	10	-5	1.025	1	1	10	-3	0	0	0	0	0	0	0	0	0	0	0;
% 	        16	0.85	-0.1095	10	-5	1.025	1	1	5	-3	0	0	0	0	0	0	0	0	0	0	0;];
% mpc.gen = vertcat(mpc.gen,Info_gen);

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

entries_pf{1} = 1:Nbus;                        % vmag
entries_pf{2} = (Nbus+1):2*Nbus;               % vang
entries_pf{3} = (2*Nbus+1):(2*Nbus+Ngen);      % Pg
entries_pf{4} = (2*Nbus+Ngen+1):2*(Nbus+Ngen); % Qg 

baseMVA     = mpc.baseMVA;              % baseMVA
cost_param  = mpc.gencost(:,5:end);     % objective coefficients
vmax        = mpc.bus(:,VMAX);                                            
vmin        = mpc.bus(:,VMIN);
Phimax      = mpc.branch(:,ANGMAX)/180*pi;
Phimin      = mpc.branch(:,ANGMIN)/180*pi;
Pgmin       = mpc.gen(:,PMIN)/baseMVA; Pgmin(id_gen_slack) = -inf;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; Qgmin(id_gen_slack) = -inf;
Pgmax       = mpc.gen(:,PMAX)/baseMVA; Pgmax(id_gen_slack) = inf;
Qgmax       = mpc.gen(:,QMAX)/baseMVA; Qgmax(id_gen_slack) = inf;
Fmax        = mpc.branch(:,RATE_A)/baseMVA;


Pd          = mpc.bus(:,PD)/baseMVA;   
Qd          = mpc.bus(:,QD)/baseMVA;
Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);

%% lower & upper bounds
% original 
lbx         = [vmin;-pi*ones(Nbus,1);Pgmin;Qgmin];       
ubx         = [vmax; pi*ones(Nbus,1);Pgmax;Qgmax];

% acitve power bounds: Pgmax actictivated
% lbx         = [-inf(size(vmin));-pi*ones(Nbus,1);Pgmin;-inf(size(Qgmin))];       
% ubx         = [inf(size(vmax)); pi*ones(Nbus,1);Pgmax;inf(size(Qgmax))];

% bound 2
% vmin(:) = -inf;
% vmax(:) = inf;
% vmin(6) = 1.1;
% vmax(6) = 1.1;
% vmin(8) = 1.1;
% vmax(8) = 1.1;
% 
% Pgmin(:) = -inf;
% Pgmax(:) = inf;
% Pgmin(3) = 0.8;
% Pgmax(3) = 0.8;
% 
% lbx         = [vmin;-pi*ones(Nbus,1);Pgmin;-inf(size(Qgmin))];       
% ubx         = [vmax;pi*ones(Nbus,1);Pgmax;inf(size(Qgmax))];

% bound 3
% vmin(:) = -inf;
% vmax(:) = inf;
% vmin(7) = 1.1;
% vmax(7) = 1.1;
% vmin(8) = 1.1;
% vmax(8) = 1.1;
% 
% Pgmin(:) = -inf;
% Pgmax(:) = inf;
% 
% Qgmin(:) = -inf;
% Qgmax(:) = inf;
% Qgmin(2) = 1;
% Qgmax(2) = 1;
% 
% lbx         = [vmin;-pi*ones(Nbus,1);Pgmin;Qgmin];       
% ubx         = [vmax;pi*ones(Nbus,1);Pgmax;Qgmax];

% bound 4
vmin(:) = -inf;
vmax(:) = inf;
vmin(8) = -inf;
vmax(8) = 1.1;

Pgmin(:) = -inf;
Pgmax(:) = inf;

Qgmin(:) = -inf;
Qgmax(:) = inf;
Qgmin(3) = -inf;
Qgmax(3) = 1.5;

lbx         = [vmin;-pi*ones(Nbus,1);Pgmin;Qgmin];       
ubx         = [vmax; pi*ones(Nbus,1);Pgmax;Qgmax];

% bound 5 
% vmin(:) = -inf;
% vmax(:) = inf;
% 
% Pgmin(:) = -inf;
% Pgmax(:) = inf;
% 
% Qgmin(:) = -inf;
% Qgmax(:) = inf;
% Qgmin(3) = 1.5;
% Qgmax(3) = 1.5;
% 
% lbx         = [vmin;-pi*ones(Nbus,1);Pgmin;Qgmin];       
% ubx         = [vmax; pi*ones(Nbus,1);Pgmax;Qgmax];

% bound 6
% vmin(:) = -inf;
% vmax(:) = inf;
% 
% Pgmin(:) = -inf;
% Pgmax(:) = inf;
% 
% Qgmin(:) = -inf;
% Qgmax(:) = inf;
% 
% lbx         = [vmin;-pi*ones(Nbus,1);Pgmin;Qgmin];       
% ubx         = [vmax; pi*ones(Nbus,1);Pgmax;Qgmax];


% voltage lower bounds: b8
% vmin(:) = -inf;
% vmax(:) = inf;
% vmin(6) = 0.9;
% vmax(6) = 0.9;
% vmin(31) = 0.9;
% vmax(31) = 0.9;
% 
% Pgmin(:) = -inf;
% Pgmax(:) = inf;
% Pgmin(3) = -0.8;
% Pgmax(3) = -0.8;
% 
% lbx         = [vmin;-pi*ones(Nbus,1);Pgmin;-inf(size(Qgmin))];       
% ubx         = [vmax;pi*ones(Nbus,1);Pgmax;inf(size(Qgmax))];

% voltage lower bounds: b9
% vmin(:) = -inf;
% vmax(:) = inf;
% vmin(27) = 0.9;
% vmax(27) = 0.9;
% vmin(32) = 0.9;
% vmax(32) = 0.9;
% 
% Pgmin(:) = -inf;
% Pgmax(:) = inf;
% 
% Qgmin(:) = -inf;
% Qgmax(:) = inf;
% Qgmin(2) = -1;
% Qgmax(2) = -1;
% 
% lbx         = [vmin;-pi*ones(Nbus,1);Pgmin;Qgmin];       
% ubx         = [vmax;pi*ones(Nbus,1);Pgmax;inf(size(Qgmax))];
%% initial state x0
vang0       = mpc.bus(:,VA)/180*pi;
vmag0       = mpc.bus(:,VM);
U0          = vmag0.*cos(vang0);
W0          = vmag0.*sin(vang0);
Pg0         = mpc.gen(:,PG)/baseMVA;
Qg0         = mpc.gen(:,QG)/baseMVA;
x0      = vertcat(vmag0, vang0, Pg0, Qg0);
%% equality & inequality constraints - current balance constraints
% create Ybus Yf Yt
[Ybus, Yf, Yt] = makeYbus(mpc);
Gbus           = real(Ybus);
Bbus           = imag(Ybus);
Gf             = real(Yf);
Bf             = imag(Yf);
Gt             = real(Yt);
Bt             = imag(Yt);
from_bus       = mpc.branch(:, F_BUS);                           %% list of "from" buses
to_bus         = mpc.branch(:, T_BUS);  
Cf             = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
Ct             = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
C              = Cf - Ct;

% reference for slack bus
vmag_ref = entries_pf{1}(id_slack);
vang_ref = entries_pf{2}(id_slack);
v_ref    = vertcat(vang_ref,vmag_ref);
eq_ref = @(x)x(v_ref);

% power flow equation
eq_pf          = @(x)create_local_power_flow_equation_pol(x(entries_pf{1}),x(entries_pf{2}),...
    x(entries_pf{3}),x(entries_pf{4}),Gbus,Bbus,Pd,Qd,Cg);
test_pf        = eq_pf(x0);
Npf            = numel(eq_pf(x0));


% upper & lower bound for voltage angle
ineq_angle     = @(x) C*x(entries_pf{2});


%ineq_genpq = @(x)create_bound_PQ(x(entries_pf{3}),x(entries_pf{4}),alpha,id_gen_nslack,id_gen_slack);
% line Limit
if any(Fmax)
    ineq_line = @(x)create_local_branch_limit_rec(x(entries_pf{1}),...
        x(entries_pf{2}), Gf, Bf, Gt, Bt, Fmax, from_bus,to_bus);
    idx_limit = find(Fmax);
    Nlimit = numel(ineq_line(x0));
    %g = @(x)vertcat(ineq_voltage(x),ineq_genp(x),ineq_genq(x),ineq_line(x),ineq_genpq(x),eq_pf(x),eq_ref(x));
    g = @(x)vertcat(ineq_angle(x),ineq_line(x),eq_pf(x),eq_ref(x)); %without generator bound
else
    ineq_line = [];
    idx_limit = [];
    Nlimit   = 0;
    %g   = @(x)vertcat(ineq_voltage(x),ineq_genp(x),ineq_genq(x),ineq_genpq(x),eq_pf(x),eq_ref(x));
    g   = @(x)vertcat(ineq_angle(x),eq_pf(x),eq_ref(x));
end
%% linie power flow for PCC 
% find slack bus and connected bus
connected_buses = unique([mpc.branch(mpc.branch(:, 1) == id_slack, 2);...
    mpc.branch(mpc.branch(:, 2) == id_slack, 1)]);
id_cline = find(mpc.branch(:, 1) == id_slack | mpc.branch(:, 2) == id_slack); % connecting line
% feasible P Q at connecting line
% obj_p = @(x)create_coupling_branch_limit_p(x(entries_pf{1}),...
%      x(entries_pf{2}),id_slack,connected_buses,id_cline, Gf, Bf, Gt, Bt);
% obj_q = @(x)create_coupling_branch_limit_q(x(entries_pf{1}),...
%      x(entries_pf{2}),id_slack,connected_buses,id_cline, Gf, Bf, Gt, Bt);

obj_p = @(x) x(entries_pf{3}(id_gen_slack));
obj_q = @(x) x(entries_pf{4}(id_gen_slack));

% obj_p = @(x)create_obj_p(x(entries_pf{1}),x(entries_pf{2}),...
%     x(entries_pf{3}),x(entries_pf{4}),Gbus,Bbus,Pd,Qd,Cg,id_slack);
% obj_q = @(x)create_obj_q(x(entries_pf{1}),x(entries_pf{2}),...
%     x(entries_pf{3}),x(entries_pf{4}),Gbus,Bbus,Pd,Qd,Cg,id_slack);

% original 
lbg = vertcat(Phimin, -inf*ones(Nlimit,1), zeros(Npf+1,1),1);
ubg = vertcat(Phimax, zeros(Nlimit,1),zeros(Npf+1,1),1);
% selected line limit
% lbg = vertcat(Phimin, 0,-inf, zeros(Npf+1,1),1);
% ubg = vertcat(Phimax, zeros(Nlimit,1),zeros(Npf+1,1),1);

%% solver options
import casadi.*
% tolerance
tol        = 1e-6;
% options.ipopt.tol             = tol;
% options.ipopt.constr_viol_tol = tol;
% options.ipopt.compl_inf_tol   = tol;
% options.ipopt.acceptable_tol  = tol;
% options.ipopt.acceptable_constr_viol_tol = tol;
% options.ipopt.print_level = 5;
% options.ipopt.grad_f = fgrad;
options.print_time        = 5;
% options.ipopt.max_iter    = 500;

Nx         =  numel(x0);
x_SX       =   SX.sym('x',Nx,1);
constraint = g(x_SX);
%% visualization 
% resolution = 100;
% s_qmin = -1;
% s_qmax = 2.2;
% s_q_step = linspace(s_qmin,s_qmax,resolution);
% b4_l= zeros(resolution,1); % store the solution under gridding
% eq_sq = @(x) obj_q(x);
% g_ext = @(x) vertcat(g(x),eq_sq(x));
% for i = 1:resolution
%     lbg_ext = vertcat(lbg,s_q_step(i));
%     ubg_ext = vertcat(ubg,s_q_step(i));
%     f_grid = @(x) obj_p(x); % p_{k,l} of PCC
%     constraint_grid = g_ext(x_SX);
%     objective = f_grid(x_SX);
%     nlp = struct('x',x_SX,'f',objective,'g',constraint_grid);
%     S   = nlpsol('solver','ipopt', nlp,options);
%     sol = S('x0', x0,'lbg', lbg_ext,'ubg', ubg_ext,...
%             'lbx', lbx, 'ubx', ubx);
%     stats = S.stats();
%     if strcmp(stats.return_status, 'Solve_Succeeded')
%         xopt2= full(sol.x);
%         obj_p_opt = xopt2(2*Nbus+id_gen_slack);
%         b4_l(i) = obj_p_opt;
%     end
% end
% s_q_step(b4_l == 0) =[];
% b4_l(b4_l == 0) =[];
% b4_l = [b4_l,s_q_step'];
%% visulization part2
resolution = 100;
s_pmin = -2;
s_pmax = 2;
s_p_step = linspace(s_pmin,s_pmax,resolution);
b4= zeros(resolution,1); % store the solution under gridding
eq_sp = @(x) obj_p(x);
g_ext = @(x) vertcat(g(x),eq_sp(x));
for i = 1:resolution
    lbg_ext = vertcat(lbg,s_p_step(i));
    ubg_ext = vertcat(ubg,s_p_step(i));
    f_grid = @(x) obj_q(x); % p_{k,l} of PCC
    constraint_grid = g_ext(x_SX);
    objective = f_grid(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint_grid);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg_ext,'ubg', ubg_ext,...
            'lbx', lbx, 'ubx', ubx);
    stats = S.stats();
    if strcmp(stats.return_status, 'Solve_Succeeded')
        xopt2= full(sol.x);
        obj_q_opt = xopt2(2*Nbus+Ngen+id_gen_slack);
        b4(i) = obj_q_opt;
    end
end
s_p_step(b4 == 0) =[];
b4(b4 == 0) =[];
b4 = [s_p_step',b4];
%% sampling points
i = 1;
Points = zeros(8,2);
obj_values = zeros(8,1);
for c1 = -1:1
    for c2 = -1:1
        if c1== 0 && c2 ==0
            continue;
        else
            f_samp   = @(x)c1*obj_p(x)+c2*obj_q(x); % p_{k,l}, q_{k,l} of PCC
            objective = f_samp(x_SX);
            nlp = struct('x',x_SX,'f',objective,'g',constraint);
            S   = nlpsol('solver','ipopt', nlp,options);
            sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                    'lbx', lbx, 'ubx', ubx);
            xopt= full(sol.x);
            obj_values(i)= full(sol.f);
            obj_p_opt = obj_p(xopt);
            obj_q_opt = obj_q(xopt);
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
gridding_p= zeros(resolution,2); % store the solution under gridding
eq_sp = @(x) obj_p(x);
g_ext = @(x)vertcat(g(x),eq_sp(x));
for c3 = [-1, 1]
    for i = 1:resolution
        lbg_ext = vertcat(lbg,s_p_step(i));
        ubg_ext = vertcat(ubg,s_p_step(i));
        f_grid = @(x)c3*obj_q(x); % q_{k,l} of PCC
        x_SX = SX.sym('x',Nx,1);
        constraint_grid = g_ext(x_SX);
        objective = f_grid(x_SX);
        nlp = struct('x',x_SX,'f',objective,'g',constraint_grid);
        S   = nlpsol('solver','ipopt', nlp,options);
        sol = S('x0', x0,'lbg', lbg_ext,'ubg', ubg_ext,...
                'lbx', lbx, 'ubx', ubx);
        xopt2= full(sol.x);
        obj_q_opt = obj_q(xopt2);
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
eq_sq = @(x) obj_q(x);
g_ext = @(x)vertcat(g(x),eq_sq(x));
for c3 = [-1, 1]
    for i = 1:resolution
        lbg_ext = vertcat(lbg,s_q_step(i));
        ubg_ext = vertcat(ubg,s_q_step(i));
        f_grid = @(x)c3*obj_p(x); % p_{k,l} of PCC
        x_SX = SX.sym('x',Nx,1);
        constraint_grid = g_ext(x_SX);
        objective = f_grid(x_SX);
        nlp = struct('x',x_SX,'f',objective,'g',constraint_grid);
        S   = nlpsol('solver','ipopt', nlp,options);
        sol = S('x0', x0,'lbg', lbg_ext,'ubg', ubg_ext,...
                'lbx', lbx, 'ubx', ubx);
        xopt2= full(sol.x);
        obj_p_opt = obj_p(xopt2);
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
title('Feasible Region of Slack Bus(exact)'); % 图像标题
grid on;            % 显示网格
toc
%% cost plot