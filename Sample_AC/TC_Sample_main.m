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

mpc = ext2int(loadcase('case33_modified'));
mpc = ext2int(mpc);
% Tiem periode
T = 6; 
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

entries_pf{1} = 1:Nbus*T;                            % vmag
entries_pf{2} = (Nbus*T+1):2*Nbus*T;% vang
entries_pf{3} = (2*Nbus*T+1):(2*Nbus*T+Ngen*T);      % Pg
entries_pf{4} = (2*Nbus*T+Ngen*T+1):2*(Nbus+Ngen)*T; % Qg 

baseMVA     = mpc.baseMVA;              % baseMVA
cost_param  = mpc.gencost(:,5:end);     % objective coefficients
vmax        = mpc.bus(:,VMAX);                                            
vmin        = mpc.bus(:,VMIN);
Phimax      = mpc.branch(:,ANGMAX)/180*pi;
Phimin      = mpc.branch(:,ANGMIN)/180*pi;
Pgmin       = mpc.gen(:,PMIN)/baseMVA; %Pgmin(id_gen_slack) = -inf;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; %Qgmin(id_gen_slack) = -inf;
Pgmax       = mpc.gen(:,PMAX)/baseMVA; %Pgmax(id_gen_slack) = inf;
Qgmax       = mpc.gen(:,QMAX)/baseMVA; %Qgmax(id_gen_slack) = inf;
Fmax        = mpc.branch(:,RATE_A)/baseMVA;


Pd          = mpc.bus(:,PD)/baseMVA; Pd = repmat(Pd,1,T); 
Qd          = mpc.bus(:,QD)/baseMVA; Qd = repmat(Qd,1,T);
Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);

%% lower & upper bounds
lbx         = [repmat(vmin,T,1);-pi*ones(Nbus*T,1);repmat(Pgmin,T,1);repmat(Qgmin,T,1)];       
ubx         = [repmat(vmax,T,1); pi*ones(Nbus*T,1);repmat(Pgmax,T,1);repmat(Qgmax,T,1)];
%% initial state x0
vang0       = mpc.bus(:,VA)/180*pi;
vmag0       = mpc.bus(:,VM);
U0          = vmag0.*cos(vang0);
W0          = vmag0.*sin(vang0);
Pg0         = mpc.gen(:,PG)/baseMVA;
Qg0         = mpc.gen(:,QG)/baseMVA;
x0      = vertcat(repmat(vmag0,T,1),repmat(vang0,T,1), repmat(Pg0,T,1), repmat(Qg0,T,1));
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
id_slack_set   = id_slack+linspace(0,Nbus*(T-1),T);
vmag_ref       = entries_pf{1}(id_slack_set)';
vang_ref       = entries_pf{2}(id_slack_set)';
v_ref          = vertcat(vang_ref,vmag_ref);
eq_ref         = @(x)x(v_ref);

% power flow equation
pflow = cell(T, 1); 

% The first loop builds a series of anonymous functions
for i = 1:T
    eq_pf_subset = @(x) create_local_power_flow_equation_pol(...
        x(entries_pf{1}((i-1)*Nbus+1:i*Nbus)),...
        x(entries_pf{2}((i-1)*Nbus+1:i*Nbus)),...
        x(entries_pf{3}((i-1)*Ngen+1:i*Ngen)),...
        x(entries_pf{4}((i-1)*Ngen+1:i*Ngen)),...
        Gbus, Bbus, Pd(:,i), Qd(:,i), Cg);
    pflow{i} = eq_pf_subset;
end

% Initialize eq_pf to an empty anonymous function
eq_pf = @(x) [];  

% In the second loop, all anonymous functions are spliced together
for i = 1:T
    eq_pf = @(x) vertcat(eq_pf(x), pflow{i}(x));
end

% Test and get dimensions
test_pf = eq_pf(x0);
Npf = numel(test_pf);


% upper & lower bound for voltage angle
C_T = [];
for i = 1:T
    C_T = blkdiag(C_T,C);
end
ineq_angle     = @(x) C_T*x(entries_pf{2});

test_angle= ineq_angle(x0);
test_ref  = eq_ref(x0);

% bound between pg_k and qg_k
% alpha = 0.8;  % based on experience

%ineq_genpq = @(x)create_bound_PQ(x(entries_pf{3}),x(entries_pf{4}),alpha,id_gen_nslack,id_gen_slack);
% line Limit
if any(Fmax)
    S_limit = cell(T, 1); 
    for i = 1:T
        ineq_line = @(x)create_local_branch_limit_rec(x(entries_pf{1}),...
        x(entries_pf{2}), Gf, Bf, Gt, Bt, Fmax, from_bus,to_bus);
        S_limit{i} = ineq_line;
    end
    ineq_Slimit = @(x) [];
    for i = 1:T
        ineq_S = @(x) vertcat(ineq_S(x), S_limit{i}(x));
    end

    Nlimit = numel(ineq_S(x0));
    %g = @(x)vertcat(ineq_voltage(x),ineq_genp(x),ineq_genq(x),ineq_line(x),ineq_genpq(x),eq_pf(x),eq_ref(x));
    g = @(x)vertcat(ineq_angle(x),ineq_line(x),eq_pf(x),eq_ref(x)); %without generator bound
else
    ineq_S = [];
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

%lbg = vertcat(vmin, Pgmin, Qgmin, -inf*ones(Nlimit+2*Ngen-2*Nslack,1), zeros(Npf+1,1));
lbg = vertcat(repmat(Phimin,T,1), -inf*ones(Nlimit,1), zeros(Npf+T,1),ones(T,1));
%ubg = vertcat(vmax, Pgmax, Qgmax, zeros(Npf+Nlimit+2*Ngen+1-2*Nslack,1));
ubg = vertcat(repmat(Phimax,T,1), zeros(Npf+Nlimit+T,1), ones(T,1));
% based on eqauation(19): but no line limits,
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
options.ipopt.max_iter    = 100;

Nx         =  numel(x0);
%% Searching Process
results = cell(1,T);
for t = 1:T
    x_SX       =  SX.sym('x',Nx,1);
    constraint =  g(x_SX);
    obj_p = @(x) x(entries_pf{3}(id_gen_slack+(t-1)*Ngen));
    obj_q = @(x) x(entries_pf{4}(id_gen_slack+(t-1)*Ngen));
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
                S   = nlpsol('solver','ipopt', nlp, options);
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
    resolution = 5;
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
    %% output for each periode
    Points_tot = [Points;[s_p_step';s_p_step'],[gridding_p(:,1);gridding_p(:,2)];...
        [gridding_q(:,1);gridding_q(:,2)],[s_q_step';s_q_step']];
    centroid = mean(Points_tot, 1);
    angles = atan2(Points_tot(:,2) - centroid(2), Points_tot(:,1) - centroid(1));
    [~, order] = sort(angles);
    results{t} = Points_tot(order, :)';
end
%% output 画出的图像中，P/Q > 0 表示向DSO中注入功率
figure;

% 设置颜色地图，以便每个多边形有不同的颜色
colors = jet(length(results)); % 'jet' 是 MATLAB 的一个颜色地图

% 循环通过每个 cell 绘制多边形
for i = 1:length(results)
    % 获取当前 cell 的坐标点
    xy = (results{i})';
    centroid = mean(xy, 1);
    angles = atan2(xy(:,2) - centroid(2), xy(:,1) - centroid(1));
    [~, order] = sort(angles);
    sortedPoints = xy(order, :);
    % 假设 xy 是一个 n x 2 的矩阵，其中第一列是 x 坐标，第二列是 y 坐标
    x = sortedPoints(:, 1);
    y = sortedPoints(:, 2);
    
    % 在三维空间中添加多边形。Z 坐标由 i 控制，以将它们堆叠起来
    z = ones(size(x)) * i;
    patch(z, x, y, colors(i, :), 'EdgeColor', 'none');
    
    % 可以选择添加一些透明度
    alpha(0.5);

    hold on; % 保持当前图形，以便在上面绘制
    plot3(z, x, y, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i, :));
end

% 调整视角
ylim([-2 2]);
zlim([-2 2]);

set(gca, 'YDir','reverse')
view(3);

% 添加轴标签
xlabel('T(time period)');
ylabel('P(p.u.)');
zlabel('Q(p.u.)');
gamma = 0;
title(['Feasible region with temporal coupling (\gamma=', num2str(gamma), ')']);

% 优化图形显示
grid on;
axis normal;
rotate3d on; % 允许使用鼠标旋转视图


% plot(Points(:,1),Points(:,2),'rx');
% hold on;
% Points_tot = [Points;[s_p_step';s_p_step'],[gridding_p(:,1);gridding_p(:,2)];...
%     [gridding_q(:,1);gridding_q(:,2)],[s_q_step';s_q_step']];
% centroid = mean(Points_tot, 1);
% angles = atan2(Points_tot(:,2) - centroid(2), Points_tot(:,1) - centroid(1));
% [~, order] = sort(angles);
% sortedPoints = Points_tot(order, :);
% plot([sortedPoints(:,1); sortedPoints(1,1)], [sortedPoints(:,2); sortedPoints(1,2)], '-bo');
% xlabel('P/p.u.');   % X 轴标签
% ylabel('Q/p.u.');   % Y 轴标签
% title('Feasible Region of Slack Bus(exact)'); % 图像标题
% grid on;            % 显示网格

%% cost plot
