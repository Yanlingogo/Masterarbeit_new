function [Points,time] = solver_optimization(mpc, delta)
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
vmax        = mpc.bus(:,VMAX);                                            
vmin        = mpc.bus(:,VMIN);
Phimax      = mpc.branch(:,ANGMAX)/180*pi;
Phimin      = mpc.branch(:,ANGMIN)/180*pi;
Pgmin       = mpc.gen(:,PMIN)/baseMVA;   
Qgmin       = mpc.gen(:,QMIN)/baseMVA; 
Pgmax       = mpc.gen(:,PMAX)/baseMVA; 
Qgmax       = mpc.gen(:,QMAX)/baseMVA; 
Fmax        = mpc.branch(:,RATE_A)/baseMVA;


Pd          = mpc.bus(:,PD)/baseMVA;   
Qd          = mpc.bus(:,QD)/baseMVA;
Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);

%% lower & upper bounds
lbx         = [vmin;-pi*ones(Nbus,1);Pgmin;Qgmin];       
ubx         = [vmax; pi*ones(Nbus,1);Pgmax;Qgmax];
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


% line Limit
if any(Fmax)
    ineq_line = @(x)create_local_branch_limit_rec(x(entries_pf{1}),...
        x(entries_pf{2}), Gf, Bf, Gt, Bt, Fmax, from_bus,to_bus);
    idx_limit = find(Fmax);
    Nlimit = numel(ineq_line(x0));
    g = @(x)vertcat(ineq_angle(x),ineq_line(x),eq_pf(x),eq_ref(x)); %without generator bound
else
    ineq_line = [];
    idx_limit = [];
    Nlimit   = 0;
    g   = @(x)vertcat(ineq_angle(x),eq_pf(x),eq_ref(x));
end

obj_p = @(x) x(entries_pf{3}(id_gen_slack));
obj_q = @(x) x(entries_pf{4}(id_gen_slack));


lbg = vertcat(Phimin, -inf*ones(Nlimit,1), zeros(Npf+1,1),vmin(id_slack));
ubg = vertcat(Phimax, zeros(Npf+Nlimit+1,1),vmin(id_slack));
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
options.print_time        = 0;
options.ipopt.max_iter    = 100;

Nx         =  numel(x0);
x_SX       =   SX.sym('x',Nx,1);
constraint = g(x_SX);
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
%% adaptive searching angle 
sigma = delta;
% Points_c = [0 0];
Points_c = mean(Points);
Points = Orderpoints(Points,Points_c);
[Points, Dist] = distance_calculation(Points);

max_iter = 50;
iter = 0;
while any(Dist > sigma)
    index = find (Dist > sigma);
    for i = index'
        delta_p = ((Points(i,1)+Points(i+1,1))/2-Points_c(1));
        delta_q = ((Points(i,2)+Points(i+1,2))/2-Points_c(2));
        Phi = atan2(delta_q,delta_p);
        c1 = cos(Phi);
        c2 = sin(Phi); 

        angle_cons = @(x) [c2 -c1]*([obj_p(x) obj_q(x)]-Points_c)';
        g_ext = @(x) vertcat(g(x),angle_cons(x));
        cons_ext = g_ext(x_SX);
        lbg_ext = vertcat(lbg,0);
        ubg_ext = vertcat(ubg,0);

        % Obj_ext = @(x) [c1 c2]*([obj_p(x) obj_q(x)]-Points_c)';
        Obj_ext = @(x) (obj_p(x)-Points_c(1))/c1+(obj_q(x)-Points_c(2))/c2;
        objective = Obj_ext(x_SX);

        nlp = struct('x',x_SX,'f',-objective,'g',cons_ext);
        S   = nlpsol('solver','ipopt', nlp,options);
        sol = S('x0', x0,'lbg', lbg_ext,'ubg', ubg_ext,...
                'lbx', lbx, 'ubx', ubx);
        xopt= full(sol.x);
        obj_values(i)= full(sol.f);
        obj_p_opt = obj_p(xopt);
        obj_q_opt = obj_q(xopt);
        Points = [Points;[obj_p_opt,obj_q_opt]];
    end
    Points = Orderpoints(Points,Points_c);
    [Points,Dist] = distance_calculation(Points);
    iter = iter+1;
    if iter == max_iter
        disp('Threshold may not be met');
        break;     
    end
end
time = toc;
end

