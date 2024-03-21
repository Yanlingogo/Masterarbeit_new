function [result, margin_opt] = margin_search(mpc)
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
    
    id_bus      = mpc.bus(:,BUS_I);
    id_gen      = mpc.gen(:,GEN_BUS);
    Nbus        = size(mpc.bus,1);
    Ngen        = numel(id_gen);
    Nbranch     = size(mpc.branch,1);
    
    % idx
    id_slack =  find(mpc.bus(:,BUS_TYPE) == REF);
    id_nslack = find(mpc.bus(:,BUS_TYPE) ~= REF);
    Nslack = numel(id_slack);
    id_gen_nslack = find(id_gen ~= id_slack);
    id_gen_slack = find(id_gen == id_slack);
    
    entries_pf{1} = 1:Nbus;                                             % U
    entries_pf{2} = (Nbus+1):Nbus+Nbranch;                              % Pij
    entries_pf{3} = Nbus+Nbranch+1:Nbus+2*Nbranch;                      % Qij
    entries_pf{4} = Nbus+2*Nbranch+1:Nbus+2*Nbranch+Ngen;               % Pg
    entries_pf{5} = Nbus+2*Nbranch+Ngen+1:Nbus+2*Nbranch+2*Ngen;        % Qg 
    entries_pf{6} = Nbus+2*Nbranch+2*Ngen+1:Nbus+2*Nbranch+2*Ngen+3;    % slack_U/P/Q
    entries_pf{7} = Nbus+2*Nbranch+2*Ngen+4;    % margin
    
    baseMVA     = mpc.baseMVA;              % baseMVA
    %cost_param  = mpc.gencost(:,5:end);     % objective coefficients
    vmax        = mpc.bus(:,VMAX);                                            
    vmin        = mpc.bus(:,VMIN);
    Umin        = vmin.^2;
    Umax        = vmax.^2;
    Phimax      = mpc.branch(:,ANGMAX)/180*pi;
    Phimin      = mpc.branch(:,ANGMIN)/180*pi;
    Pgmin       = mpc.gen(:,PMIN)/baseMVA;   
    Qgmin       = mpc.gen(:,QMIN)/baseMVA; 
    Pgmax       = mpc.gen(:,PMAX)/baseMVA; 
    Qgmax       = mpc.gen(:,QMAX)/baseMVA; 
    Fmax        = mpc.branch(:,RATE_A)/baseMVA;
    if Fmax ==0 
        mpc.branch(:,RATE_A) = 10e3*ones(size(Fmax));
        Fmax = 10e3*ones(size(Fmax));
    end
    
    Pd          = mpc.bus(:,PD)/baseMVA;   
    Qd          = mpc.bus(:,QD)/baseMVA;
    Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);

    %% lower & upper bounds
    lbx         = -inf(Nbus+2*Nbranch+2*Ngen+4,1);       
    ubx         = inf(Nbus+2*Nbranch+2*Ngen+4,1);
    % lbx = vertcat(0.8*ones(size(Umin)),-inf(2*Nbranch+Ngen,1),Qgmin,-inf(4,1));
    % ubx = vertcat(2*ones(size(Umax)),inf(2*Nbranch+Ngen,1),Qgmax,inf(4,1));
    % lbx = vertcat(0.8*ones(size(Umin)),-inf(2*Nbranch,1),2*Pgmin,-inf(Ngen+4,1));
    % ubx = vertcat(1.5*ones(size(Umax)),inf(2*Nbranch,1),2*Pgmax,inf(Ngen+4,1));
   
    %% initial state x0
    vmag0       = (mpc.bus(:,VM)).^2;
    Pg0         = mpc.gen(:,PG)/baseMVA;
    Qg0         = mpc.gen(:,QG)/baseMVA;
    x0          = vertcat(vmag0, zeros(2*Nbranch,1),Pg0, Qg0, zeros(4,1));

    %% equality & inequality constraints 
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
    % branch info
    branch_r   = mpc.branch(:,BR_R);
    branch_x   = mpc.branch(:,BR_X);


    %% LinDistFlow
    % voltage constraints
    % volt_eq    = C*U -2*(branch_r .* Pij + branch_x .* Qij);
    volt_eq    = @(x) C*x(entries_pf{1}) - 2*(branch_r .* x(entries_pf{2})+branch_x .* x(entries_pf{3}));
    % power balance
    % pf_p_eq    = Cg*Pg - Pd - C'*Pij;
    % pf_q_eq    = Cg*Qg - Qd - C'*Qij;
    pf_p_eq    = @(x) Cg*x(entries_pf{4}) - Pd - C'*x(entries_pf{2});
    pf_q_eq    = @(x) Cg*x(entries_pf{5}) - Qd - C'*x(entries_pf{3});
    % ref bus
    ref_eq     = @(x) x(entries_pf{1}(id_slack)) - mpc.bus(id_slack, VM).^2;
    % ref PQ 
    ref_pq     = @(x) create_PQtest(x(entries_pf{4}),x(entries_pf{5}),Pg0,Qg0,id_gen_slack);

    % slack inequality constraints 
    U_slack     =@(x)create_U_slack(x(entries_pf{1}),x(entries_pf{6}(1)),Umax,Umin,id_nslack); % 2*Nbus

    gen_slack      =@(x)create_gen_slack(x(entries_pf{4}),x(entries_pf{5}),... 4*Ngen
        x(entries_pf{6}(2)),x(entries_pf{6}(3)),Pgmax,Pgmin,Qgmax,Qgmin);

    margin_cons    =@(x)create_margin(x(entries_pf{6}(1)),x(entries_pf{6}(2)),... 5
        x(entries_pf{6}(3)),x(entries_pf{7}));

    
    g = @(x)vertcat(volt_eq(x),pf_p_eq(x),pf_q_eq(x),ref_eq(x),ref_pq(x),U_slack(x),gen_slack(x),margin_cons(x));
    % upper and lower bound on g
    g_test = g(x0);
    num_g = numel(g_test);

    lbg = vertcat(zeros(Nbranch+2*Nbus+3,1),-inf((num_g-(Nbranch+2*Nbus+3)),1));
    ubg = vertcat(zeros(num_g,1));
    %% solver options
    import casadi.*
    % tolerance
    tol        = 1e-6;
    options.ipopt.tol             = tol;
    options.ipopt.constr_viol_tol = tol;
    options.ipopt.compl_inf_tol   = tol;
    options.ipopt.acceptable_tol  = tol;
    options.ipopt.acceptable_constr_viol_tol = tol;
    options.ipopt.print_level = 0;
    % options.ipopt.grad_f = fgrad;
    options.print_time        = 5;
    options.ipopt.max_iter    = 100;

    Nx         =  numel(x0);
    x_SX       =  SX.sym('x',Nx,1);
    constraint = g(x_SX);

    max_margin = @(x) -x(entries_pf{7});
    obj = max_margin(x_SX);
    nlp = struct('x',x_SX,'f',obj,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    xopt= full(sol.x);
    
    margin_opt = xopt(end);
    result = struct();
    result.U = xopt(1:Nbus);
    result.pg = xopt(Nbus+2*Nbranch+1:Nbus+2*Nbranch+Ngen);
    result.qg = xopt(Nbus+2*Nbranch+Ngen+1:Nbus+2*Nbranch+2*Ngen);

end
    

