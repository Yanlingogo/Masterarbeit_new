function [mpc, margin_opt] = margin_search(mpc)
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
    
    entries_pf{1} = 1:Nbus;                        % vmag
    entries_pf{2} = (Nbus+1):2*Nbus;               % vang
    entries_pf{3} = (2*Nbus+1):(2*Nbus+Ngen);      % Pg
    entries_pf{4} = (2*Nbus+Ngen+1):2*(Nbus+Ngen); % Qg 
    entries_pf{5} = 2*(Nbus+Ngen)+1:2*(Nbus+Ngen)+5;% slack_v/a/p/q/line
    entries_pf{6} = 2*(Nbus+Ngen)+6;               % margin
    
    baseMVA     = mpc.baseMVA;              % baseMVA
    %cost_param  = mpc.gencost(:,5:end);     % objective coefficients
    vmax        = mpc.bus(:,VMAX);                                            
    vmin        = mpc.bus(:,VMIN);
    Phimax      = mpc.branch(:,ANGMAX)/180*pi;
    Phimin      = mpc.branch(:,ANGMIN)/180*pi;
    Pgmin       = mpc.gen(:,PMIN)/baseMVA; %Pgmin(id_gen_slack) = -inf;  
    Qgmin       = mpc.gen(:,QMIN)/baseMVA; %Qgmin(id_gen_slack) = -inf;
    Pgmax       = mpc.gen(:,PMAX)/baseMVA; %Pgmax(id_gen_slack) = inf;
    Qgmax       = mpc.gen(:,QMAX)/baseMVA; %Qgmax(id_gen_slack) = inf;
    Fmax        = mpc.branch(:,RATE_A)/baseMVA;
    if Fmax ==0 
        mpc.branch(:,RATE_A) = 10e3*ones(size(Fmax));
        Fmax = 10e3*ones(size(Fmax));
    end
    
    Pd          = mpc.bus(:,PD)/baseMVA;   
    Qd          = mpc.bus(:,QD)/baseMVA;
    Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);

    %% lower & upper bounds
    lbx         = vertcat(-inf(2*Nbus+2*Ngen+5,1),-inf);     
    ubx         =  inf(2*Nbus+2*Ngen+6,1);
    %% initial state x0
    vang0       = mpc.bus(:,VA)/180*pi;
    vmag0       = mpc.bus(:,VM);
    Pg0         = mpc.gen(:,PG)/baseMVA;
    Qg0         = mpc.gen(:,QG)/baseMVA;
    x0          = vertcat(vmag0, vang0, Pg0, Qg0, zeros(6,1));

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
    % reference for slack bus
    vmag_ref = entries_pf{1}(id_slack);
    vang_ref = entries_pf{2}(id_slack);
    v_ref    = vertcat(vang_ref,vmag_ref);
    eq_ref   = @(x)x(v_ref);

    % power flow equation
    eq_pf          = @(x)create_local_power_flow_equation_pol(x(entries_pf{1}),x(entries_pf{2}),...
        x(entries_pf{3}),x(entries_pf{4}),Gbus,Bbus,Pd,Qd,Cg);
    %test_pf        = eq_pf(x0);
    Npf            = numel(eq_pf(x0));

    % slack inequality constraints 
    vmag_slack     =@(x)create_vmag_slack(x(entries_pf{1}),x(entries_pf{5}(1)),vmax,vmin,id_nslack); % 2*Nbus
    phi_slack      =@(x)create_phi_slack(x(entries_pf{2}),x(entries_pf{5}(2)),C,Phimax,Phimin); % 2*Nbranch
    gen_slack      =@(x)create_gen_slack(x(entries_pf{3}),x(entries_pf{4}),... 4*Ngen
        x(entries_pf{5}(3)),x(entries_pf{5}(4)),Pgmax,Pgmin,Qgmax,Qgmin);
    line_slack     =@(x)create_line_slack(x(entries_pf{1}), x(entries_pf{2}),... 2*Nbranch
        x(entries_pf{5}(5)), Gf, Bf, Gt, Bt, Fmax, from_bus,to_bus);
    margin_cons    =@(x)create_margin(x(entries_pf{5}(1)),x(entries_pf{5}(2)),... 5
        x(entries_pf{5}(3)),x(entries_pf{5}(4)),x(entries_pf{5}(5)),...
        x(entries_pf{6}));
    % slack bus test 
    Test_P = entries_pf{3}(id_gen_slack);
    Test_Q = entries_pf{4}(id_gen_slack);
    slack_test = vertcat(Test_P,Test_Q);
    eq_slack = @(x)x(slack_test);
    
    g = @(x)vertcat(vmag_slack(x),phi_slack(x),gen_slack(x),line_slack(x),margin_cons(x),...
        eq_pf(x),eq_ref(x),eq_slack(x));
    % upper and lower bound on g
    lbg = vertcat(-inf(2*Nbus+4*Nbranch+4*Ngen+3,1), zeros(Npf+1,1),1.0,Pg0(id_gen_slack),Qg0(id_gen_slack));
    ubg = vertcat(zeros(2*Nbus+4*Nbranch+4*Ngen+3,1),zeros(Npf+1,1),1.0,Pg0(id_gen_slack),Qg0(id_gen_slack));
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
    x_SX       =   SX.sym('x',Nx,1);
    constraint = g(x_SX);

    max_margin = @(x) -x(entries_pf{6});
    obj = max_margin(x_SX);
    nlp = struct('x',x_SX,'f',obj,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    xopt= full(sol.x);
    margin_opt = xopt(end);

    mpc.bus(:,VM) = xopt(1:Nbus);
    mpc.bus(:,VA) = xopt(Nbus+1:2*Nbus);
    mpc.gen(:,PG) = xopt(2*Nbus+1:2*Nbus+Ngen);
    mpc.gen(:,QG) = xopt(2*Nbus+Ngen+1:2*Nbus+2*Ngen);
end
    

