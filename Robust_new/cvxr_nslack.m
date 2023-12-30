function [margin_opt] = cvxr(mpc,mpc_target)
    if nargin<3
        phase_shift = false;
    end
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
%% constant paramets 
    % number of all elements
    Nbus = size(mpc.bus,1);
    Ngen = size(mpc.gen,1);
    idx_load = find((mpc.bus(:,PD).^2+mpc.bus(:,QD).^2) ~= 0);
    Nload = numel(idx_load);
    Npq = sum(mpc.bus(:,BUS_TYPE)==1);
    Nbranch = size(mpc.branch,1);
    % parameters 
    bus = mpc.bus;
    gen = mpc.gen;
    branch = mpc.branch;

    bus_type = bus(:,BUS_TYPE);
    idx_slack = find(bus_type == 3);
    idx_nslack = find(bus_type ~= 3);Nnslack = numel(idx_nslack);
    idx_sgen = find(gen(:,GEN_BUS)==idx_slack);

    idx_fr = branch(:,F_BUS);
    idx_to = branch(:,T_BUS);
    id_gen = gen(:,GEN_BUS);
    %idx_load mentioned before
    idx_pq = find(bus_type == 1);
    onoff_pq = zeros(Nbus,1); onoff_pq(idx_pq) = 1;
    idx_pvs = find(bus_type ~=1);

    %gencost = mpc.gencost(:,COST:end);

    E_fr = sparse(idx_fr, 1:Nbranch, 1, Nbus, Nbranch);
    E_to = sparse(idx_to, 1:Nbranch, 1, Nbus, Nbranch);
    E = E_fr - E_to;
    Cg = sparse(id_gen,1:Ngen,1,Nbus,Ngen);
    Cl = sparse(idx_load,1:Nload,1,Nbus,Nload);

    v0 = bus(:,VM);
    theta0 = bus(:,VA);
    Phi0 = E'*theta0;
    % alpha0 = gen(:,end);
    % delta0 = mpc.delta;
    pg0 = (gen(:,GEN_STATUS).*gen(:,PG))/mpc.baseMVA;
    qg0 = (gen(:,GEN_STATUS).*gen(:,QG))/mpc.baseMVA;
    pl0 = bus(idx_load,PD)/mpc.baseMVA;
    ql0 = bus(idx_load,QD)/mpc.baseMVA;
    p_inj0 = Cg*pg0-Cl*pl0;
    q_inj0 = Cg*qg0-Cl*ql0;
    % power_factor = ql0./(pl0+1e-3); 注意考虑这三项是否有用
    % limitation of the network
    vmax = bus(:,VMAX);
    vmin = bus(:,VMIN);
    Phi_max = deg2rad(branch(:,ANGMAX));
    Phi_min = deg2rad(branch(:,ANGMIN));
    pg_max = (gen(:,GEN_STATUS).*gen(:,PMAX))/mpc.baseMVA;
    pg_min = (gen(:,GEN_STATUS).*gen(:,PMIN))/mpc.baseMVA;
    qg_max = (gen(:,GEN_STATUS).*gen(:,QMAX))/mpc.baseMVA;
    qg_min = (gen(:,GEN_STATUS).*gen(:,QMIN))/mpc.baseMVA;
    alpha_p = gen(:,end-1);
    alpha_q = gen(:,end);
    s_line_max = branch(:,RATE_A)/mpc.baseMVA;

    Sigma0 = mpc.uncertainty.Sigma; 
    gamma0 = mpc.uncertainty.gamma0;
    % Jacobian 
    J1 = -diag(v0(idx_fr).*v0(idx_to).*sin(E'*theta0))*E'; %Nbranch*Nbus
    J2 = diag(v0(idx_to).*cos(E'*theta0))*E_fr'+diag(v0(idx_fr).*cos(E'*theta0))*E_to';
    J3 = diag(v0(idx_fr).*v0(idx_to).*cos(E'*theta0))*E';
    J4 = diag(v0(idx_to).*sin(E'*theta0))*E_fr'+diag(v0(idx_fr).*sin(E'*theta0))*E_to';
    D = diag(v0);
    J_psi0 = [zeros(2*Nbus,2*size(idx_nslack,1)),[Cg*alpha_p zeros(Nbus,1);zeros(Nbus,1) Cg*alpha_q];
              J1(:,idx_nslack), J2(:,idx_nslack),zeros(Nbranch,2);
              J3(:,idx_nslack), J4(:,idx_nslack),zeros(Nbranch,2);
              zeros(Nbus,size(idx_nslack,1)), 2*D(:,idx_nslack),zeros(Nbus,2)];% (3*Nbus+2Nbranch)*(2Nslack+2)
    % original value of residues and basis function
    g_vvcos0 = v0(idx_fr).*v0(idx_to).*cos(Phi0) - onoff_pq(idx_fr).*v0(idx_fr).*v0(idx_to).*cos(Phi0)...
        - v0(idx_fr).*onoff_pq(idx_to).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi0;
    g_vvsin0 = v0(idx_fr).*v0(idx_to).*sin(Phi0) - onoff_pq(idx_fr).*v0(idx_fr).*v0(idx_to).*sin(Phi0)...
        - v0(idx_fr).*onoff_pq(idx_to).*v0(idx_to).*sin(Phi0) - v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi0;
    g_vv0 = v0.^2-2*onoff_pq.*v0.*v0;
    Psi_vvcos0 = v0(idx_fr).*v0(idx_to).*cos(Phi0);
    Psi_vvsin0 = v0(idx_fr).*v0(idx_to).*sin(Phi0);
    Psi_vv0 = v0.^2;

    [~,~,~,M,M_line] = makeYbus_cvxr(mpc,phase_shift); % M: 2Nbus*(3Nbus+2Nbranch)
    M_eq = M(:,:);% Nbus+Npq*(3Nbus+2Nbranch)
    M_ineq = M([idx_pvs',idx_pvs'+Nbus],:);

    J_inv_M = (M_eq*J_psi0)\M_eq; % (Nbus+Npq)*(3Nbus+2Nbranch) 

    %R = -[Cl,zeros(size(Cl(idx_pq,:)));zeros(size(Cl)),Cl(idx_pq,:)];% 2Npq*2Nload
    %J_inv_R = (M_eq*J_psi0)\R; % 2Npq*2Nload

    C = [E(idx_nslack,:)' zeros(Nbranch,Nnslack+2);
        zeros(Nnslack,Nnslack) eye(Nnslack) zeros(Nnslack,2);
        zeros(1,2*Nnslack) 1 0;
        zeros(1,2*Nnslack+1) 1];

    A = [eye(Nbranch+Nnslack+2);-eye(Nbranch+Nnslack+2)]*C; % (Nbranch+Nbus)*(2Npq)
    K = A*J_inv_M; % (2Nbranch+2Nbus)*(3Nbus+2Nbranch)
    %T = A*J_inv_R; % (Nbranch+Nbus)*2Nload

    K_plus  = max(K,zeros(size(K))).*(K>0);
    K_minus = min(K,zeros(size(K))).*(K<0);

    M_ineq_plus  = max(M_ineq, zeros(size(M_ineq))).*(M_ineq>0);
    M_ineq_minus = min(M_ineq, zeros(size(M_ineq))).*(M_ineq<0);
    M_line_plus  = max(M_line, zeros(size(M_line))).*(M_line>0);
    M_line_minus = min(M_line, zeros(size(M_line))).*(M_line<0);
    % meaning of zeta and xi, idx_pvs, power_factor ?
    %zeta = sqrt(sum((Cl(idx_pvs,:)*diag(power_factor)*chol(Sigma0,'upper')).^2,2));
    %xi = sqrt(sum((K(:,1:2*Nbus)*[Cl;Cl]*chol(Sigma0,'upper')).^2,2));
     xi0 = -T*[pl0;ql0];
     xi_gamma = sqrt(sum((-T*chol(Sigma0,'upper')).^2,2));

    Cl_zeta = blkdiag(Cl(idx_pvs,:),Cl(idx_pvs,:));
%     zeta = sqrt(sum((Cl_zeta*chol(Sigma0,'upper')).^2,2)); %==0, pvs nodes have loads
    zeta =0;
%% optimization problem
    entries{1}  = 1:Nbus;        %v_u
    entries{2}  = 1:Nbus;        %v_l
    entries{3}  = 1:Nbus;        %Delta_v_u
    entries{4}  = 1:Nbus;        %Delta_v_l
    entries{5}  = 1:Nbranch;     %Phi_u
    entries{6}  = 1:Nbranch;     %Phi_l
    entries{7}  = 1:Nbranch;     %Delta_Phi_u
    entries{8}  = 1:Nbranch;     %Delta_Phi_l
    %entries{9} = 1:4; % set of delta
    entries{9}  = 1:Nbus;        %g_u_pinj
    entries{10} = 1:Nbus;        %g_l_pinj
    entries{11} = 1:Nbus;        %g_u_qinj
    entries{12} = 1:Nbus;        %g_l_qinj  
    entries{13} = 1:Nbranch;     %g_u_vvcos %定义顺序与 julia不同
    entries{14} = 1:Nbranch;     %g_l_vvcos
    entries{15} = 1:Nbranch;     %g_u_vvsin
    entries{16} = 1:Nbranch;     %g_l_vvsin
    entries{17} = 1:Nbus;        %g_u_vv
    entries{18} = 1:Nbus;        %g_l_vv

    entries{19} = 1:Nbranch;     %Psi_u_vvcos
    entries{20} = 1:Nbranch;     %Psi_l_vvcos
    entries{21} = 1:Nbranch;     %Psi_u_vvsin
    entries{22} = 1:Nbranch;     %Psi_l_vvsin
    entries{23} = 1:Nbus;        %Psi_u_vv
    entries{24} = 1:Nbus;        %Psi_l_vv

    entries{25} = 1:8*Nbranch;   %D_vv_u_uu...
    entries{26} = 1:8*Nbranch;   %D_cos_u_u...
    entries{27} = 1:4*Nbranch;   %p_line_fr_u...
    
    entries{28} = 1:2*Ngen;      %pg_opt, qg_opt
    entries{29} = 1:2*Nload;     %pl_opt, ql_opt
    entries{30} = 1:5;             %gamma_opt, slack_pq, slack_qg, slack_line, margin_opt
    %% lower & upper bounds for the variables 
    lbx = [vmin;vmin;vmin-v0;vmin-v0;Phi_min;Phi_min;Phi_min-Phi0;Phi_max-Phi0;-inf(8*Nbus+28*Nbranch+2*Ngen+2*Nload+5,1)];
    ubx = [vmax;vmax;vmax-v0;vmax-v0;Phi_max;Phi_max;Phi_min-Phi0;Phi_max-Phi0; inf(8*Nbus+28*Nbranch+2*Ngen+2*Nload+5,1)];
    %% initial state x0
    x0 = vertcat(v0,v0,zeros(2*Nbus,1),Phi0,Phi0,zeros(8*Nbus+30*Nbranch,1),pg0,qg0,pl0,ql0,zeros(5,1));
    %% equality/inequality constraints
    % eq checked
    eq_cons = @(x) create_eq_constraints(v0,Phi0,id_gen,x(entries{1}),x(entries{2}),...
        x(entries{3}),x(entries{4}),x(entries{5}),x(entries{6}),x(entries{7}),x(entries{8})); % = 0
    % constraints from target network(equality constraints) checked
    eq_target = @(x) create_eq_target(x(entries{28}(1:Ngen)),x(entries{28}(Ngen+1:2*Ngen)),...
        x(entries{29}(1:Nload)),x(entries{29}(Nload+1:2*Nload)),x(entries{1}),mpc_target,idx_load,id_gen);

    % ineq/ convex restriction checked
    ineq_vv = @(x) create_ineq_vv(x(entries{25}(1:Nbranch)),x(entries{25}(Nbranch+1:2*Nbranch)),x(entries{25}(2*Nbranch+1:3*Nbranch)),x(entries{25}(3*Nbranch+1:4*Nbranch)),...
        x(entries{25}(4*Nbranch+1:5*Nbranch)),x(entries{25}(5*Nbranch+1:6*Nbranch)),x(entries{25}(6*Nbranch+1:7*Nbranch)),x(entries{25}(7*Nbranch+1:8*Nbranch)),...
        x(entries{3}),x(entries{4}),v0,idx_fr,idx_to);
    % checked
    ineq_cossin = @(x) create_ineq_cossin(x(entries{26}(1:Nbranch)),x(entries{26}(Nbranch+1:2*Nbranch)),x(entries{26}(2*Nbranch+1:3*Nbranch)),x(entries{26}(3*Nbranch+1:4*Nbranch)),...
        x(entries{26}(4*Nbranch+1:5*Nbranch)),x(entries{26}(5*Nbranch+1:6*Nbranch)),x(entries{26}(6*Nbranch+1:7*Nbranch)),x(entries{26}(7*Nbranch+1:8*Nbranch)),...
        x(entries{7}),x(entries{8}),Phi0); 
    % checked
    ineq_vvcos_sin = @(x)create_ineq_vvcos(x(entries{13}),x(entries{14}),x(entries{15}),x(entries{16}),...
        x(entries{25}(1:Nbranch)),x(entries{25}(Nbranch+1:2*Nbranch)),x(entries{25}(2*Nbranch+1:3*Nbranch)),x(entries{25}(3*Nbranch+1:4*Nbranch)),...
        x(entries{25}(4*Nbranch+1:5*Nbranch)),x(entries{25}(5*Nbranch+1:6*Nbranch)),x(entries{25}(6*Nbranch+1:7*Nbranch)),x(entries{25}(7*Nbranch+1:8*Nbranch)),...
        x(entries{26}(1:Nbranch)),x(entries{26}(Nbranch+1:2*Nbranch)),x(entries{26}(2*Nbranch+1:3*Nbranch)),x(entries{26}(3*Nbranch+1:4*Nbranch)),...
        x(entries{26}(4*Nbranch+1:5*Nbranch)),x(entries{26}(5*Nbranch+1:6*Nbranch)),x(entries{26}(6*Nbranch+1:7*Nbranch)),x(entries{26}(7*Nbranch+1:8*Nbranch)),...
        x(entries{1}),x(entries{2}),x(entries{5}),x(entries{6}),Psi_vvcos0,Psi_vvsin0,Phi0,v0,idx_fr,idx_to,onoff_pq);
    % checked
    ineq_gvv = @(x) create_ineq_gvv(x(entries{17}),x(entries{18}),...
        x(entries{1}),x(entries{2}),v0,onoff_pq,Psi_vv0); 
    % checked
    ineq_Psivv = @(x) create_ineq_Psivv(x(entries{19}),x(entries{20}),x(entries{21}),x(entries{22}),x(entries{23}),x(entries{24}),...
        x(entries{25}(1:Nbranch)),x(entries{25}(4*Nbranch+1:5*Nbranch)),x(entries{25}(5*Nbranch+1:6*Nbranch)),...
        x(entries{26}(1:Nbranch)),x(entries{26}(Nbranch+1:2*Nbranch)),x(entries{26}(2*Nbranch+1:3*Nbranch)),x(entries{26}(3*Nbranch+1:4*Nbranch)),...cos set
        x(entries{26}(4*Nbranch+1:5*Nbranch)),x(entries{26}(7*Nbranch+1:8*Nbranch)),...sin set
        x(entries{1}),x(entries{2}),x(entries{3}),x(entries{4}),... v_u/l set
        Psi_vvcos0,Psi_vvsin0,Psi_vv0,Phi0,v0,idx_fr,idx_to);

    ineq_inj = @(x)create_ineq_inj(x(entries{9}),x(entries{10}),x(entries{11}),x(entries{12}),...
        x(entries{27}(1:Nbranch)),x(entries{27}(Nbranch+1:2*Nbranch)),x(entries{27}(2*Nbranch+1:3*Nbranch)),x(entries{27}(3*Nbranch+1:4*Nbranch)),...
        x(entries{28}(1:Ngen)),x(entries{28}(Ngen+1:2*Ngen)),x(entries{29}(1:Nload)),x(entries{29}(Nload+1:2*Nload)),...
        x(entries{30}(1)),x(entries{30}(2)),x(entries{30}(3)),x(entries{30}(4)),x(entries{30}(5)),...
        Cg,Cl,pg_max,pg_min,qg_max,qg_min,s_line_max);

    ineq_residues = @(x)create_ineq_residues(x(entries{9}),x(entries{10}),x(entries{11}),x(entries{12}),x(entries{13}),x(entries{14}),x(entries{15}),x(entries{16}),x(entries{17}),x(entries{18}),...
        x(entries{30}(1)),K_plus,K_minus,xi0,xi_gamma,x(entries{5}),x(entries{6}),x(entries{1}),x(entries{2}),idx_nslack);

    ineq_linelimit = @(x)create_ineq_linelimit(x(entries{19}),x(entries{20}),x(entries{21}),x(entries{22}),x(entries{23}),x(entries{24}),...
        x(entries{27}(1:Nbranch)),x(entries{27}(Nbranch+1:2*Nbranch)),x(entries{27}(2*Nbranch+1:3*Nbranch)),x(entries{27}(3*Nbranch+1:4*Nbranch)),...
        M_line_plus,M_line_minus,Nbus);
    
    ineq_genboundes = @(x)create_ineq_genbounds(x(entries{19}),x(entries{20}),x(entries{21}),x(entries{22}),x(entries{23}),x(entries{24}),...
        x(entries{29}(1:Nload)),x(entries{29}(Nload+1:2*Nload)),x(entries{30}(2)),x(entries{30}(3)),x(entries{30}(1)),...
        M_ineq_plus,M_ineq_minus,Cg,Cl,pg_max,qg_max,pg_min,qg_min,pl0,ql0,zeta);
    g = @(x)vertcat(eq_cons(x),eq_target(x),ineq_vv(x),ineq_cossin(x),ineq_vvcos_sin(x),ineq_gvv(x),ineq_Psivv(x),ineq_inj(x),ineq_residues(x),ineq_linelimit(x),ineq_genboundes(x));
    %% bound on equality/inequality constraints
    lbg = vertcat(zeros(4*Ngen+2*Nbus+2*Nbranch+2*Nload,1),-inf(12*Nbus+84*Nbranch+8*Ngen+2*size(idx_nslack,1)+4,1));
    ubg = vertcat(zeros(4*Ngen+2*Nbus+2*Nbranch+2*Nload,1),zeros(12*Nbus+84*Nbranch+8*Ngen+2*size(idx_nslack,1)+4,1));

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
    x_SX       =   SX.sym('x',Nx,1);
    constraint = g(x_SX);

    %% solving problem
    % objective 
    obj = @(x) x(entries{30}(5));

    objective = -obj(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    xopt= full(sol.x);
    status = S.stats().return_status;
    if strcmp(status, 'Solve_Succeeded')
        margin_opt = xopt(end);
    else
        margin_opt = -1;
    end
end 

