function [mpc, sanity_check] = cvxrs(mpc,option,target_mpc,phase_shift)
    if nargin < 3
        target_mpc = [];
        phase_shift = false;
    elseif nargin < 4
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

    Nbus        = size(mpc.bus,1); 
    id_gen      = mpc.gen(:,GEN_BUS);
    Ngen        = numel(id_gen);                 
    Nbranch     = size(mpc.branch,1);
    % added idx
    ptc_factor = 26;
    %don't use the number of load
    
    bus = mpc.bus;
    slack_bus = find(bus(:,BUS_TYPE) == 3);
    idx_nslack = bus(bus(:,BUS_TYPE) ~= 3,BUS_I);
    branch = mpc.branch;
    idx_fr = branch(:,F_BUS);
    idx_to = branch(:,T_BUS);
    gen = mpc.gen;

    %load not used
    idx_load = find((mpc.bus(:,PD).^2+mpc.bus(:,QD).^2) ~= 0);
    idx_pq = find(mpc.bus(:,BUS_TYPE)==1); 
    Npq = sum(mpc.bus(:,BUS_TYPE)==1);
    Nload = numel(idx_load);
    onoff_pq = zeros(Nbus,1); onoff_pq(idx_pq) = 1;

    gencost = mpc.gencost(:,COST:end);
    
    E_fr = sparse(idx_fr, 1:Nbranch, 1, Nbus, Nbranch);
    E_to = sparse(idx_to, 1:Nbranch, 1, Nbus, Nbranch);
    E = E_fr - E_to;
    Cg = sparse(id_gen,1:Ngen,1,Nbus,Ngen);
    Cl = sparse(idx_load,1:Nload,1,Nbus,Nload);

    v0 = bus(:,VM);
    theta0 = deg2rad(bus(:,VA));
    Phi0 = E'*theta0;
    alpha0 = gen(:,ptc_factor);
    delta0 = mpc.delta;

    pg0 = (gen(:,GEN_STATUS).*gen(:,PG))/mpc.baseMVA;
    qg0 = (gen(:,GEN_STATUS).*gen(:,QG))/mpc.baseMVA;
    pl0 = bus(idx_load,PD)/mpc.baseMVA;
    ql0 = bus(idx_load,QD)/mpc.baseMVA;

    %get the injection at each bus
    p_inj0 = Cg*(pg0+alpha0*delta0)-Cl*pl0;
    q_inj0 = Cg*qg0-Cl*ql0;
    power_factor = ql0./(pl0+1e-3);

    vmax = bus(:,VMAX);
    vmin = bus(:,VMIN);
    Phi_max = deg2rad(branch(:,ANGMAX));
    Phi_min = deg2rad(branch(:,ANGMIN));
    pg_max = (gen(:,GEN_STATUS).*gen(:,PMAX))/mpc.baseMVA;
    pg_min = (gen(:,GEN_STATUS).*gen(:,PMIN))/mpc.baseMVA;
    qg_max = (gen(:,GEN_STATUS).*gen(:,QMAX))/mpc.baseMVA;
    qg_min = (gen(:,GEN_STATUS).*gen(:,QMIN))/mpc.baseMVA;
    s_line_max = branch(:,RATE_A)/mpc.baseMVA;

    Sigma0 = mpc.uncertainty.Sigma0;
    gamma0 = mpc.uncertainty.gamma0;

    J1 = -diag(v0(idx_fr)).*v0(idx_to).*sin(E'*theta0)*E';
    J2 = diag(v0(idx_to).*cos(E'*theta0))*E_fr'+diag(v0(idx_fr).*cos(E'*theta0))*E_to';
    J3 = diag(v0(idx_fr).*v0(idx_to).*cos(E'*theta0))*E';
    J4 = diag(v0(idx_to).*sin(E'*theta0))*E_fr'+diag(v0(idx_fr).*sin(E'*theta0))*E_to';
    D = diag(v0);
    J_psi0 = [zeros(2*Nbus,size(idx_nslack,1)+Npq), [Cg * alpha0; zeros(Nbus,1)];
        J1(:,idx_nslack), J2(:,idx_pq), zeros(Nbranch,1);
        J3(:,idx_nslack), J4(:,idx_pq), zeros(Nbranch,1);
        zeros(Nbus,size(idx_nslack,1)), 2*D(:,idx_pq), zeros(Nbus,1)];
    
    g_vvcos0 = v0(idx_fr).*v0(idx_to).*cos(Phi0) - onoff_pq(idx_fr).*v0(idx_fr).*v0(idx_to).*cos(Phi0)...
        - v0(idx_fr).*onoff_pq(idx_to).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi0;
    g_vvsin0 = v0(idx_fr).*v0(idx_to).*sin(Phi0) - onoff_pq(idx_fr).*v0(idx_fr).*v0(idx_to).*sin(Phi0)...
        - v0(idx_fr).*onoff_pq(idx_to).*v0(idx_to).*sin(Phi0) - v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi0;
    g_vv0 = v0.^2-2*onoff_pq.*v0.*v0;
    Psi_vvcos0 = v0(idx_fr).*v0(idx_to).*cos(Phi0);
    Psi_vvsin0 = v0(idx_fr).*v0(idx_to).*sin(Phi0);
    Psi_vv0 = v0.^2;
    
    %[Y,Yf,Yt] = makeYbus_cvxr(mpc,phase_shift);
    [Y,Yf,Yt] = makeYbus(mpc,phase_shift);
    % Auxiliary variables
    y_line = branch(:,BR_STATUS)./(branch(:,BR_R)+1i*branch(:,BR_X));
    b_c = branch(:,BR_STATUS).*branch(:,BR_B);
    tap_ratio = zeros(Nbranch,1);
    for i = 1: Nbranch
        if branch(i,TAP) == 0
            tap_ratio(i) = 1;
        else
            tap_ratio(i) = branch(i,TAP);
        end
    end
    tap = tap_ratio.*exp(1i*branch(:,SHIFT));
    y_sh = bus(:,GS)+bus(:,BS);
    % Csh = not necessary
    y_tt = y_line+1i*b_c/2;
    y_ff = y_tt./(tap.*conj(tap));
    y_ft = -y_line./conj(tap);
    y_tf = -y_line./tap;

    if phase_shift == false
        Y_cos = sparse([idx_fr; idx_to],[(1:Nbranch)'; (1:Nbranch)'],[y_ft;y_tf],Nbus,Nbranch);
        Y_sin = sparse([idx_fr; idx_to],[(1:Nbranch)'; (1:Nbranch)'],[y_ft;-y_tf],Nbus,Nbranch);
        Y_diag = E_to*diag(y_tt)*E_to'+E_fr*diag(y_ff)*E_fr'+diag(y_sh);

        M = [eye(Nbus), zeros(Nbus), -real(Y_cos), -imag(Y_sin), -real(Y_diag);...
         zeros(Nbus), eye(Nbus),  imag(Y_cos), -real(Y_sin),  imag(Y_diag)];
        M_line=[zeros(Nbranch,2*Nbus)  real(diag(y_ft)) imag(diag(y_ft))  real(diag(y_ff)*E_fr');
                zeros(Nbranch,2*Nbus)  real(diag(y_tf)) -imag(diag(y_tf))   real(diag(y_tt)*E_to');
                zeros(Nbranch,2*Nbus) -imag(diag(y_ft)) real(diag(y_ft)) -imag(diag(y_ff)*E_fr');
                zeros(Nbranch,2*Nbus) -imag(diag(y_tf)) -real(diag(y_tf))  -imag(diag(y_tt)*E_to')];
    elseif phase_shift==true
        Y_plus = E_fr*diag(y_ft.*exp.(-1i*Phi0))+E_to*diag(y_tf.*exp.(1i*Phi0));
        Y_minus = E_fr*diag(y_ft.*exp.(-1i*Phi0))-E_to*diag(y_tf.*exp.(1i*Phi0));
        Y_diag = E_to*diag(y_tt)*E_to'+E_fr*diag(y_ff)*E_fr'+diag(y_sh);

        M = [eye(Nbus), zeros(Nbus), -real(Y_plus), -imag(Y_minus), -real(Y_diag);...
        zeros(Nbus), eye(Nbus),  imag(Y_plus), -real(Y_minus),  imag(Y_diag)];
        M_line = [zeros(Nbranch, 2*Nbus),  real(diag(y_ft.*exp(-1i*Phi0))), imag(diag(y_ft.*exp(-1i*Phi0))),  real(diag(y_ff)*E_fr');
              zeros(Nbranch, 2*Nbus),  real(diag(y_tf.*exp(1i*Phi0))),  -imag(diag(y_tf.*exp(1i*Phi0))),   real(diag(y_tt)*E_to');
              zeros(Nbranch, 2*Nbus), -imag(diag(y_ft.*exp(-1i*Phi0))), real(diag(y_ft.*exp(-1i*Phi0))), -imag(diag(y_ff)*E_fr');
              zeros(Nbranch, 2*Nbus), -imag(diag(y_tf.*exp(1i*Phi0))),  -real(diag(y_tf.*exp(1i*Phi0))),  -imag(diag(y_tt)*E_to')];       
    end 

    M_eq = M([1:Nbus, Nbus+idx_pq'],:);
    idx_pvs=setdiff(1:Nbus,idx_pq');
    M_ineq = M(Nbus+idx_pvs,:);
    Nnpq = sum((bus(:, BUS_TYPE) ~= 1),1);

    J_inv_M = -(M_eq * J_psi0) \ M_eq;
    C = [E(idx_nslack,:)' zeros(Nbranch,Npq+1);
        zeros(Npq,size(idx_nslack,1)) eye(Npq) zeros(Npq,1);
        zeros(1,size(idx_nslack,1)+Npq) 1];
    A = [eye(Nbranch+Npq+1); -eye(Nbranch+Npq+1)]*C;
    K = A*J_inv_M;

    theta_test = theta0; v_test = v0; delta_test = delta0;
    theta_test(idx_nslack) = rand(size(idx_nslack,1),1); 
    v_test(idx_pq) = rand(Npq,1);
    v_cpx = v_test.*cos(theta_test) + 1i * v_test .* sin(theta_test);
    p_inj_test = real(v_cpx .* conj(Y * v_cpx));
    q_inj_test = imag(v_cpx .* conj(Y * v_cpx));
 
    if max(abs(M*Psi(mpc,v0,theta0,p_inj0,q_inj0))) > 1e-5
        error("Power flow equation is incorrect (nominal point)");
    end
    if max(abs(M*Psi(mpc,v_test,theta_test,p_inj_test,q_inj_test))) > 1e-5
        warning("Power flow equation is Incorrect (test)"); 
    end
    if max(abs(J_inv_M*gg(mpc,v0,v0,theta0,Phi0,p_inj0,q_inj0)-[theta0(idx_nslack);v0(idx_pq);delta0]))>1e-5
        error("Function g is Incorrect"); 
    end
    if max(abs(J_inv_M*gg(mpc,v0,v0,theta0,Phi0,p_inj0,-Cl*ql0)-[theta0(idx_nslack);v0(idx_pq);delta0]))>1e-5 
        error("Function g without q_pv is Incorrect"); 
    end
    if max(abs(J_inv_M*gg(mpc,v_test,v0,theta_test,Phi0,p_inj_test,q_inj_test)-[theta_test(idx_nslack);v_test(idx_pq);delta_test]))>1e-5 
        warning("Function g is Incorrect (test)"); 
    end

    tol1 = 0;
    K_plus  = max(K,zeros(size(K))).*(K>tol1);
    K_minus = min(K,zeros(size(K))).*(K<-tol1);

    M_ineq_plus  = max(M_ineq, zeros(size(M_ineq))).*(M_ineq>tol1);
    M_ineq_minus = min(M_ineq, zeros(size(M_ineq))).*(M_ineq<-tol1);
    M_line_plus  = max(M_line, zeros(size(M_line))).*(M_line>tol1);
    M_line_minus = min(M_line, zeros(size(M_line))).*(M_line<-tol1);
    
%% set of all variables
    import casadi.*

    v_u = MX.sym('v_u', Nbus, 1);
    v_l = MX.sym('v_l', Nbus, 1);
    Delta_v_u = MX.sym('Delta_v_u', Nbus, 1);
    Delta_v_l = MX.sym('Delta_v_l', Nbus, 1);
    Phi_u = MX.sym('Phi_u', Nbranch, 1);
    Phi_l = MX.sym('Phi_l', Nbranch, 1);
    Delta_Phi_u = MX.sym('Delta_Phi_u', Nbranch, 1);
    Delta_Phi_l = MX.sym('Delta_Phi_l', Nbranch, 1);
    
    delta_u = MX.sym('delta_u');
    delta_l = MX.sym('delta_l');
    Ddelta_u = MX.sym('Ddelta_u');
    Ddelta_l = MX.sym('Ddelta_l');
    
    g_u_pinj = MX.sym('g_u_pinj', Nbus, 1);
    g_l_pinj = MX.sym('g_l_pinj', Nbus, 1);
    g_u_qinj = MX.sym('g_u_qinj', Nbus, 1);
    g_l_qinj = MX.sym('g_l_qinj', Nbus, 1);
    g_u_vvsin = MX.sym('g_u_vvsin', Nbranch, 1);
    g_l_vvsin = MX.sym('g_l_vvsin', Nbranch, 1);
    g_u_vvcos = MX.sym('g_u_vvcos', Nbranch, 1);
    g_l_vvcos = MX.sym('g_l_vvcos', Nbranch, 1);
    g_u_vv = MX.sym('g_u_vv', Nbus, 1);
    g_l_vv = MX.sym('g_l_vv', Nbus, 1);
    
    Psi_u_vvcos = MX.sym('Psi_u_vvcos', Nbranch, 1);
    Psi_l_vvcos = MX.sym('Psi_l_vvcos', Nbranch, 1);
    Psi_u_vvsin = MX.sym('Psi_u_vvsin', Nbranch, 1);
    Psi_l_vvsin = MX.sym('Psi_l_vvsin', Nbranch, 1);
    Psi_u_vv = MX.sym('Psi_u_vv', Nbus, 1);
    Psi_l_vv = MX.sym('Psi_l_vv', Nbus, 1);
    
    Delta_vv_u_uu = MX.sym('Delta_vv_u_uu', Nbranch, 1);
    Delta_vv_u_ll = MX.sym('Delta_vv_u_ll', Nbranch, 1);
    Delta_vv_u_ul = MX.sym('Delta_vv_u_ul', Nbranch, 1);
    Delta_vv_u_lu = MX.sym('Delta_vv_u_lu', Nbranch, 1);
    Delta_vv_l_uu = MX.sym('Delta_vv_l_uu', Nbranch, 1);
    Delta_vv_l_ll = MX.sym('Delta_vv_l_ll', Nbranch, 1);
    Delta_vv_l_ul = MX.sym('Delta_vv_l_ul', Nbranch, 1);
    Delta_vv_l_lu = MX.sym('Delta_vv_l_lu', Nbranch, 1);
    
    Delta_cos_u_u = MX.sym('Delta_cos_u_u', Nbranch, 1);
    Delta_cos_u_l = MX.sym('Delta_cos_u_l', Nbranch, 1);
    Delta_cos_l_u = MX.sym('Delta_cos_l_u', Nbranch, 1);
    Delta_cos_l_l = MX.sym('Delta_cos_l_l', Nbranch, 1);
    Delta_sin_u_u = MX.sym('Delta_sin_u_u', Nbranch, 1);
    Delta_sin_u_l = MX.sym('Delta_sin_u_l', Nbranch, 1);
    Delta_sin_l_u = MX.sym('Delta_sin_l_u', Nbranch, 1);
    Delta_sin_l_l = MX.sym('Delta_sin_l_l', Nbranch, 1);
    
    p_line_fr_u = MX.sym('p_line_fr_u', Nbranch, 1);
    q_line_fr_u = MX.sym('q_line_fr_u', Nbranch, 1);
    p_line_to_u = MX.sym('p_line_to_u', Nbranch, 1);
    q_line_to_u = MX.sym('q_line_to_u', Nbranch, 1);
    
    pg_opt = MX.sym('pg_opt', Ngen, 1);
    pl_opt = MX.sym('pl_opt', Nload, 1);
    ql_opt = MX.sym('ql_opt', Nload, 1);
    pg_opt_u = MX.sym('pg_opt_u', Ngen, 1);
    alpha_opt = MX.sym('alpha_opt', Ngen, 1);
    Dalpha_opt = MX.sym('Dalpha_opt', Ngen, 1);
    
    cost_opt = MX.sym('cost_opt');
    gamma_opt = MX.sym('gamma_opt');
    slack_pg = MX.sym('slack_pg');
    slack_qg = MX.sym('slack_qg');
    slack_line = MX.sym('slack_line');
    margin_opt = MX.sym('margin_opt');
    x = vertcat(v_u,v_l,Delta_v_u,Delta_v_l,Phi_u,Phi_l,Delta_Phi_u,Delta_Phi_l,delta_u,delta_l,Ddelta_u,Ddelta_l,...%8*9+4 = 76:1-76
        g_u_pinj,g_l_pinj,g_u_qinj,g_l_qinj,g_u_vvsin,g_l_vvsin,g_u_vvcos,g_l_vvcos,g_u_vv,g_l_vv,...%90: 77-166
        Psi_u_vvcos,Psi_l_vvcos,Psi_u_vvsin,Psi_l_vvsin,Psi_u_vv,Psi_l_vv,...%54: 167-220
        Delta_vv_u_uu,Delta_vv_u_ll,Delta_vv_u_ul,Delta_vv_u_lu,Delta_vv_l_uu,Delta_vv_l_ll,Delta_vv_l_ul,Delta_vv_l_lu,...%72 221-292
        Delta_cos_u_u,Delta_cos_u_l,Delta_cos_l_u,Delta_cos_l_l,Delta_sin_u_u,Delta_sin_u_l,Delta_sin_l_u,Delta_sin_l_l,...%72 293-364
        p_line_fr_u,q_line_fr_u,p_line_to_u,q_line_to_u,pg_opt,pl_opt,ql_opt,pg_opt_u,alpha_opt,Dalpha_opt,...%60 365-424
        cost_opt,gamma_opt,slack_pg,slack_qg,slack_line,margin_opt);% 6: 425-430
    %% Describe the problem
    % Define initial value
    x0 = vertcat(v0,v0,zeros(2*Nbus,1),Phi0,Phi0,zeros(2*Nbranch,1),delta0,delta0,zeros(2+8*Nbus+28*Nbranch,1),...
        pg0,pl0,ql0,pg0,alpha0,zeros(Ngen,1),mpc.cost,zeros(5,1));
    % Define the lower bound of variables
%     lbx = vertcat(vmin,vmin,vmin-v0,vmin-v0,Phi_min,Phi_min,Phi_min-Phi0,Phi_min-Phi0,...
%         -inf(10+8*Nbus+28*Nbranch+4*Ngen+2*Nload,1));
    lbx = vertcat(vmin,vmin,vmin-v0,vmin-v0,Phi_min,Phi_min,Phi_min-Phi0,Phi_min-Phi0,...
        -inf(4+8*Nbus+28*Nbranch,1),pg_min,-inf(2*Nload+3*Ngen+6,1));

    % Define the upper bound of variables
%     ubx = vertcat(vmax,vmax,vmax-v0,vmax-v0,Phi_max,Phi_max,Phi_max-Phi0,Phi_max-Phi0,...
%         inf(10+8*Nbus+28*Nbranch+4*Ngen+2*Nload,1));
    ubx = vertcat(vmax,vmax,vmax-v0,vmax-v0,Phi_max,Phi_max,Phi_max-Phi0,Phi_max-Phi0,...
        inf(4+8*Nbus+28*Nbranch,1),pg_max,inf(2*Nload+3*Ngen+6,1));
    % Define gfun
    % Equality constraints
    g_v1 = v_u(id_gen)-v_l(id_gen);
    g_v2 = v_u-Delta_v_u;
    g_v3 = v_l-Delta_v_l;
    g_Phi1 = Phi_u-Delta_Phi_u;
    g_Phi2 = Phi_l-Delta_Phi_l;
    g_delta1 = delta_u-Ddelta_u;
    g_delta2 = delta_l-Ddelta_l;
    g_alpha = alpha_opt-Dalpha_opt;
    g_alpha2 = alpha_opt; % ==alpha0
    % inequality constraints 
    % >=
    g_vv1 = Delta_vv_u_uu-(Delta_v_u(idx_fr).*v0(idx_to)+v0(idx_fr).*Delta_v_u(idx_to)+1/4*(Delta_v_u(idx_fr)+Delta_v_u(idx_to)).^2);
    g_vv2 = Delta_vv_u_ll-(Delta_v_l(idx_fr).*v0(idx_to)+v0(idx_fr).*Delta_v_l(idx_to)+1/4*(Delta_v_l(idx_fr)+Delta_v_l(idx_to)).^2);
    g_vv3 = Delta_vv_u_ul-(Delta_v_u(idx_fr).*v0(idx_to)+v0(idx_fr).*Delta_v_l(idx_to)+1/4*(Delta_v_u(idx_fr)+Delta_v_l(idx_to)).^2);
    g_vv4 = Delta_vv_u_lu-(Delta_v_l(idx_fr).*v0(idx_to)+v0(idx_fr).*Delta_v_u(idx_to)+1/4*(Delta_v_l(idx_fr)+Delta_v_u(idx_to)).^2);
    % <= 
    g_vv5 = Delta_vv_l_uu-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)+1/4*(Delta_v_u(idx_fr)-Delta_v_u(idx_to)).^2;
    g_vv6 = Delta_vv_l_ll-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)+1/4*(Delta_v_l(idx_fr)-Delta_v_l(idx_to)).^2;
    g_vv7 = Delta_vv_l_ul-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)+1/4*(Delta_v_u(idx_fr)-Delta_v_l(idx_to)).^2;
    g_vv8 = Delta_vv_l_lu-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)+1/4*(Delta_v_l(idx_fr)-Delta_v_u(idx_to)).^2;
    % >= / <=
    g_cos1 = Delta_cos_u_u+sin(Phi0).*Delta_Phi_u-1/2*(Delta_Phi_u).^2;
    g_cos2 = Delta_cos_u_l+sin(Phi0).*Delta_Phi_l-1/2*(Delta_Phi_l).^2;
    g_cos3 = Delta_cos_l_u+sin(Phi0).*Delta_Phi_u+1/2*(Delta_Phi_u).^2;
    g_cos4 = Delta_cos_l_l+sin(Phi0).*Delta_Phi_l+1/2*(Delta_Phi_l).^2;
    % >= / <=
    g_sin1 = Delta_sin_u_u-cos(Phi0).*Delta_Phi_u-1/2*(Delta_Phi_u).^2;
    g_sin2 = Delta_sin_u_l-cos(Phi0).*Delta_Phi_l-1/2*(Delta_Phi_l).^2;
    g_sin3 = Delta_sin_l_u-cos(Phi0).*Delta_Phi_u+1/2*(Delta_Phi_u).^2;
    g_sin4 = Delta_sin_l_l-cos(Phi0).*Delta_Phi_l+1/2*(Delta_Phi_l).^2;
    % >= ψ_vvcos0
    g_gcos11 = g_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_uu+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
    g_gcos12 = g_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_uu+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
    g_gcos13 = g_u_vvcos - cos(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_ul+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
    g_gcos14 = g_u_vvcos - cos(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_ul+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
    g_gcos15 = g_u_vvcos - cos(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_lu+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
    g_gcos16 = g_u_vvcos - cos(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_lu+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
    g_gcos17 = g_u_vvcos - cos(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_ll+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
    g_gcos18 = g_u_vvcos - cos(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_ll+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
    % <= ψ_vvcos0
    g_gcos21 = g_l_vvcos - cos(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_uu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
    g_gcos22 = g_l_vvcos - cos(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_uu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
    g_gcos23 = g_l_vvcos - cos(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_ul-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
    g_gcos24 = g_l_vvcos - cos(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_ul-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
    g_gcos25 = g_l_vvcos - cos(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_lu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
    g_gcos26 = g_l_vvcos - cos(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_lu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
    g_gcos27 = g_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_ll-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
    g_gcos28 = g_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_ll-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
    % >= Psi_vvsin0
    g_gsin11 = g_u_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_uu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin12 = g_u_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_uu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin13 = g_u_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_ul+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin14 = g_u_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_ul+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin15 = g_u_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_lu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin16 = g_u_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_lu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin17 = g_u_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_ll+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin18 = g_u_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_ll+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    % >= Psi_vvsin0
    g_gsin21 = g_u_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_uu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin22 = g_u_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_uu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin23 = g_u_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_ul+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin24 = g_u_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_ul+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin25 = g_u_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_lu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin26 = g_u_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_lu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin27 = g_u_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_ll+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin28 = g_u_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_ll+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    % <= Psi_vvsin0
    g_gsin31 = g_l_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_uu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin32 = g_l_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_uu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin33 = g_l_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_ul-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin34 = g_l_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_ul-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin35 = g_l_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_lu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin36 = g_l_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_lu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin37 = g_l_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_ll-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin38 = g_l_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_ll-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    % <= Psi_vvsin0
    g_gsin41 = g_l_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_uu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin42 = g_l_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_uu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin43 = g_l_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_ul-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin44 = g_l_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_ul-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin45 = g_l_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_lu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin46 = g_l_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_lu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    g_gsin47 = g_l_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_ll-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
    g_gsin48 = g_l_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_ll-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
    % >= 0
    g_gvv1 = g_u_vv - v_u.^2 + 2*v0.*onoff_pq.*v_u;
    g_gvv2 = g_u_vv - v_l.^2 + 2*v0.*onoff_pq.*v_l;
    % <= Psi_vv0
    g_gvv3 = g_l_vv - 2*v0.*(v_u-v0) + 2*v0.*onoff_pq.*v_u;
    g_gvv4 = g_l_vv - 2*v0.*(v_l-v0) + 2*v0.*onoff_pq.*v_l;
    % >= Psi_vvcos0 / <=
    g_Psi1 = Psi_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_uu+Delta_cos_u_u).^2;
    g_Psi2 = Psi_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_uu+Delta_cos_u_l).^2;
    g_Psi3 = Psi_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_ll-Delta_cos_l_u).^2;
    g_Psi4 = Psi_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_ll-Delta_cos_l_l).^2;
    % >= Psi_vvsin0 / <=
    g_Psi5 = Psi_u_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_uu+Delta_sin_u_u).^2;
    g_Psi6 = Psi_u_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_uu+Delta_sin_u_u).^2;
    g_Psi7 = Psi_l_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_ll-Delta_sin_l_l).^2;
    g_Psi8 = Psi_l_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_uu-Delta_sin_l_l).^2;
    
    g_Psi9 = Psi_u_vv - v_u.^2; %>=0
    g_Psi10 = Psi_u_vv - v_l.^2;%>=0
    % <= Psi_vv0
    g_Psi11 = Psi_l_vv - 2*v0.*Delta_v_u;
    g_Psi12 = Psi_l_vv - 2*v0.*Delta_v_l;
    
    g_g1 = g_u_pinj - Cg*pg_opt + Cl*pl_opt; %>=0
    g_g2 = g_l_pinj - Cg*pg_opt + Cl*pl_opt; %<=0
    g_pgmax = pg_opt + alpha0.*delta_u + slack_pg; %<=pg_max  
    g_pgmin = pg_opt + alpha0.*delta_l - slack_pg; %>=pg_min  
    
    zeta = sqrt(sum((Cl(idx_pvs,:)*diag(power_factor)*chol(Sigma0,'upper')).^2,2));
    % Cl,ql_opt是否替换成非零节点
    g_M1 = M_ineq_minus*[zeros(Nbus,1); Cg*qg_max-Cl*ql_opt;Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]+M_ineq_plus *[zeros(Nbus,1);Cg*qg_max-Cl*ql0;Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv]-slack_qg-zeta.*gamma_opt;%>=0
    g_M2 = M_ineq_plus *[zeros(Nbus,1); Cg*qg_min-Cl*ql_opt;Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]+M_ineq_minus*[zeros(Nbus,1);Cg*qg_min-Cl*ql0;Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv]+slack_qg+zeta.*gamma_opt;%<=0

    g_M3 =  M_line_plus *[zeros(2*Nbus,1);Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]+M_line_minus*[zeros(2*Nbus,1);Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv] - [p_line_fr_u;p_line_to_u;q_line_fr_u;q_line_to_u];%<=0
    g_M4 = -M_line_minus*[zeros(2*Nbus,1);Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]-M_line_plus *[zeros(2*Nbus,1);Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv] - [p_line_fr_u;p_line_to_u;q_line_fr_u;q_line_to_u];%<=0
    g_line1 = p_line_fr_u.^2+q_line_fr_u.^2+slack_line;% <= s_line_max.^2
    g_line2 = p_line_to_u.^2+q_line_to_u.^2+slack_line;% <= s_line_max.^2
    
    xi = sqrt(sum((K(:,1:2*Nbus)*[Cl;Cl*diag(power_factor)]*chol(Sigma0,'upper')).^2,2));
    g_K1 = K_plus*[g_u_pinj;-Cl*ql_opt;g_u_vvcos;g_u_vvsin;g_u_vv]+K_minus*[g_l_pinj;-Cl*ql_opt;g_l_vvcos;g_l_vvsin;g_l_vv]+xi.*gamma_opt-[Phi_u;v_u(idx_pq);delta_u;-Phi_l;-v_l(idx_pq);-delta_l];%<=0
    g_margin1 = slack_pg - margin_opt;
    g_margin2 = slack_qg - margin_opt;
    g_margin3 = slack_line -margin_opt;
    g_margin4 = gamma_opt - margin_opt;%>=0
    g_pgopt = pg_opt+alpha0.*delta_u-pg_opt_u; %<=0 
    g_cost = gencost(:,1)'*(pg_opt_u*mpc.baseMVA).^2+gencost(:,2)'*(pg_opt_u*mpc.baseMVA)+sum(gencost(:,3))-cost_opt;%<=0
    
    gfun = vertcat(g_v1,g_v2,g_v3,g_Phi1,g_Phi2,g_delta1,g_delta2,g_alpha,g_alpha2,...% Ngen+4*Nbus+8
        g_vv1,g_vv2,g_vv3,g_vv4,g_vv5,g_vv6,g_vv7,g_vv8,g_cos1,g_cos2,g_cos3,g_cos4,g_sin1,g_sin2,g_sin3,g_sin4,...% 8*Nbranch+
        g_gcos11,g_gcos12,g_gcos13,g_gcos14,g_gcos15,g_gcos16,g_gcos17,g_gcos18,...
        g_gcos21,g_gcos22,g_gcos23,g_gcos24,g_gcos25,g_gcos26,g_gcos27,g_gcos28,...
        g_gsin11,g_gsin12,g_gsin13,g_gsin14,g_gsin15,g_gsin16,g_gsin17,g_gsin18,...
        g_gsin21,g_gsin22,g_gsin23,g_gsin24,g_gsin25,g_gsin26,g_gsin27,g_gsin28,...
        g_gsin31,g_gsin32,g_gsin33,g_gsin34,g_gsin35,g_gsin36,g_gsin37,g_gsin38,...
        g_gsin41,g_gsin42,g_gsin43,g_gsin44,g_gsin45,g_gsin46,g_gsin47,g_gsin48,...
        g_gvv1,g_gvv2,g_gvv3,g_gvv4,g_Psi1,g_Psi2,g_Psi3,g_Psi4,g_Psi5,g_Psi6,g_Psi7,g_Psi8,g_Psi9,g_Psi10,g_Psi11,g_Psi12,...% 36+108
        g_g1,g_g2,g_pgmax,g_pgmin,g_M1,g_M2,g_M3,g_M4,g_line1,g_line2,g_K1,...%2*Nbus+2*Ngen+6+72+36
        g_margin1,g_margin2,g_margin3,g_margin4,g_pgopt,g_cost);    % 6
    % Define the lower bound of constraints
    lbg = vertcat(zeros(Ngen,1), v0, v0, Phi0, Phi0, delta0, delta0, alpha0, alpha0, zeros(4*Nbranch,1), -inf(4*Nbranch,1), zeros(2*Nbranch,1), -inf(2*Nbranch,1), zeros(2*Nbranch,1), -inf(2*Nbranch,1),...
        repmat(Psi_vvcos0,8,1),-inf(8*length(Psi_vvcos0),1), repmat(Psi_vvsin0,8,1),repmat(Psi_vvsin0,8,1),-inf(8*length(Psi_vvsin0),1),-inf(8*length(Psi_vvsin0),1),...
        zeros(2*Nbus,1), -inf(2*Nbus,1),  Psi_vvcos0,Psi_vvcos0, -inf(2*Nbranch,1),   Psi_vvsin0,Psi_vvsin0, -inf(2*Nbranch,1),      zeros(2*Nbus,1), -inf(2*Nbus,1),  zeros(Nbus,1),-inf(Nbus,1),...
        -inf(Ngen,1),pg_min, zeros(Nnpq,1),-inf(Nnpq,1),-inf(8*Nbranch,1),-inf(2*Nbranch,1),        -inf(2*(Nbranch+Npq+1),1),   zeros(4,1), -inf(Ngen,1),  -inf);
    % Define the upper bound of constraints
    ubg = vertcat(zeros(Ngen,1),v0, v0, Phi0,Phi0,delta0,delta0,alpha0,alpha0,inf(4*Nbranch,1),zeros(4*Nbranch,1), inf(2*Nbranch,1),zeros(2*Nbranch,1),inf(2*Nbranch,1),zeros(2*Nbranch,1),...
        inf(8*length(Psi_vvcos0),1),repmat(Psi_vvcos0,8,1), inf(8*length(Psi_vvsin0),1),inf(8*length(Psi_vvsin0),1),repmat(Psi_vvsin0,8,1),repmat(Psi_vvsin0,8,1),...
        inf(2*Nbus,1),   Psi_vv0,Psi_vv0, inf(2*Nbranch,1),     Psi_vvcos0,Psi_vvcos0, inf(2*Nbranch,1),       Psi_vvsin0, Psi_vvsin0, inf(2*Nbus,1),   Psi_vv0,Psi_vv0, inf(Nbus,1),zeros(Nbus,1),...
        pg_max,inf(Ngen,1), inf(Nnpq,1),zeros(Nnpq,1),zeros(8*Nbranch,1),s_line_max.^2,s_line_max.^2,zeros(2*(Nbranch+Npq+1),1), inf(4,1),   zeros(Ngen,1),   0);

    if ~isempty(target_mpc)
        gen_target = target_mpc.gen;
        pg_target = gen_target(:,GEN_STATUS) .* gen_target(:,PG);
        
        bus_target = target_mpc.bus;
        v_target = bus(:,VM);

        pl_target = bus_target(idx_pq,PD);
        ql_target = bus_target(idx_pq,QD);
    end

    %% solver the problem under different inputs
    tol        = 1e-8;
    options.ipopt.tol             = tol;
    options.ipopt.constr_viol_tol = tol;
    options.ipopt.compl_inf_tol   = tol;
    options.ipopt.acceptable_tol  = tol;
    options.ipopt.acceptable_constr_viol_tol = tol;
    options.ipopt.print_level = 1;


    if option == "margin" && isempty(target_mpc)
        g_add1 = pg_opt;
        g_add2 = pl_opt;
        g_add3 = ql_opt;
        g_add4 = v_u(id_gen);
        g_add5 = margin_opt;
        gfun = vertcat(gfun,g_add1,g_add2,g_add3,g_add4,g_add5);
        lbg = vertcat(lbg,pg0,pl0,ql0,v0(id_gen),0);
        ubg = vertcat(ubg,pg0,pl0,ql0,v0(id_gen),inf);
        ffun = gamma_opt;
        nlp = struct('x',x,'f',-ffun,'g',gfun);
        S = nlpsol('solver','ipopt', nlp, options);
        sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                'lbx', lbx, 'ubx', ubx);
        status = S.stats().return_status;
        if strcmp(status, 'Solve_Succeeded')
            gamma0 = sol.x(end-4);%gamma_opt
            mpc.uncertainty.Sigma0 = Sigma0;
            mpc.uncertainty.gamma0 = gamma0;
        else
            disp('Error in solving the problem.');
        end
    elseif option == "margin" && ~isempty(target_mpc)
        g_add1 = pg_opt;
        g_add2 = pl_opt;
        g_add3 = ql_opt;
        g_add4 = v_u(id_gen);
        ffun = margin_opt;
        gfun = vertcat(gfun,g_add1,g_add2,g_add3,g_add4);
        lbg = vertcat(lbg,pg_target,pl_target,ql_target,v_target);
        ubg = vertcat(ubg,pg_target,pl_target,ql_target,v_target);
        nlp = struct('x',x,'f',-ffun,'g',gfun);
        S = nlpsol('solver','ipopt', nlp, options);
        sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                'lbx', lbx, 'ubx', ubx);
        status = S.stats().return_status;
        if strcmp(status, 'Solve_Succeeded')
            gamma0 = sol.x(end);%margin_opt
            mpc.uncertainty.Sigma0 = Sigma0;
            mpc.uncertainty.gamma0 = gamma0;
        else
            disp('Error in solving the problem.');
        end
    elseif option == "obj"
        g_add1 = pl_opt;
        g_add2 = ql_opt;
        g_add3 = margin_opt;
        g_add4 = gamma_opt;  
        ffun = cost_opt;
        gfun = vertcat(gfun,g_add1,g_add2,g_add3,g_add4);
        lbg = vertcat(lbg,pl0,ql0,0,gamma0);
        ubg = vertcat(ubg,pl0,ql0,inf,inf);
        nlp = struct('x',x,'f',ffun,'g',gfun);
        S = nlpsol('solver','ipopt', nlp, options);  
        sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                'lbx', lbx, 'ubx', ubx);
        status = S.stats().return_status;
        if strcmp(status, 'Solve_Succeeded')
            %mpc.gen(:,PG) = full(sol.x(end-5-4*Ngen-2*Npq:end-5-3*Ngen-2*Npq-1))*mpc.baseMVA; 
            mpc.gen(:,PG) = full(sol.x(end-5-3*Ngen:end-5-2*Ngen-1))*mpc.baseMVA;
            mpc.gen(:,end)= full(sol.x(end-5-2*Ngen:end-5-Ngen-1));%alpha_opt in x
            sol_v_u = full(sol.x(1:Nbus));
            mpc.gen(:,VG) = sol_v_u(id_gen);
            fprintf('%s: obj = %.2f\n', status,full(sol.x(end-5)));
        else
            fprintf('Error in solving the problem: %s\n', status);
        end
    end
    result =struct();
    result.v_u = full(sol.x(1:Nbus));
    result.v_l = full(sol.x(Nbus+1:2*Nbus));
    result.Delta_v_u = full(sol.x(2*Nbus+1:3*Nbus));
    result.Delta_v_l = full(sol.x(3*Nbus+1:4*Nbus));
    result.Phi_u = full(sol.x(4*Nbus+1:4*Nbus+Nbranch));
    result.Phi_l = full(sol.x(4*Nbus+Nbranch+1:4*Nbus+2*Nbranch));
    result.Delta_Phi_u = full(sol.x(4*Nbus+2*Nbranch+1:4*Nbus+3*Nbranch));
    result.Delta_Phi_l = full(sol.x(4*Nbus+3*Nbranch+1:4*Nbus+4*Nbranch));
    result.delta = full(sol.x(4*Nbus+4*Nbranch+1:4*Nbus+4*Nbranch+4));
    result.g_u_pinj = full(sol.x(4*Nbus+4*Nbranch+5:5*Nbus+4*Nbranch+4));
    result.g_l_pinj = full(sol.x(5*Nbus+4*Nbranch+5:6*Nbus+4*Nbranch+4));
    result.g_u_qinj = full(sol.x(6*Nbus+4*Nbranch+5:7*Nbus+4*Nbranch+4));
    result.g_l_qinj = full(sol.x(7*Nbus+4*Nbranch+5:8*Nbus+4*Nbranch+4));
    result.g_u_vvsin = full(sol.x(8*Nbus+4*Nbranch+5:8*Nbus+5*Nbranch+4));
    result.g_l_vvsin = full(sol.x(8*Nbus+5*Nbranch+5:8*Nbus+6*Nbranch+4));
    result.g_u_vvcos = full(sol.x(8*Nbus+6*Nbranch+5:8*Nbus+7*Nbranch+4));
    result.g_l_vvcos = full(sol.x(8*Nbus+7*Nbranch+5:8*Nbus+8*Nbranch+4));
    result.g_u_vv = full(sol.x(8*Nbus+8*Nbranch+5:9*Nbus+8*Nbranch+4));
    result.g_l_vv = full(sol.x(9*Nbus+8*Nbranch+5:10*Nbus+8*Nbranch+4));
    result.Psi_u_vvcos = full(sol.x(10*Nbus+8*Nbranch+5:10*Nbus+9*Nbranch+4));
    result.Psi_l_vvcos = full(sol.x(10*Nbus+9*Nbranch+5:10*Nbus+10*Nbranch+4));
    result.Psi_u_vvsin = full(sol.x(10*Nbus+10*Nbranch+5:10*Nbus+11*Nbranch+4));
    result.Psi_l_vvsin = full(sol.x(10*Nbus+11*Nbranch+5:10*Nbus+12*Nbranch+4));
    result.Psi_u_vv = full(sol.x(10*Nbus+12*Nbranch+5:11*Nbus+12*Nbranch+4));
    result.Psi_l_vv = full(sol.x(11*Nbus+12*Nbranch+5:12*Nbus+12*Nbranch+4));
    result.Delta_vv_u_uu = full(sol.x(12*Nbus+12*Nbranch+5:12*Nbus+13*Nbranch+4));
    result.Delta_vv_u_ll = full(sol.x(12*Nbus+13*Nbranch+5:12*Nbus+14*Nbranch+4));
    result.Delta_vv_u_ul = full(sol.x(12*Nbus+14*Nbranch+5:12*Nbus+15*Nbranch+4));
    result.Delta_vv_u_lu = full(sol.x(12*Nbus+15*Nbranch+5:12*Nbus+16*Nbranch+4));
    result.Delta_vv_l_uu = full(sol.x(12*Nbus+16*Nbranch+5:12*Nbus+17*Nbranch+4));
    result.Delta_vv_l_ll = full(sol.x(12*Nbus+17*Nbranch+5:12*Nbus+18*Nbranch+4));
    result.Delta_vv_l_ul = full(sol.x(12*Nbus+18*Nbranch+5:12*Nbus+19*Nbranch+4));
    result.Delta_vv_l_lu = full(sol.x(12*Nbus+19*Nbranch+5:12*Nbus+20*Nbranch+4));
    result.Delta_cos_u_u = full(sol.x(12*Nbus+20*Nbranch+5:12*Nbus+21*Nbranch+4));
    result.Delta_cos_u_l = full(sol.x(12*Nbus+21*Nbranch+5:12*Nbus+22*Nbranch+4));
    result.Delta_cos_l_u = full(sol.x(12*Nbus+22*Nbranch+5:12*Nbus+23*Nbranch+4));
    result.Delta_cos_l_l = full(sol.x(12*Nbus+23*Nbranch+5:12*Nbus+24*Nbranch+4));
    result.Delta_sin_u_u = full(sol.x(12*Nbus+24*Nbranch+5:12*Nbus+25*Nbranch+4));
    result.Delta_sin_u_l = full(sol.x(12*Nbus+25*Nbranch+5:12*Nbus+26*Nbranch+4));
    result.Delta_sin_l_u = full(sol.x(12*Nbus+26*Nbranch+5:12*Nbus+27*Nbranch+4));
    result.Delta_sin_l_l = full(sol.x(12*Nbus+27*Nbranch+5:12*Nbus+28*Nbranch+4));
    result.p_line_fr_u = full(sol.x(12*Nbus+28*Nbranch+5:12*Nbus+29*Nbranch+4));
    result.q_line_fr_u = full(sol.x(12*Nbus+29*Nbranch+5:12*Nbus+30*Nbranch+4));
    result.p_line_to_u = full(sol.x(12*Nbus+30*Nbranch+5:12*Nbus+31*Nbranch+4));
    result.q_line_to_u = full(sol.x(12*Nbus+31*Nbranch+5:12*Nbus+32*Nbranch+4));
    result.pg_opt = full(sol.x(12*Nbus+32*Nbranch+5:12*Nbus+32*Nbranch+Ngen+4));
    result.pl_opt = full(sol.x(12*Nbus+32*Nbranch+Ngen+5:12*Nbus+32*Nbranch+Ngen+Nload+4));
    result.ql_opt = full(sol.x(12*Nbus+32*Nbranch+Ngen+Nload+5:12*Nbus+32*Nbranch+Ngen+2*Nload+4));
    result.pg_opt_u = full(sol.x(12*Nbus+32*Nbranch+Ngen+2*Nload+5:12*Nbus+32*Nbranch+2*Ngen+2*Nload+4));
    result.alpha_opt = full(sol.x(12*Nbus+32*Nbranch+2*Ngen+2*Nload+5:12*Nbus+32*Nbranch+3*Ngen+2*Nload+4));
    result.Dalpha_opt = full(sol.x(12*Nbus+32*Nbranch+3*Ngen+2*Nload+5:12*Nbus+32*Nbranch+4*Ngen+2*Nload+4));
    result.rest = full(sol.x(12*Nbus+32*Nbranch+4*Ngen+2*Nload+5:12*Nbus+32*Nbranch+4*Ngen+2*Nload+10));

    %% create the structure to store the values of results
    sanity_check = struct();

    sanity_check.v0 = v0;
    sanity_check.Phi0 = Phi0;
    sanity_check.delta0 = delta0;
    sanity_check.alpha0 = alpha0;

    % solution 
    sanity_check.delta_u = result.delta(1);
    sanity_check.delta_l = result.delta(2);
    sanity_check.v_u = result.v_u;
    sanity_check.v_l = result.v_l;
    sanity_check.Phi_u = result.Phi_u;
    sanity_check.Phi_l = result.Phi_l;
    sanity_check.Psi_u_vvcos = result.Psi_u_vvcos;
    sanity_check.Psi_l_vvcos = result.Psi_l_vvcos;
    sanity_check.Psi_u_vvsin = result.Psi_l_vvsin;
    sanity_check.Psi_l_vvsin = result.Psi_l_vvsin;
    sanity_check.Psi_u_vv    = result.Psi_u_vv;
    sanity_check.Psi_l_vv    = result.Psi_l_vv;
    sanity_check.g_u_vvsin = result.g_u_vvsin;
    sanity_check.g_l_vvsin = result.g_l_vvsin;
    sanity_check.g_u_vvcos = result.g_u_vvcos;
    sanity_check.g_l_vvcos = result.g_l_vvcos;
    sanity_check.g_u_vv    = result.g_u_vv;
    sanity_check.g_l_vv    = result.g_l_vv;
    sanity_check.obj = result.rest(1);

    %sanity_check.solve_time = sol.solvertime;
    sanity_check.status = strcmp(status, 'Solve_Succeeded');
end