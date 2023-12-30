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

    %ppq0 = bus(idx_pq,PD)/mpc.baseMVA; 
    %qpq0 = bus(idx_pq,QD)/mpc.baseMVA;
    %get the injection at each bus
    p_inj0 = Cg*(pg0+alpha0*delta0)-Cl*pl0;
    q_inj0 = Cg*qg0-Cl*ql0;
    power_factor = ql0./(pl0+1e-3);

    vmax = bus(:,VMAX);
    vmin = bus(:,VMIN);
%     Phi_max = deg2rad(branch(:,ANGMAX));
%     Phi_min = deg2rad(branch(:,ANGMIN));
    Phi_min =  -1.0472*eye(9,1);
    Phi_max =  1.0472*eye(9,1);
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

    zeta = sqrt(sum((Cl(idx_pvs,:)*diag(power_factor)*chol(Sigma0,'upper')).^2,2));
    xi = sqrt(sum((K(:,1:2*Nbus)*[Cl;Cl*diag(power_factor)]*chol(Sigma0,'upper')).^2,2));
    
%% set of all variables 
    v_u = sdpvar(Nbus, 1);
    v_l = sdpvar(Nbus, 1);
    Delta_v_u = sdpvar(Nbus, 1);
    Delta_v_l = sdpvar(Nbus, 1);
    Phi_u = sdpvar(Nbranch, 1);
    Phi_l = sdpvar(Nbranch, 1);
    Delta_Phi_u = sdpvar(Nbranch, 1);
    Delta_Phi_l = sdpvar(Nbranch, 1);

    delta_u = sdpvar(1);
    delta_l = sdpvar(1);
    Ddelta_u = sdpvar(1);
    Ddelta_l = sdpvar(1);
    
    g_u_pinj = sdpvar(Nbus, 1);
    g_l_pinj = sdpvar(Nbus, 1);
    %g_u_qinj = sdpvar(Nbus, 1);
    %g_l_qinj = sdpvar(Nbus, 1);
    g_u_vvsin = sdpvar(Nbranch, 1);
    g_l_vvsin = sdpvar(Nbranch, 1);
    g_u_vvcos = sdpvar(Nbranch, 1);
    g_l_vvcos = sdpvar(Nbranch, 1);
    g_u_vv = sdpvar(Nbus, 1);
    g_l_vv = sdpvar(Nbus, 1);
    
    Psi_u_vvcos = sdpvar(Nbranch, 1);
    Psi_l_vvcos = sdpvar(Nbranch, 1);
    Psi_u_vvsin = sdpvar(Nbranch, 1);
    Psi_l_vvsin = sdpvar(Nbranch, 1);
    Psi_u_vv = sdpvar(Nbus, 1);
    Psi_l_vv = sdpvar(Nbus, 1);

    Delta_vv_u_uu = sdpvar(Nbranch, 1);
    Delta_vv_u_ll = sdpvar(Nbranch, 1);
    Delta_vv_u_ul = sdpvar(Nbranch, 1);
    Delta_vv_u_lu = sdpvar(Nbranch, 1);
    Delta_vv_l_uu = sdpvar(Nbranch, 1);
    Delta_vv_l_ll = sdpvar(Nbranch, 1);
    Delta_vv_l_ul = sdpvar(Nbranch, 1);
    Delta_vv_l_lu = sdpvar(Nbranch, 1);
    
    Delta_cos_u_u = sdpvar(Nbranch, 1);
    Delta_cos_u_l = sdpvar(Nbranch, 1);
    Delta_cos_l_u = sdpvar(Nbranch, 1);
    Delta_cos_l_l = sdpvar(Nbranch, 1);
    Delta_sin_u_u = sdpvar(Nbranch, 1);
    Delta_sin_u_l = sdpvar(Nbranch, 1);
    Delta_sin_l_u = sdpvar(Nbranch, 1);
    Delta_sin_l_l = sdpvar(Nbranch, 1);
  
    p_line_fr_u = sdpvar(Nbranch, 1);
    q_line_fr_u = sdpvar(Nbranch, 1);
    p_line_to_u = sdpvar(Nbranch, 1);
    q_line_to_u = sdpvar(Nbranch, 1);
    
    pg_opt = sdpvar(Ngen, 1);
    pl_opt = sdpvar(Nload, 1);
    ql_opt = sdpvar(Nload, 1);
    pg_opt_u = sdpvar(Ngen, 1);
    alpha_opt = sdpvar(Ngen, 1);
    Dalpha_opt = sdpvar(Ngen, 1);
    
    cost_opt = sdpvar;
    gamma_opt = sdpvar;
    slack_pg = sdpvar;
    slack_qg = sdpvar;
    slack_line = sdpvar;
    margin_opt = sdpvar;

%% constraints
    Constraints = [v_u(id_gen)-v_l(id_gen) == 0, 
                v_u-Delta_v_u == v0,
                v_l-Delta_v_l == v0,
                Phi_u-Delta_Phi_u == Phi0,
                Phi_l-Delta_Phi_l == Phi0,
                delta_u-Ddelta_u == delta0,
                delta_l-Ddelta_l == delta0,
                alpha_opt-Dalpha_opt == alpha0,
                alpha_opt == alpha0
                ];%pg_opt_u <= pg_max,pg_opt_u>=pg_min

    Constraints = [Constraints, 
                    Delta_vv_u_uu-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)-1/4*(Delta_v_u(idx_fr)+Delta_v_u(idx_to)).^2>=0,
                    Delta_vv_u_ll-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)-1/4*(Delta_v_l(idx_fr)+Delta_v_l(idx_to)).^2>=0,
                    Delta_vv_u_ul-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)-1/4*(Delta_v_u(idx_fr)+Delta_v_l(idx_to)).^2>=0,
                    Delta_vv_u_lu-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)-1/4*(Delta_v_l(idx_fr)+Delta_v_u(idx_to)).^2>=0];

    Constraints = [Constraints, 
                    Delta_vv_l_uu-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)+1/4*(Delta_v_u(idx_fr)-Delta_v_u(idx_to)).^2<=0,
                    Delta_vv_l_ll-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)+1/4*(Delta_v_l(idx_fr)-Delta_v_l(idx_to)).^2<=0,
                    Delta_vv_l_ul-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)+1/4*(Delta_v_u(idx_fr)-Delta_v_l(idx_to)).^2<=0,
                    Delta_vv_l_lu-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)+1/4*(Delta_v_l(idx_fr)-Delta_v_u(idx_to)).^2<=0];
    Constraints = [Constraints, 
        Delta_cos_u_u+sin(Phi0).*Delta_Phi_u-1/2*(Delta_Phi_u).^2>=0,
        Delta_cos_u_l+sin(Phi0).*Delta_Phi_l-1/2*(Delta_Phi_l).^2>=0,
        Delta_cos_l_u+sin(Phi0).*Delta_Phi_u+1/2*(Delta_Phi_u).^2<=0,
        Delta_cos_l_l+sin(Phi0).*Delta_Phi_l+1/2*(Delta_Phi_l).^2<=0];
    Constraints = [Constraints, 
        Delta_sin_u_u-cos(Phi0).*Delta_Phi_u-1/2*(Delta_Phi_u).^2>=0,
        Delta_sin_u_l-cos(Phi0).*Delta_Phi_l-1/2*(Delta_Phi_l).^2>=0,
        Delta_sin_l_u-cos(Phi0).*Delta_Phi_u+1/2*(Delta_Phi_u).^2<=0,
        Delta_sin_l_l-cos(Phi0).*Delta_Phi_l+1/2*(Delta_Phi_l).^2<=0];
    Constraints = [Constraints,
        g_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_uu+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u>=Psi_vvcos0,
        g_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_uu+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l>=Psi_vvcos0,
        g_u_vvcos - cos(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_ul+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u>=Psi_vvcos0,
        g_u_vvcos - cos(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_ul+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l>=Psi_vvcos0,
        g_u_vvcos - cos(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_lu+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u>=Psi_vvcos0,
        g_u_vvcos - cos(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_lu+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l>=Psi_vvcos0,
        g_u_vvcos - cos(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_ll+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u>=Psi_vvcos0,
        g_u_vvcos - cos(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_ll+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l>=Psi_vvcos0];
    Constraints = [Constraints,
        g_l_vvcos - cos(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_uu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u<=Psi_vvcos0,
        g_l_vvcos - cos(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_uu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l<=Psi_vvcos0,
        g_l_vvcos - cos(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_uu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u<=Psi_vvcos0,
        g_l_vvcos - cos(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_uu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l<=Psi_vvcos0,
        g_l_vvcos - cos(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_uu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u<=Psi_vvcos0,
        g_l_vvcos - cos(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_uu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l<=Psi_vvcos0,
        g_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_uu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u<=Psi_vvcos0,
        g_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_uu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l<=Psi_vvcos0];
    Constraints = [Constraints,
        g_u_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_uu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_uu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_ul+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_ul+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_lu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_lu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_ll+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_ll+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l>= Psi_vvsin0];
    Constraints = [Constraints,
        g_u_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_uu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_uu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_ul+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_ul+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_lu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_lu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_ll+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u>= Psi_vvsin0,
        g_u_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_ll+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l>= Psi_vvsin0];
    % <= Psi_vvsin0
    Constraints = [Constraints,
        g_l_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_uu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_uu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_ul-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_ul-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_lu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_lu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_ll-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_ll-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l<= Psi_vvsin0];
    % <= Psi_vvsin0
    Constraints = [Constraints,
        g_l_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_uu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_uu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_ul-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_ul-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_lu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_lu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_ll-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u<= Psi_vvsin0,
        g_l_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_ll-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l<= Psi_vvsin0];

    % >= 0
    Constraints = [Constraints,
        g_u_vv - v_u.^2 + 2*v0.*onoff_pq.*v_u>=0,
        g_u_vv - v_l.^2 + 2*v0.*onoff_pq.*v_l>=0];
    % <= Psi_vv0
    Constraints = [Constraints,
        g_l_vv - 2*v0.*(v_u-v0) + 2*v0.*onoff_pq.*v_u<= Psi_vv0,
        g_l_vv - 2*v0.*(v_l-v0) + 2*v0.*onoff_pq.*v_l<= Psi_vv0];

    % >= Psi_vvcos0 / <=
    Constraints = [Constraints,
        Psi_u_vvcos >= cos(Phi0).*Delta_vv_u_uu + v0(idx_fr).*v0(idx_to).*Delta_cos_u_u + 1/4*(Delta_vv_u_uu+Delta_cos_u_u).^2+Psi_vvcos0,
        Psi_u_vvcos >= cos(Phi0).*Delta_vv_u_uu + v0(idx_fr).*v0(idx_to).*Delta_cos_u_l + 1/4*(Delta_vv_u_uu+Delta_cos_u_l).^2+Psi_vvcos0,
        Psi_l_vvcos <= cos(Phi0).*Delta_vv_l_ll + v0(idx_fr).*v0(idx_to).*Delta_cos_l_u - 1/4*(Delta_vv_l_ll-Delta_cos_l_u).^2+Psi_vvcos0,
        Psi_l_vvcos <= cos(Phi0).*Delta_vv_l_ll + v0(idx_fr).*v0(idx_to).*Delta_cos_l_l - 1/4*(Delta_vv_l_ll-Delta_cos_l_l).^2+Psi_vvcos0];
    % >= Psi_vvsin0 / <=
    Constraints = [Constraints,
        Psi_u_vvsin >= sin(Phi0).*Delta_vv_u_uu + v0(idx_fr).*v0(idx_to).*Delta_sin_u_u + 1/4*(Delta_vv_u_uu+Delta_sin_u_u).^2+ Psi_vvsin0,
        Psi_u_vvsin >= sin(Phi0).*Delta_vv_l_uu + v0(idx_fr).*v0(idx_to).*Delta_sin_u_u + 1/4*(Delta_vv_l_uu+Delta_sin_u_u).^2+ Psi_vvsin0,
        Psi_l_vvsin <= sin(Phi0).*Delta_vv_l_ll + v0(idx_fr).*v0(idx_to).*Delta_sin_l_l - 1/4*(Delta_vv_l_ll-Delta_sin_l_l).^2+ Psi_vvsin0,
        Psi_l_vvsin <= sin(Phi0).*Delta_vv_u_uu + v0(idx_fr).*v0(idx_to).*Delta_sin_l_l - 1/4*(Delta_vv_u_uu-Delta_sin_l_l).^2+ Psi_vvsin0];
    %>=0
    Constraints = [Constraints,Psi_u_vv>=v_u.^2, Psi_u_vv>=v_l.^2];
    % <= Psi_vv0
    Constraints = [Constraints,
        Psi_l_vv <= 2*v0.*Delta_v_u+Psi_vv0,
        Psi_l_vv <= 2*v0.*Delta_v_l+Psi_vv0];

    Constraints = [Constraints,
        g_u_pinj >= Cg*pg_opt - Cl*pl_opt, %>=0
        g_l_pinj <= Cg*pg_opt - Cl*pl_opt, %<=0
        pg_opt + alpha0.*delta_u + slack_pg<=pg_max, %<=pg_max  
        pg_opt + alpha0.*delta_l - slack_pg>=pg_min];

    Constraints = [Constraints,M_ineq_minus*[zeros(Nbus,1); Cg*qg_max-Cl*ql_opt;Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]+M_ineq_plus *[zeros(Nbus,1);Cg*qg_max-Cl*ql0;Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv]-slack_qg-zeta.*gamma_opt>=0];%>=0
    Constraints = [Constraints,M_ineq_plus *[zeros(Nbus,1); Cg*qg_min-Cl*ql_opt;Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]+M_ineq_minus*[zeros(Nbus,1);Cg*qg_min-Cl*ql0;Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv]+slack_qg+zeta.*gamma_opt<=0];%<=0
    Constraints = [Constraints, M_line_plus *[zeros(2*Nbus,1);Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]+M_line_minus*[zeros(2*Nbus,1);Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv] <= [p_line_fr_u;p_line_to_u;q_line_fr_u;q_line_to_u]];%<=0
    Constraints = [Constraints,-M_line_minus*[zeros(2*Nbus,1);Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]-M_line_plus *[zeros(2*Nbus,1);Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv] <= [p_line_fr_u;p_line_to_u;q_line_fr_u;q_line_to_u]];%<=0
    
    Constraints = [Constraints,p_line_fr_u.^2+q_line_fr_u.^2+slack_line<= s_line_max.^2];% <= s_line_max.^2
    Constraints = [Constraints,p_line_to_u.^2+q_line_to_u.^2+slack_line<= s_line_max.^2];% <= s_line_max.^2

    Constraints = [Constraints,K_plus*[g_u_pinj;-Cl*ql_opt;g_u_vvcos;g_u_vvsin;g_u_vv]+K_minus*[g_l_pinj;-Cl*ql_opt;g_l_vvcos;g_l_vvsin;g_l_vv]+xi.*gamma_opt<=[Phi_u;v_u(idx_pq);delta_u;-Phi_l;-v_l(idx_pq);-delta_l]];%<=0
    Constraints = [Constraints,
        slack_pg - margin_opt>=0,
        slack_qg - margin_opt>=0,
        slack_line -margin_opt>=0,
        gamma_opt - margin_opt>=0,%>=0
        pg_opt+alpha0.*delta_u-pg_opt_u<=0]; %<=0 
    Constraints = [Constraints,gencost(:,1)'*(pg_opt_u*mpc.baseMVA).^2+gencost(:,2)'*(pg_opt_u*mpc.baseMVA)+sum(gencost(:,3))-cost_opt<=0];%<=0

    if option == "margin" && isempty(target_mpc)
        g_add1 = pg_opt;
        g_add2 = pl_opt;
        g_add3 = ql_opt;
        g_add4 = v_u(id_gen);
        g_add5 = margin_opt;
        gfun = vertcat(gfun,g_add1,g_add2,g_add3,g_add4,g_add5);
        lbg = vertcat(lbg,pg0,ppq0,qpq0,v0(id_gen),0);
        ubg = vertcat(ubg,pg0,ppq0,qpq0,v0(id_gen),inf);
        ffun = gamma_opt;
        nlp = struct('x',x,'f',-ffun,'g',gfun);
        S = nlpsol('solver','ipopt', nlp,options);
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
        S = nlpsol('solver','ipopt', nlp,options);
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
        Constraints = [Constraints,pl_opt == pl0];
        Constraints = [Constraints,ql_opt == ql0];
        Constraints = [Constraints,margin_opt>=0];
        Constraints = [Constraints,gamma_opt>=gamma0];
        Objective = cost_opt;
        options = sdpsettings('solver', 'gurobi');
        result = optimize(Constraints, Objective, options);
        rdelta_u = value(delta_u);
        rpg_opt = value(pg_opt);
        rpg_opt_u = value(pg_opt_u);
    end
end