function [pcc_grid_1,pcc_grid_2] = compensation(vert,mpc)
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
    %% Compensation with z_0, second
    % grid the feasible region
    resolution =200;
    pcc_grid = segmentation(vert,resolution);
    
    % parameters
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
    
    B              = (Cf + Ct)';
    A              = [zeros(Nbranch,1) eye(Nbranch)]*B - eye(Nbranch);
    IN              = (eye(Nbranch) - A)\eye(size(A));
    DR             = (eye(Nbranch) - A)\(A*R);
    DX             = (eye(Nbranch) - A)\(A*X);
    DX_p = max(DX,zeros(size(DX))).*(DX>0);
    DX_m = min(DX,zeros(size(DX))).*(DX<0);
    
    Mp = IN'*R*IN;
    Mq = IN'*X*IN;
    
    H_l   = IN'*(2*(R*DR+X*DX) + Z2);
    
    % mpc = runpf(mpc); %pg = 0, qg = 0
    DER_save = mpc.gen(2:end,:);
    mpc.gen(2:end,:) = [];
    mpc = runpf(mpc);
    mpc.gen(2:3,:) = DER_save;
    
    U_0 = mpc.bus(:,8).^2;
    p_pcc0 = mpc.gen(id_gen_slack,PG)/mpc.baseMVA;
    q_pcc0 = mpc.gen(id_gen_slack,QG)/mpc.baseMVA;
    s_0 = [p_pcc0,q_pcc0];
    P = mpc.branch(:,14)/mpc.baseMVA;
    Q = mpc.branch(:,15)/mpc.baseMVA;
    Pg_0 = mpc.gen(gen_nslack,PG)/mpc.baseMVA;
    Qg_0 = mpc.gen(gen_nslack,QG)/mpc.baseMVA;
    L = (P.^2+Q.^2)./U_0(from_bus);
    z_0 = [U_0; P; Q; p_pcc0; q_pcc0; L];
    
    e_st   = sparse(id_slack,1,1,Nbus,1);
    P_p = -Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
    Q_p = -Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
    H =blkdiag(P_p,Q_p); 
    
    D    = full([1, zeros(1,Nbus*2+Nbranch);
            C, -2*diag(branch_r),-2*diag(branch_x), zeros(Nbranch,2);
            zeros(Nbus,Nbus),C',zeros(Nbus,Nbranch),-e_st,zeros(Nbus,1);
            zeros(Nbus,Nbus+Nbranch),C',zeros(Nbus,1),-e_st]);
    RR  = sparse(diag(branch_r));
    XX  = sparse(diag(branch_x));
    PP  = sparse(diag(P));
    QQ  = sparse(diag(Q));
    Ga  = sparse(diag(L.^(-1)));
    
    E   = [sparse(1,Nbranch);RR^2+XX^2;Ct' *RR;Ct'*XX];
    % F   = [Cf, - 2*PP*Ga, -2*QQ*Ga,sparse(Nbranch,2)];
    % G   = (PP^2+QQ^2) * Ga^2;
    F = [L.*Cf,  - 2*PP, - 2*QQ, sparse(Nbranch,2)];
    G = diag(Cf*U_0);

    
    jac_z = [D E;F G]; % jacobian 
    jac_u = [zeros(Nbus,2*(Ngen-1)); 
             -Cg_ns zeros(Nbus,Ngen-1);
             zeros(Nbus,Ngen-1) -Cg_ns;
             zeros(Nbranch,2*(Ngen-1))];
    
    delta_s = pcc_grid - s_0;
    u_0 = [Pg_0;Qg_0];
    e_st     = sparse(id_slack,1,1,Nbus,1);
    pcc_grid_1 = [];
    pcc_grid_2 = [];

    for i = 1:size(delta_s,1)

        p_step_iter = P_p.*delta_s(i,1);
        q_step_iter = Q_p.*delta_s(i,2);


        power_step = [p_step_iter;q_step_iter];
        
        d_z_pred = - jac_z\jac_u*(power_step-u_0);

        z_pred   = z_0 + d_z_pred;
        U_pred   = z_pred(1:Nbus);
        Pij_pred = z_pred(Nbus+1:Nbus+Nbranch);
        Qij_pred = z_pred(Nbus+Nbranch+1:Nbus+2*Nbranch);
        p_pcc    = z_pred(Nbus+2*Nbranch+1); 
        q_pcc    = z_pred(Nbus+2*Nbranch+2);
        L_pred   = z_pred(3*Nbus+1:end);
        e_st     = sparse(id_slack,1,1,Nbus,1);
        g_pred = [e_st'*U_pred-1;
                  C*U_pred - 2*(R*Pij_pred+X*Qij_pred) + Z2*L_pred;
                  C'*Pij_pred + Ct'*R*L_pred - e_st*p_pcc-Cg_ns*p_step_iter+Pd;
                  C'*Qij_pred + Ct'*X*L_pred - e_st*q_pcc-Cg_ns*q_step_iter+Qd;
                    L_pred .* Cf*U_pred - Pij_pred.^2  -  Qij_pred .^2];
                  % L_pred - (Pij_pred.^2+Qij_pred.^2)./U_pred(from_bus,:)];
        d_z_cor = -jac_z\g_pred;
        z_comp = z_0 + d_z_pred + d_z_cor;
    
        U_comp = z_comp(1:Nbus,:);
        Pij_comp = z_comp(Nbus+1:Nbus+Nbranch,:);
        Qij_comp = z_comp(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
        L_comp = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:);
        PQ_loss = ([branch_r';branch_x']*L_comp)';
        PCC = pcc_grid(i,:) + PQ_loss;
    
        % PCC = z_comp(Nbus+2*Nbranch+1:Nbus+2*Nbranch+2)';
        
        pcc_grid_1 = [pcc_grid_1;PCC];

        if all(U_comp >= 0.81 & U_comp <= 1.21)
            pcc_grid_2 = [pcc_grid_2;PCC];
        end
    end

end

