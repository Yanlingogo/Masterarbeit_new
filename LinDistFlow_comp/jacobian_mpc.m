function [parameters] = jacobian_mpc(mpc)
%% index
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
%% parameters
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
    
    mpc = runpf(mpc); %pg = 0, qg = 0
    
    
    U_0 = mpc.bus(:,8).^2;
    p_pcc0 = mpc.gen(id_gen_slack,PG)/mpc.baseMVA;
    q_pcc0 = mpc.gen(id_gen_slack,QG)/mpc.baseMVA;
    s_0 = [p_pcc0,q_pcc0];
    P = mpc.branch(:,14)/mpc.baseMVA;
    Q = mpc.branch(:,15)/mpc.baseMVA;
    L = (P.^2+Q.^2)./U_0(from_bus);
    z_0 = [U_0; P; Q; p_pcc0; q_pcc0; L];
    
    e_st   = sparse(id_slack,1,1,Nbus,1);
    P_p = Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
    Q_p = Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
    H =blkdiag(P_p,Q_p); 
    
    D    = full([1, zeros(1,Nbus*2+Nbranch);
            C, -2*diag(branch_r),-2*diag(branch_x), zeros(Nbranch,2);
            zeros(Nbus,Nbus),C',zeros(Nbus,Nbranch),-e_st,zeros(Nbus,1);
            zeros(Nbus,Nbus+Nbranch),C',zeros(Nbus,1),-e_st]);

    PP  = sparse(diag(P));
    QQ  = sparse(diag(Q));
    Ga  = sparse(diag(L.^(-1)));
    
    E   = [sparse(1,Nbranch);R^2+X^2;Ct' *R;Ct'*X];
    F   = [Cf, - 2*PP*Ga, -2*QQ*Ga,sparse(Nbranch,2)];
    G   = (PP^2+QQ^2) * Ga^2;
    
    jac_z = [D E;F G]; % jacobian 
    jac_u = [zeros(Nbus,2*(Ngen-1)); 
             -Cg_ns zeros(Nbus,Ngen-1);
             zeros(Nbus,Ngen-1) -Cg_ns;
             zeros(Nbranch,2*(Ngen-1))];

    parameters = struct();
    parameters.pcc_0 = s_0;
    parameters.z_0 = z_0;
    parameters.jac_z = jac_z;
    parameters.jac_u = jac_u;
    parameters.branch_r = branch_r;
    parameters.branch_x = branch_x;
end

