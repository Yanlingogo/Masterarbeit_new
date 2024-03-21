function [intersections,time] = solver_projection(mpc)
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

    tic;

    baseMVA     = mpc.baseMVA;
    baseKV      = mpc.bus(1,BASE_KV);
    Umax        = mpc.bus(:,VMAX).^2;                                            
    Umin        = mpc.bus(:,VMIN).^2;
    Pgmin       = mpc.gen(:,PMIN)/baseMVA;  
    Qgmin       = mpc.gen(:,QMIN)/baseMVA; 
    Pgmax       = mpc.gen(:,PMAX)/baseMVA; 
    Qgmax       = mpc.gen(:,QMAX)/baseMVA; 
    
    % load data from case file
    Pd          = mpc.bus(:,PD)/baseMVA;
    Qd          = mpc.bus(:,QD)/baseMVA;
    
    id_gen      = mpc.gen(:,GEN_BUS);
    id_slack    =  find(mpc.bus(:,BUS_TYPE) == REF);
    id_gen_slack  = find(id_gen == id_slack);
    id_gen_nslack = find(id_gen ~= id_slack);
    Ngen_nslack = numel(id_gen_nslack);
    Nbus        = size(mpc.bus,1);
    Ngen        = numel(id_gen);
    Nbranch     = size(mpc.branch,1);
    
    branch_r   = mpc.branch(:,BR_R);
    branch_x   = mpc.branch(:,BR_X);
    
    from_bus       = mpc.branch(:, F_BUS);                         
    to_bus         = mpc.branch(:, T_BUS);  
    Cf             = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
    Ct             = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
    C              = Cf - Ct;
    Cg             = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);
    Cg_nslack      = Cg(:,id_gen_nslack);
    Cg_slack       = Cg(:,id_gen_slack);
    
    % beta_V
    A_v = [eye(Nbus) zeros(Nbus,2*(Ngen+Nbranch));
           -eye(Nbus) zeros(Nbus,2*(Ngen+Nbranch))];
    % beta_Generator bus 
    A_gen_ns = [zeros(Ngen,Nbus) eye(Ngen) zeros(Ngen,Ngen+2*(Nbranch));
                zeros(Ngen,Nbus) -eye(Ngen) zeros(Ngen,Ngen+2*(Nbranch));
                zeros(Ngen,(Nbus+Ngen)) eye(Ngen) zeros(Ngen,+2*(Nbranch));
                zeros(Ngen,(Nbus+Ngen)) -eye(Ngen) zeros(Ngen,+2*(Nbranch));];
                % bounds for P Q
    
    % beta_voltage_constraints for injection power
    C_comb = C;
    BR_comb = diag(branch_r);
    BX_comb = diag(branch_x);

    A_inj = [C_comb zeros(Nbranch,2*Ngen) -2*BR_comb -2*BX_comb;
             -C_comb zeros(Nbranch,2*Ngen) 2*BR_comb  2*BX_comb];

    Cg_comb = Cg;
    C_comb2 = C'; % constant matrix for flow variables 
 
    A_eq = [zeros(Nbus) Cg_comb zeros(Nbus, Ngen) -C_comb2 zeros(Nbus, Nbranch) ;
            zeros(Nbus) -Cg_comb zeros(Nbus, Ngen) C_comb2 zeros(Nbus, Nbranch) ;
            zeros(Nbus, (Nbus+Ngen)) Cg_comb zeros(Nbus, Nbranch) -C_comb2 ;
            zeros(Nbus, (Nbus+Ngen)) -Cg_comb zeros(Nbus, Nbranch) C_comb2 ;];

    A = vertcat(A_v,A_gen_ns,A_inj,A_eq);

    tol_cons = 1e-6;
    % lower and upper bounds on U
    b0_U = vertcat(Umax, -Umin);
    
    % lower and upper bounds on P/Q
    b0_P = vertcat(Pgmax,-Pgmin);
    b0_Q = vertcat(Qgmax,-Qgmin);
    
    % power flow equation
    b0_inj = tol_cons*ones(2*Nbranch,1);
    Pd = Pd(:);
    Qd = Qd(:);
    b0_pf = vertcat(Pd+tol_cons,-Pd+tol_cons,Qd+tol_cons,-Qd+tol_cons);
    
    % Integrate all vectors
    
    b0 = vertcat(b0_U,b0_P,b0_Q,b0_inj,b0_pf);

    A_u = A;
    b_u = b0;
%% find the feasible region

    idx_pqs = [1 1+Ngen];
    B = A_u(:,Nbus+idx_pqs);
    remainingIndices = setdiff(1:size(A_u,2), Nbus+idx_pqs);
    A_de = A_u(:, remainingIndices);
    % generate a large enough space
    D0 = [ 1 0;
     -1 0;
      0 1;
      0 -1;];
    v = [Pgmax(id_gen_slack);
         -Pgmin(id_gen_slack);
         Qgmax(id_gen_slack);
         -Qgmin(id_gen_slack)];
    % parameter for M-method
    M = 100;
    [D0,v] = feas_cut(B,A_de,b_u,D0,v,M);
    [D0,v] = remove_redundancy(D0,v);
    intersections = intersection(D0, v);
    time = toc;
end



