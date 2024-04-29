function [valid_vertecex] = fun_sample_vertex(mpc)
%% Index setting
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

mpc = ext2int(mpc);

%% parameters
tic

% Time period
T = 1;
% network parameters
baseMVA     = mpc.baseMVA
Pd          = mpc.bus(:,PD)/baseMVA;
Qd          = mpc.bus(:,QD)/baseMVA;

id_gen          = mpc.gen(:,GEN_BUS);
id_slack        =  find(mpc.bus(:,BUS_TYPE) == REF);
id_gen_slack    = find(id_gen == id_slack);
id_gen_nslack   = find(id_gen ~= id_slack);
Ngen_nslack     = numel(id_gen_nslack);
Nbus            = size(mpc.bus,1);
Ngen            = numel(id_gen);
Nbranch         = size(mpc.branch,1);

baseMVA     = mpc.baseMVA;
baseKV      = mpc.bus(1,BASE_KV);
Umax        = mpc.bus(:,VMAX).^2;                                            
Umin        = mpc.bus(:,VMIN).^2;
Pgmin       = mpc.gen(:,PMIN)/baseMVA; Pgmin(id_gen_slack) = -inf;  
Qgmin       = mpc.gen(:,QMIN)/baseMVA; Qgmin(id_gen_slack) = -inf;
Pgmax       = mpc.gen(:,PMAX)/baseMVA; Pgmax(id_gen_slack) = inf;
Qgmax       = mpc.gen(:,QMAX)/baseMVA; Qgmax(id_gen_slack) = inf;

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

Ness = 0;
Nwind = 0;
Npv = 0;

%% Problem formulation
% beta_V
A_v = [eye(Nbus*T) zeros(Nbus*T,2*T*(Ngen+Nbranch+Ness)+T*(Nwind+Npv));
       -eye(Nbus*T) zeros(Nbus*T,2*T*(Ngen+Nbranch+Ness)+T*(Nwind+Npv))];
% beta_Generator bus 
A_gen_ns = [zeros(Ngen*T,Nbus*T) eye(Ngen*T) zeros(Ngen*T,(Nwind+Npv+Ngen)*T+2*T*(Nbranch+Ness));
            zeros(Ngen*T,Nbus*T) -eye(Ngen*T) zeros(Ngen*T,(Nwind+Npv+Ngen)*T+2*T*(Nbranch+Ness));
            zeros(Ngen*T,(Nbus+Ngen)*T) eye(Ngen*T) zeros(Ngen*T,(Nwind+Npv)*T+2*T*(Nbranch+Ness));
            zeros(Ngen*T,(Nbus+Ngen)*T) -eye(Ngen*T) zeros(Ngen*T,(Nwind+Npv)*T+2*T*(Nbranch+Ness));];
            % bounds for P Q

% beta_voltage_constraints for injection power
C_comb = C;
BR_comb = diag(branch_r);
BX_comb = diag(branch_x);
for i = 2:T
    C_comb = blkdiag(C_comb,C);
    BR_comb = blkdiag(BR_comb,diag(branch_r));
    BX_comb = blkdiag(BX_comb,diag(branch_x));
end
A_inj = [C_comb zeros(Nbranch*T,2*T*Ngen) -2*BR_comb -2*BX_comb zeros(Nbranch*T,2*T*Ness+T*(Nwind+Npv));
         -C_comb zeros(Nbranch*T,2*T*Ngen) 2*BR_comb  2*BX_comb zeros(Nbranch*T,2*T*Ness+T*(Nwind+Npv))];

%% beta_power flow balance
Cg_comb = Cg;
C_comb2 = C'; % constant matrix for flow variables 
Cess_comb = [];
Cw_comb = [];
Cp_comb = [];
for i = 2:T
    Cg_comb = blkdiag(Cg_comb,Cg);
    C_comb2 = blkdiag(C_comb2,C');
    Cess_comb = blkdiag(Cess_comb,C_ess);
    Cw_comb = blkdiag(Cw_comb,C_wind);
    Cp_comb = blkdiag(Cp_comb,C_pv);
end
A_eq = [zeros(Nbus*T) Cg_comb zeros(Nbus*T, Ngen*T) -C_comb2 zeros(Nbus*T, Nbranch*T) Cw_comb Cp_comb Cess_comb zeros(Nbus*T,Ness*T);
        zeros(Nbus*T) -Cg_comb zeros(Nbus*T, Ngen*T) C_comb2 zeros(Nbus*T, Nbranch*T) -Cw_comb -Cp_comb -Cess_comb zeros(Nbus*T,Ness*T);
        zeros(Nbus*T, (Nbus+Ngen)*T) Cg_comb zeros(T*Nbus, Nbranch*T) -C_comb2 zeros(T*Nbus, (Nwind+Npv+2*Ness)*T);
        zeros(Nbus*T, (Nbus+Ngen)*T) -Cg_comb zeros(T*Nbus, Nbranch*T) C_comb2 zeros(T*Nbus, (Nwind+Npv+2*Ness)*T);];
% introduce the uncertainty for pg


% Integrate all matrices
%A = vertcat(A_v,A_gen_ns,A_ramp,A_inj,A_ess,A_DERs,A_eq);
A = vertcat(A_v,A_gen_ns,A_inj,A_eq);
%% bounds on constraints
tol_cons = 1e-6;
% lower and upper bounds on U
b0_U = vertcat(repmat(Umax,T,1), -repmat(Umin,T,1));

% lower and upper bounds on P/Q
b0_P = vertcat(repmat(Pgmax,T,1),-repmat(Pgmin,T,1));
b0_Q = vertcat(repmat(Qgmax,T,1),-repmat(Qgmin,T,1));

% power flow equation
b0_inj = tol_cons*ones(2*Nbranch*T,1);
Pd = Pd(:);
Qd = Qd(:);
b0_pf = vertcat(Pd+tol_cons,-Pd+tol_cons,Qd+tol_cons,-Qd+tol_cons);

% Integrate all vectors

b0 = vertcat(b0_U,b0_P,b0_Q,b0_inj,b0_pf);
%% condense model with umbrella constraints
% [A_u, b_u] = E_UCI(A, b0);
A_u = A;
b_u = b0;

idx_pqs = [1 1+Ngen];
C = A_u(:,Nbus+idx_pqs);
remainingIndices = setdiff(1:size(A_u,2), Nbus+idx_pqs);
B = A_u(:, remainingIndices);
%% solver with gurobi
var_y = sdpvar(Nbus+2*Ngen+2*Nbranch-2,1);
var_z = sdpvar(2,1);
% initialization 
cons = B*var_y+C*var_z <= b0;

sigma_0 = [1 0; -1 0; 0 1; 0 -1];
vert = zeros(4,2);
for i = 1:4
    obj = sigma_0(i,:)*var_z;
    options = sdpsettings('solver', 'gurobi');
    sol_1 = optimize(cons, -obj, options);

    vert(i,:) = value(var_z)';
end
vert = unique_coor(vert, 1e-4);
% order the vertices
vert = order_points(vert);
% find the normal vectors
Normals = Norm_vector(vert);

H = 100;
vert_c = sum(vert,1)/size(vert,1);
while any(H >= 1e-4)
    H = zeros(size(Normals,1),1);
    for i = 1:size(Normals,1)
        obj = Normals(i,:)*var_z;
        cons_2 = B*var_y+C*(var_z+vert_c') <= b0;
        sol_2 = optimize(cons_2, -obj, options);
        
        H(i) = Normals(i,:)*(value(var_z)'+vert_c)'-1;
        if  H(i) >= 1e-4 % no movement outwards
            vert(end+1,:) = value(var_z)'+vert_c;
        end
    end
    % update the central point and normal vector
    vert_c = sum(vert,1)/size(vert,1);
    vert = unique_coor(vert, 1e-4);
    vert = order_points(vert);
    Normals_new = Norm_vector(vert);
    Normals = unique_normal(Normals, Normals_new,1e-4);
end
toc
%% delete the reduntant points
valid_vertecex = removeRedun(vert, 1e-6);
plot(valid_vertecex(:,1),valid_vertecex(:,2),'bo');
end