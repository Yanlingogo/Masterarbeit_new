clc
clear
close all
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

mpc = case9();
% mpc.bus(:,BUS_I) = [3;5;9;2;4;6;7;1;8];
% mpc.branch(:,F_BUS) = [3;2;4;9;6;7;1;1;8];
% mpc.branch(:,T_BUS) = [2;4;6;6;7;1;5;8;2];
% mpc.gen(:,GEN_BUS) = [3;5;9];


id_bus      = mpc.bus(:,BUS_I);
id_gen      = mpc.gen(:,GEN_BUS);
Nbus        = size(mpc.bus,1);        %size(,1) get the rows of matrix
Ngen        = numel(id_gen);                   %nume1 calculate the sum of the elements
Nbranch     = size(mpc.branch,1);
Nx          = Nbus*2 + Ngen*2; 
Ncore       = numel(id_bus);
sigma       = ones(Nx,1);

%% different aims

user_input = input('Please enter a number (1, or 2): ');
if user_input == 1
    %sequential convex restriction 
    mpc = opf_initialization(mpc,0.05);
    % idx_slack = find(mpc.bus(:,BUS_TYPE)==3);
    % slack_gen = find(mpc.gen(:,GEN_BUS) == idx_slack);
    % mpc.gen(slack_gen,end) = 1;
    mpopt = mpoption('verbose', 0);
    mpc = runpf(mpc,mpopt);
    fprintf('Initial OPF objective solution: %.2f\n', round(mpc.cost, 2));
    max_iter_SCRS = 10;
    [mpc, result_cvxr] = scrs(mpc, max_iter_SCRS);
else 
% Uncertainty
    mpc = opf_initialization(mpc, 0.0);
    gamma_req = 0.2;
    Sigma_0 = mpc.uncertainty.Sigma0;
    gamma_0 = mpc.uncertainty.gamma0;
    mpc.uncertainty.gamma0 = gamma_req;
    max_iter_SCRS = 10;
    [mpc, result_cvxr] = scrs(mpc,max_iter_SCRS);
end
