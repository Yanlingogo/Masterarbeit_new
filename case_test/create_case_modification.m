function [mpc] = create_case_modification(mpc)
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
    Nslack = numel(id_slack);
    id_gen_nslack = find(id_gen ~= id_slack);
    id_gen_slack = find(id_gen == id_slack);

    % bounds modification
    mpc.gen(id_gen_slack,PMAX) = mpc.gen(id_gen_slack,PMAX);
    mpc.gen(id_gen_slack,PMIN) = mpc.gen(id_gen_slack,PMAX);
end
