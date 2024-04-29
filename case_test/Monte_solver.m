function [PCC, time] = Monte_Carlo(mpc)
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
Max_iter = 5000;

bus = mpc.bus;
gen = mpc.gen;

id_gen      = mpc.gen(:,GEN_BUS);
id_slack    =  find(mpc.bus(:,BUS_TYPE) == REF);
id_gen_slack  = find(id_gen == id_slack);
id_gen_nslack = find(id_gen ~= id_slack);
Ngen_nslack = numel(id_gen_nslack);
Nbus        = size(mpc.bus,1);
Ngen        = numel(id_gen);
Nbranch     = size(mpc.branch,1);
vmin = bus(id_gen_nslack,VMIN);
vmax = bus(id_gen_nslack,VMAX);
Pgmin = gen(2:end,PMIN);
Pgmax = gen(2:end,PMAX);

PCC = [];
for i = 1: Max_iter
    % voltage setup
    mpc_rand = mpc;
    mpc_rand.gen(2:end,6) = vmin + (vmax-vmin).*rand(Ngen_nslack,1);
    % active power setup
    baseMVA = mpc.baseMVA;
    mpc_rand.gen(2:end,PG) = (Pgmin + (Pgmax-Pgmin).*rand(Ngen-1,1));

    test_rand = Monte_runpf(mpc_rand);
    if test_rand.success == 1
        if Monte_checkviolation(test_rand)
            PCC =  [PCC; test_rand.gen(1,2:3)/baseMVA];
        end
    end
end

X = PCC(:,1);
Y = PCC(:,2);

k = convhull(X,Y);

convexX = X(k);
convexY = Y(k);

PCC = [convexX, convexY];


time = toc;