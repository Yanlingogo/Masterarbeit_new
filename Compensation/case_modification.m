function [mpc] = case_modification(mpc, N_der, ratio)
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

%% add the generators
% P_sum  = min(sum(mpc.bus(:,PD)),mpc.gen(1,PMAX))/2;
% Q_sum  = min(sum(mpc.bus(:,QD)),mpc.gen(1,QMAX))/2;

P_sum  = sum(mpc.bus(:,PD))*ratio/N_der;
Q_sum  = sum(mpc.bus(:,QD))*ratio/N_der;

% number of the generators
Nbus = size(mpc.bus,1);
baseMVA = mpc.baseMVA;

mpc.bus([floor(Nbus/3),floor(2*Nbus/3),Nbus-1],BUS_TYPE) = 2;


gen_addition = [
	floor(Nbus/4)	0	0	Q_sum	-Q_sum	1	baseMVA	1	P_sum	-P_sum	0	0	0	0	0	0	0	0	0	0	0;
	floor(2*Nbus/4)	0	0	Q_sum	-Q_sum	1	baseMVA	1	P_sum	-P_sum	0	0	0	0	0	0	0	0	0	0	0;
	floor(3*Nbus/4)	0	0	Q_sum	-Q_sum	1	baseMVA	1	P_sum	-P_sum	0	0	0	0	0	0	0	0	0	0	0;    
    Nbus-1	0	0	Q_sum	-Q_sum  1	baseMVA	1	P_sum	-P_sum	0	0	0	0	0	0	0	0	0	0	0;
    ];
mpc.gen = [mpc.gen;gen_addition];

mpc.gencost =  [
	2	0	0	3	0	20	0;
	2	0	0	3	0	20	0;
	2	0	0	3	0	20	0;
    2	0	0	3	0	20	0;
];


%% parameters modification
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
    mpc.gen(id_gen_slack,PMAX) = 100;
    mpc.gen(id_gen_slack,PMIN) = -100;
    mpc.gen(id_gen_slack,QMAX) = 100;
    mpc.gen(id_gen_slack,QMIN) = -100;

    % voltage at slack bus
    mpc.bus(id_slack,VMAX) = 1;
    mpc.bus(id_slack,VMIN) = 1;

    mpc = rmfield(mpc, 'order');
    mpc.bus(2:end,VMAX) = 1.1*ones(Nbranch,1);
    mpc.bus(2:end,VMIN) = 0.9*ones(Nbranch,1);
end

