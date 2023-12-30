function [] = test_runpf_cvxr(mpc)
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

    Nbus = size(mpc.bus,1);
    vm = mpc.bus(:,VM);
    va = mpc.bus(:,VA);
    
    mpopt = mpoption();
    mpopt.out.all = 0; 

    mpc_test = mpc;
    mpc_test.gen(:,PG) = mpc_test.gen(:,PG) + mpc_test.gen(:,end)*mpc_test.delta*mpc.baseMVA;
    mpc_test.bus(:,VA) = rad2deg(mpc_test.bus(:,VA));
    result_test = runpf(mpc_test,mpopt);
    vm_test = result_test.bus(:,VM);
    va_test = deg2rad(result_test.bus(:,VA));

    if max([max(abs(vm - vm_test)), max(abs(va - va_test))]) > 1e-6
    error('PowerFlow Failed');
    end
end

