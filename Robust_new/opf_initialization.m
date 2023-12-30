function [mpc] = opf_initialization(mpc,restriction_level)
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
    
    mpc_opf = mpc;
    id_gen = mpc.gen(:,GEN_BUS);

    if strcmp(restriction_level, 'uniform')  
        mpc_opf.gencost(:,COST:COST+1) = [1.0, 0.0];
        mpc_opf.gen(:,NCOST) = 2;  
        restriction_level = 0;
    elseif restriction_level == -1
        mpc_opf.gencost(:,COST:COST+1) = [1.0, 0.0];
        mpc_opf.gen(:,NCOST) = 2;
    end
    vmax = mpc_opf.bus(:,VMAX);
    vmin = mpc_opf.bus(:,VMIN);
    mpc_opf.bus(:, VMAX) = vmax - restriction_level * (vmax - vmin);
    mpc_opf.bus(:, VMIN) = vmin + restriction_level * (vmax - vmin);

    mpc_opf.gen(:,PMAX) = mpc_opf.gen(:,GEN_STATUS) .* mpc_opf.gen(:,PMAX) - restriction_level *...
        (mpc_opf.gen(:,GEN_STATUS) .* mpc_opf.gen(:,PMAX)-mpc_opf.gen(:,PMIN));
    mpc_opf.gen(:,PMIN) = mpc_opf.gen(:,GEN_STATUS) .* mpc_opf.gen(:,PMIN) + restriction_level *...
        (mpc_opf.gen(:,GEN_STATUS) .* mpc_opf.gen(:,PMAX)-mpc_opf.gen(:,PMIN));
    mpc_opf.gen(:,QMAX) = mpc_opf.gen(:,GEN_STATUS) .* mpc_opf.gen(:,QMAX) - restriction_level *...
        (mpc_opf.gen(:,GEN_STATUS) .* mpc_opf.gen(:,QMAX)-mpc_opf.gen(:,QMIN));
    mpc_opf.gen(:,QMIN) = mpc_opf.gen(:,GEN_STATUS) .* mpc_opf.gen(:,QMIN) + restriction_level *...
        (mpc_opf.gen(:,GEN_STATUS) .* mpc_opf.gen(:,QMAX)-mpc_opf.gen(:,QMIN));

    mpc_opf.branch(:,ANGMAX) = 60;
    mpc_opf.branch(:,ANGMIN) = -60; % operation based on PowerModels in Julia
    mpc_opf.branch(:,ANGMAX) = mpc_opf.branch(:,ANGMAX) - restriction_level * (mpc_opf.branch(:,ANGMAX) - mpc_opf.branch(:, ANGMIN));
    mpc_opf.branch(:,ANGMIN) = mpc_opf.branch(:,ANGMIN) + restriction_level * (mpc_opf.branch(:,ANGMAX) - mpc_opf.branch(:, ANGMIN));
    mpc_opf.branch(:,RATE_A) = mpc_opf.branch(:,RATE_A) * (1-restriction_level);
    mpc_opf.branch(:,RATE_B) = mpc_opf.branch(:,RATE_B) * (1-restriction_level);
    mpc_opf.branch(:,RATE_C) = mpc_opf.branch(:,RATE_C) * (1-restriction_level);
  
    opts = mpoption;
    opts.opf.violation   = 1e-12;
    opts.mips.costtol    = 1e-12;
    opts.mips.gradtol    = 1e-12;
    opts.mips.comptol    = 1e-12;
    opts.opf.ignore_angle_lim = true;
    opts.out.all = 0; 

    %mpopt = mpoption('verbose', 0);
    result_opf = runopf(mpc_opf,opts);

    if result_opf.success ~=1
        error(['OPF Initialization failed. (', result_opf.et, ')']);
    end
    %error string
    %update
    mpc.bus = result_opf.bus;
    mpc.bus(:,VA) = deg2rad(mpc.bus(:,VA));
    mpc.gen = result_opf.gen;
    mpc.branch = result_opf.branch;
    mpc.cost = result_opf.f;

    for i = 1:size(mpc.gen,1)
        if isnan(mpc.gen(i,PG))  % 检查pg字段是否为NaN
        mpc.gen(i,PG) = 0;  % 如果是，将其设置为0
        end
    end
    
    mpc.gen(:,VG) = mpc.bus(id_gen,VM);

    id_load = find((mpc.bus(:,PD).^2+mpc.bus(:,QD).^2) ~= 0);
    pl0 = mpc.bus(id_load,3)/mpc.baseMVA;
    ql0 = mpc.bus(id_load,4)/mpc.baseMVA;
    Sigma0 = diag([pl0.^2;ql0.^2]);
    gamma_0 = 0;
    mpc.uncertainty.Sigma = Sigma0; 
    mpc.uncertainty.gamma0 = gamma_0; 
end