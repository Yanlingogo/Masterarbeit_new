function [U1_plot,U2_plot,exact_plot,cvxrs_plot,u_plot0] = plot2D(mpc,plot_bus,plot_rng,option,resolution)
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

    if strcmp(option,'pg')
        u_plot0 = [mpc.gen(plot_bus(1),PG),mpc.gen(plot_bus(2),PG)]/mpc.baseMVA;
    elseif strcmp(option,'pd')
        u_plot0 = [mpc.bus(plot_bus(1),PD),mpc.bus(plot_bus(2),PD)]/mpc.baseMVA;
    end

    u1_plot = linspace(plot_rng(1), plot_rng(2), resolution);
    u2_plot = linspace(plot_rng(3), plot_rng(4), resolution);
    U1_plot = repmat(u1_plot, length(u2_plot), 1);
    U2_plot = repmat(u2_plot', 1, length(u1_plot));
    exact_plot = zeros(size(U1_plot));
    cvxrs_plot = zeros(size(U1_plot));
    Sigma0 = mpc.uncertainty.Sigma0;
    gamma_plot = mpc.uncertainty.gamma0;

    Ngen = size(mpc.gen,1);
    pg_plot = mpc.gen(:,GEN_STATUS).*mpc.gen(:,PG)/mpc.baseMVA;
    for i = 1:size(U1_plot, 1)
        fprintf('progress: %d/%d\n', i, size(U1_plot, 1));
        for j = 1:size(U1_plot, 2)
            mpc_plot = mpc; 
            if strcmp(option, 'pg')
                mpc_plot.gen(plot_bus(1),PG) = U1_plot(i, j);
                mpc_plot.gen(plot_bus(2),PG) = U2_plot(i, j);
            elseif strcmp(option, 'pd')
                mpc_plot.bus(plot_bus(1),PD) = U1_plot(i, j)*mpc.baseMVA;
                mpc_plot.bus(plot_bus(2),PD) = U2_plot(i, j)*mpc.baseMVA;
                pf1 = mpc.bus(plot_bus(1),QD) / mpc.bus(plot_bus(1),PD);
                pf2 = mpc.bus(plot_bus(2),QD) / mpc.bus(plot_bus(2),PD);
                mpc_plot.bus(plot_bus(1),QD) = pf1 * U1_plot(i, j)*mpc.baseMVA;
                mpc_plot.bus(plot_bus(2),QD) = pf2 * U2_plot(i, j)*mpc.baseMVA;
            end
            mpc_plot = runpf_cvxr(mpc_plot);
            % test_runpf(network_data_plot);
%             [violation_status, margin_plot] = check_violation(mpc_plot, 1);
%             exact_plot(i, j) = margin_plot;
            [mpc_plot, sanity_check] = cvxrs2(mpc, 'margin', mpc_plot);
            Sigma0 = mpc_plot.uncertainty.Sigma0;
            gamma_plot = mpc_plot.uncertainty.gamma0;
            cvxrs_plot(i, j) = gamma_plot;
        end
    end

end

