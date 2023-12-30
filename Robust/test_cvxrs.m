function [] = test_cvxrs(mpc,sanity_check)
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
    limit_tol = 1e-5;
    digit_spec = 6;

    if (mpc.delta > sanity_check.delta_u+limit_tol) || (mpc.delta < sanity_check.delta_l-limit_tol)
        format_str = sprintf('delta out of range: %%.%df ≤ %%.%df ≤ %%.%df', digit_spec, digit_spec, digit_spec);
        warningMessage = sprintf(format_str, ...
            round(sanity_check.delta_l, digit_spec), ...
            round(mpc.delta, digit_spec), ...
            round(sanity_check.delta_u, digit_spec));
    end
    bus = mpc.bus;
    for i = 1:size(mpc.bus,1)
        if (bus(i,VM) > sanity_check.v_u(i)+limit_tol) || (bus(i,VM) < sanity_check.v_u(i)+limit_tol)
            warningMessage = sprintf('Vm out of range at bus');
        end

        pq_i = bus(i,BUS_TYPE)==1;
        Psi_vv_i = bus(i,VM).^2;
        g_vv_i = bus(i,VM).^2-2*pq_i.*sanity_check.v0(i).*bus(i,VM);

        if (Psi_vv_i > sanity_check.Psi_u_vv(i)+limit_tol) || (Psi_vv_i<sanity_check.Psi_l_vv(i)-limit_tol)
            warningMessage = sprintf('vv out of range at bus');
        end
        if (g_vv_i > sanity_check.g_u_vv(i)+limit_tol) || (g_vv_i < sanity_check.g_l_vv(i)-limit_tol)
            warningMessage = sprintf('vv out of range at bus');
        end
    end
    
    branch = mpc.branch;
    for i = 1: size(branch,1)
        angle_diff(i) = bus(branch(i,F_BUS),VA) - bus(branch(i,T_BUS),VA);
        angle_diff(i) = deg2rad(angle_diff(i));
        if (angle_diff(i)>sanity_check.Phi_u(i)+limit_tol) || (angle_diff(i)<sanity_check.Phi_l(i)-limit_tol)
            warningMessage = sprintf('Angle out of range at line');
        end

        v_fr = bus(branch(:,F_BUS),VM);
        v_to = bus(branch(:,T_BUS),VM);
        pq_fr = bus(branch(:,F_BUS),BUS_TYPE)==1;
        pq_to = bus(branch(:,F_BUS),BUS_TYPE)==1;
        v0_fr = sanity_check.v0(branch(:,F_BUS));
        v0_to = sanity_check.v0(branch(:,T_BUS));

        Psi_vvcos_i = v_fr(i)*v_to(i)*cos(angle_diff(i));
        Psi_vvsin_i = v_fr(i)*v_to(i)*sin(angle_diff(i));
        g_vvcos_i = v_fr(i)*v_to(i)*cos(angle_diff(i))-pq_fr(i)*v_fr(i)*v0_to(i)...
            *cos(sanity_check.Phi0(i))-v0_fr(i)*pq_to(i)*v_to(i)*cos(sanity_check.Phi0(i))...
            +v0_fr(i)*v0_to(i)*sin(sanity_check.Phi0(i))*angle_diff(i);
        g_vvsin_i = v_fr(i)*v_to(i)*sin(angle_diff(i))-pq_fr(i)*v_fr(i)*v0_to(i)...
            *sin(sanity_check.Phi0(i))-v0_fr(i)*pq_to(i)*v_to(i)*sin(sanity_check.Phi0(i))...
            -v0_fr(i)*v0_to(i)*cos(sanity_check.Phi0(i))*angle_diff(i);

        if (Psi_vvcos_i>sanity_check.Psi_u_vvcos(i)+limit_tol)||(Psi_vvcos_i<sanity_check.Psi_l_vvcos(i)-limit_tol)
            warningMessage = sprintf('vvcos out of range at line');
        end
        if (g_vvcos_i>sanity_check.g_u_vvcos(i)+limit_tol) || (g_vvcos_i<sanity_check.g_l_vvcos(i)-limit_tol)
            warningMessage = sprintf('g_vvcos out of range at line');
        end
        if (Psi_vvsin_i>sanity_check.Psi_u_vvsin(i)+limit_tol)||(Psi_vvsin_i<sanity_check.Psi_l_vvsin(i)-limit_tol)
            warningMessage = sprintf('vvsin out of range at line');
        end
        if (g_vvsin_i>sanity_check.g_u_vvsin(i)+limit_tol) || (g_vvsin_i<sanity_check.g_l_vvsin(i)-limit_tol)
            warningMessage = sprintf('g_vvsin out of range at line');
        end
    end
end

