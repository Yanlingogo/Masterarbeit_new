function [violation_status,margin] = check_violation(mpc,suppress_warning)
    if nargin < 2
        suppress_warning = 0;
    else
        suppress_warning = suppress_warning;
    end
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
    Nbranch = size(mpc.branch,1);
    Ngen =size(mpc.gen,1);
    Npq = sum(mpc.bus(:,BUS_TYPE)==1);
    bus = mpc.bus;
    branch = mpc.branch;
    gen = mpc.gen;
    idx_fr = branch(:,F_BUS);
    idx_to = branch(:,T_BUS);
    id_gen = mpc.gen(:,GEN_BUS);
    idx_pq = find(mpc.bus(:,BUS_TYPE)==1);
    idx_load = find((mpc.bus(:,PD).^2+mpc.bus(:,QD).^2) ~= 0);
    Nload = numel(idx_load);

    Cg = sparse(id_gen,1:Ngen,1,Nbus,Ngen);
    Cl = sparse(idx_load,1:Nload,1,Nbus,Nload);

    v0 = bus(:,VM);
    theta0 = bus(:,VA);

    v_cplx0 = v0.*exp(1i*theta0);
    [Y,Yf,Yt] = makeYbus(mpc);
    
    gen(:,PG) = gen(:,PG)/mpc.baseMVA;
    gen(:,PMAX) = gen(:,PMAX)/mpc.baseMVA;
    gen(:,PMIN) = gen(:,PMIN)/mpc.baseMVA;
    ql0 = bus(idx_load,QD)/mpc.baseMVA;
    qg_inj0 = imag(v_cplx0.*conj(Y*v_cplx0))+Cl*ql0;
    s_line_fr0 = v_cplx0(idx_fr).*conj(Yf*v_cplx0);
    s_line_to0 = v_cplx0(idx_to).*conj(Yt*v_cplx0);

    qg_max = (gen(:,GEN_STATUS).*gen(:,QMAX))/mpc.baseMVA;
    qg_min = (gen(:,GEN_STATUS).*gen(:,QMIN))/mpc.baseMVA;

    qinj_max = Cg*qg_max;
    qinj_min = Cg*qg_min;

    limit_tolerence = 1e-5;
    margin = 0.1;
    digits_spec = 2;

    violation_status = struct();
    violation_status.vmag = [];
    violation_status.qg = [];
    violation_status.pg = [];
    violation_status.angle = [];
    violation_status.line = [];


    if isnan(mpc.delta)&&suppress_warning ==0
        error('Base point did not converge!');
    end

    for i = 1:Nbus
        if bus(i,VM)>(bus(i,VMAX)+limit_tolerence) || bus(i,VM)<(bus(i,VMIN)-limit_tolerence)
            violation_status.vmag(end+1) = i;
            if suppress_warning == 0
                warning(['Vm violation at bus ', num2str(i), ': ', num2str(round(bus(i,VMIN), digits_spec)), ' ≤ ', num2str(round(bus(i,VM),digits_spec)), ' ≤ ', num2str(round(bus(i,VMAX), digits_spec))]);
            end
            margin = min(margin, min(bus(i,VMAX) - bus(i,VM), bus(i,VM) - bus(i,VMIN)));
        end
        if (qg_inj0(i)>qinj_max(i)+limit_tolerence || qg_inj0(i)<qinj_min(i)-limit_tolerence)
            violation_status.qg(end+1) = i;
            if suppress_warning == 0
                warning(['Qg violation at bus ', num2str(i), ': ', num2str(round(qinj_min(i), digits_spec)), ' ≤ ', num2str(round(qg_inj0(i), digits_spec)), ' ≤ ', num2str(round(qinj_max(i), digits_spec))]);
            end
            margin = min(margin, min(qinj_max(i) - qg_inj0(i), qg_inj0(i) - qinj_min(i)));
        end
    end
    p_dispatch = zeros(Ngen,1);
    for i = 1:Ngen
        p_dispatch(i) = gen(i,PG)+gen(i,end)*mpc.delta;
        if p_dispatch(i)>gen(i,GEN_STATUS)*gen(i,PMAX)+limit_tolerence || p_dispatch(i)<gen(i,GEN_STATUS)*gen(i,PMIN)-limit_tolerence
            violation_status.pg(end+1) = i;
            if suppress_warning == 0
                warning(['Pg violation at gen ', num2str(i), ': ', num2str(round(gen(i,PMIN), digits_spec)), ' ≤ ', num2str(round(p_dispatch(i), digits_spec)), ' ≤ ', num2str(round(gen(i,PMAX), digits_spec))]);
            end
            margin = min(margin, min(gen(i,PMAX) - p_dispatch(i), p_dispatch(i) - gen(i,PMIN)));
        end
    end
    angle_diff = zeros(Nbranch,1);
    for i = 1:Nbranch
        angle_diff(i) = bus(branch(i,F_BUS),VA) - bus(branch(i,T_BUS),VA);
        if angle_diff(i)>(branch(i,ANGMAX)+limit_tolerence) || angle_diff(i)<(branch(i,ANGMIN)-limit_tolerence)
            violation_status.angle(end+1) = i;
            if suppress_warning == 0
                warning(['Angle violation at line ', num2str(i), ': ', num2str(round(branch(i,ANGMIN), digits_spec)), ' ≤ ', num2str(round(angle_diff(i), digits_spec)), ' ≤ ', num2str(round(branch(i,ANGMAX), digits_spec))]);
            end
            margin = min(margin,min(deg2rad(branch(i,ANGMAX))-angle_diff(i),angle_diff(i)-deg2rad(branch(i,ANGMIN))));
        end
        if abs(s_line_fr0(i))>branch(i,RATE_A)+limit_tolerence || abs(s_line_to0(i))>branch(i,RATE_A)+limit_tolerence
            violation_status.line(end+1) = i;
            if suppress_warning == 0
                warning(['Line Flow violation at line ', num2str(i), ': ', num2str(round(abs(s_line_fr0(branch_index)), digits_spec)), ', ', num2str(round(abs(s_line_to0(branch_index)), digits_spec)), ' ≤ ', num2str(round(branch(i,RATE_A), digits_spec))]);
            end
            margin = min(margin, min(branch(i,rate_a) - abs(s_line_fr0(i)), branch(i,rate_a) - abs(s_line_to0(i))));
        end
    end
end

