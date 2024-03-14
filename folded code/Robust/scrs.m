function [mpc,result_cvxr] = scrs(mpc, max_iter_SCRS, varargin)
    % Set default value for phase_shift
    if nargin < 3  % Less than 3 arguments provided
        phase_shift = false;
    else
        phase_shift = varargin; % or simply use phase_shift = varargin{1}; if you expect more optional inputs
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

    Ngen = size(mpc.gen,1);

    gen = mpc.gen;
    bus = mpc.bus;

    result_cvxr = struct();
    result_cvxr.worst_cost  = zeros(1, max_iter_SCRS + 1);
    result_cvxr.solver_time = zeros(1, max_iter_SCRS + 1);
    result_cvxr.issue       = 0;
    result_cvxr.max_iter    = max_iter_SCRS;
    result_cvxr.obj         = mpc.cost;
    result_cvxr.pg          = (gen(:,GEN_STATUS).*gen(:,PG))/mpc.baseMVA;
    result_cvxr.vg          = gen(:,VG);
    result_cvxr.alpha       = gen(:,end);

    for iter = 2:result_cvxr.max_iter+1
        [mpc,sanity_check] = cvxrs2(mpc,"obj",[],phase_shift);
        %[mpc,sanity_check] = cvxrs_gurobi2(mpc,"obj",[],phase_shift);
        mpc = runpf_cvxr(mpc);
        test_runpf_cvxr(mpc);
        test_cvxrs(mpc,sanity_check);
        [violation_status,margin]=check_violation(mpc);

        result_cvxr.obj(:,iter) = mpc.cost;
        result_cvxr.worst_cost(iter) = sanity_check.obj;
        result_cvxr.pg(:,iter) = (mpc.gen(:,GEN_STATUS).*mpc.gen(:,PG))/mpc.baseMVA;
        result_cvxr.vg(:,iter) = mpc.gen(:,VG);
        result_cvxr.alpha(:,iter) = mpc.gen(:,end);
        %result_cvxr.solver_time(iter) = sanity_check.solve_time;

        dpg_step = norm(result_cvxr.pg(:,iter)-result_cvxr.pg(:,iter-1));
        dvg_step = norm(result_cvxr.vg(:,iter)-result_cvxr.vg(:,iter-1));
        if sanity_check.status
            status_str = 'true';
        else
            status_str = 'false';
        end
        fprintf('Iteration %d: %.2f %.3f %s %.4f %.4f\n', iter-1, result_cvxr.obj(iter),...
            sanity_check.obj, status_str, dpg_step, dvg_step);
        if dpg_step<1e-4 && dvg_step<1e-4
            result_cvxr.max_iter = iter; 
            break; 
        end
    end

end

