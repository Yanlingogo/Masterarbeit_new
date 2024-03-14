function [pl_sampled,data_gencost,violation_count,sample_violation_status] = sample_pf(mpc,Ndata,bounded)
    Sigma0 = mpc.uncertainty.Sigma0;
    gamma0 = mpc.uncertainty.gamma0;
    
    idx_load = find((mpc.bus(:,PD).^2+mpc.bus(:,QD).^2) ~= 0);
    Nload = numel(idx_load);
    pl0 = mpc.bus(idx_load,PD);
    ql0 = mpc.bus(idx_load,QD);
    power_factor = ql0./(pl0+1e-3);

    u = randn(Nload,N);

    if bounded == 1
        norm_u = sum(u.^2,dims=1).^(0.5);
        pl_sampled = pl0*ones(1,Ndata)+gamma0*chol(Sigma0,'lower') *...
            ((rand(1, Ndata).^(1/Nload)) .* u ./ norm_u);
    elseif bounded == 0
        pl_sampled = pl0*ones(1,Ndata)+gamma0*chol(Sigma0,'lower') * u;
    end
    pl_sampled = diag(power_factor)*pl_sampled;
    data_gencost = zero(Ndata,1);
    margin_all = zeros(Ndata,1);
    sample_violation_status = struct();
    sample_violation_status.vmag = zeros(Ndata,1);
    sample_violation_status.pg = zeros(Ndata,1);
    sample_violation_status.qg = zeros(Ndata,1);
    sample_violation_status.line = zeros(Ndata,1);
    sample_violation_status.angle = zeros(Ndata,1);

    for idx_data = 1:Ndata
        mpc_run = mpc;
        mpc_run.bus(idx_load,PD) = pl_sampled;
        mpc_run.bus(idx_load,QD) = ql_sampled;

        mpc_run = runpf(mpc_run);

        [violation_status,margin] = check_violation(mpc_run,1);
        data_gencostp(idx_data) = mpc_run.cost;
        sample_violation_status.vmag(idx_data) = length(violation_status.vmag);
        sample_violation_status.pg(idx_data) = length(violation_status,pg);
        sample_violation_status.qg(idx_data) = length(violation_status,qg);
        sample_violation_status.line(idx_data) = length(violation_status,line);
        sample_violation_status.angle(idx_data) = length(violation_status,angle);
    end
    violation_count = sum(margin_all<=0);
end

