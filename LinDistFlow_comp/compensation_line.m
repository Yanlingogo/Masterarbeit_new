function [compensated_points] = compensation_line(base_point, endpoints, power_dispatch, resolution)
    [x_power, ~] = size(power_dispatch);
    p_dispatch = power_dispatch(1:x_power/2,:);
    q_dispatch = power_dispatch(x_power/2:end,:);
    factor_endpoint_p = p_dispatch./sum(p_dispatch);
    factor_endpoint_q = q_dispatch./sum(q_dispatch);

    pcc_0 = base_point.pcc_0;
    z_0 = base_point.z_0;
    jac_z = base_point.jac_z;
    jac_u = base_point.jac_u;
    branch_r = base_point.branch_r;
    branch_x = base_point.branch_x;

    compensated_points = zeros(resolution,2);
    for i = 0: resolution
        p_factor_iter = factor_endpoint_p(:,1) + (i/resolution)*(factor_endpoint_p(:,2)-factor_endpoint_p(:,1));
        q_factor_iter = factor_endpoint_q(:,1) + (i/resolution)*(factor_endpoint_q(:,2)-factor_endpoint_q(:,1));
        compensation_factor = blk_diag(p_factor_iter,q_factor_iter);

        pcc_iter = endpoints(:,1) + (i/resolution)*(endpoints(:,2) - endpoints(:,1));

        delta_s = pcc_iter - pcc_0;
        z_comp = z_0 + jac_z\jac_u*compensation_factor*delta_s;
        U_comp = z_comp(1:Nbus,:);
        Pij_comp = z_comp(Nbus+1:Nbus+Nbranch,:);
        Qij_comp = z_comp(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
        L_comp = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:);
        PQ_loss = ([branch_r';branch_x']*L_comp)';

        compensated_points(i+1,:) = pcc_grid + PQ_loss;
    end
end

