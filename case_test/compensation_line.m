function [compensated_points] = compensation_line(base_point, endpoints, power_dispatch, resolution)
    [x_power, ~] = size(power_dispatch);
    p_dispatch = power_dispatch(1:x_power/2,:);
    q_dispatch = power_dispatch(x_power/2+1:end,:);

    pcc_0    = base_point.pcc_0;
    z_0      = base_point.z_0;
    jac_z    = base_point.jac_z;
    jac_u    = base_point.jac_u;
    branch_r = base_point.branch_r;
    branch_x = base_point.branch_x;
    from_bus = base_point.from_bus;
    Nbranch  = numel(branch_r);
    Nbus     = Nbranch+1;


    compensated_points = zeros(resolution+1,2);
    pcc_grid = vector_segment(endpoints, resolution+1);
    for i = 0: resolution
        p_step_iter = p_dispatch(:,1) + (i/resolution)*(p_dispatch(:,2)-p_dispatch(:,1));
        q_step_iter = q_dispatch(:,1) + (i/resolution)*(q_dispatch(:,2)-q_dispatch(:,1));

        % delta_s = pcc_iter - pcc_0;
        % z_comp = z_0 + jac_z\jac_u*compensation_factor*delta_s;
        power_step = [p_step_iter;q_step_iter];
        z_comp     = z_0 - jac_z\jac_u*power_step;
        U_comp     = z_comp(1:Nbus,:);
        Pij_comp   = z_comp(Nbus+1:Nbus+Nbranch,:);
        Qij_comp   = z_comp(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
        L_comp     = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:);
        PQ_loss    = ([branch_r';branch_x']*L_comp);

        if all(U_comp >= 0.81)
            compensated_points(i+1,:) = pcc_grid(:,i+1) + PQ_loss;
        end
    end
end

