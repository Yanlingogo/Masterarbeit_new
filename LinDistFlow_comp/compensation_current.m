function [compensated_points] = compensation_current(base_point, endpoints, power_dispatch, states, resolution)
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
    id_slack = base_point.id_slack;
    C        = base_point.C;
    Ct       = base_point.Ct;
    Cg_ns    = base_point.Cg_ns;
    Pd       = base_point.Pd;
    Qd       = base_point.Qd;
    u_0      = base_point.u_0;

    Nbranch  = numel(branch_r);
    Nbus     = Nbranch+1;
    Ngen     = size(Cg_ns,2)+1;
    R = diag(branch_r);
    X = diag(branch_x);
    Z2       = R.^2 + X.^2;





    % compensation without correction
    V_pred = states(1:Nbus,:);
    Pf_pred = states(Nbus+1:Nbus+Nbranch,:);
    Qf_pred = states(2*Nbus:Nbus+2*Nbranch,:);
    L_pred = (Pf_pred.^2+Qf_pred.^2)./V_pred(from_bus);
    pcc_grid = vector_segment(endpoints, resolution);
    L_segment = vector_segment(L_pred,resolution);
    PQ_loss = ([branch_r';branch_x']*L_segment);
    compensated_points = (pcc_grid + PQ_loss)';

    % compensation with correction
    % pcc_grid = vector_segment(endpoints, resolution);
    % states_segment = vector_segment(states, resolution);
    % compensated_points = zeros(resolution,2);
    % for i = 1: resolution
    %     U_pred   = states_segment(1:Nbus,i);
    %     if any(U_pred < 0.81)
    %         error('infeasible voltage')
    %     end
    %     Pij_pred = states_segment(Nbus+1:Nbus+Nbranch,i);
    %     Qij_pred = states_segment(2*Nbus:Nbus+2*Nbranch,i);
    %     p_pcc    = states_segment(Nbus+2*Nbranch+1,i);
    %     q_pcc    = states_segment(Nbus+2*Nbranch+Ngen+1,i);
    %     pg_pred  = states_segment(Nbus+2*Nbranch+2:Nbus+2*Nbranch+Ngen,i);
    %     qg_pred  = states_segment(Nbus+2*Nbranch+Ngen+2:Nbus+2*Nbranch+2*Ngen,i);
    %     L_pred   = (Pij_pred.^2 + Qij_pred.^2)./U_pred(from_bus);
    % 
    %     e_st     = sparse(id_slack,1,1,Nbus,1);
    %     g_pred   = [e_st'*U_pred-1;
    %                   C*U_pred - 2*(R*Pij_pred+X*Qij_pred) + Z2*L_pred;
    %                   C'*Pij_pred + Ct'*R*L_pred - e_st*p_pcc-Cg_ns*pg_pred+Pd;
    %                   C'*Qij_pred + Ct'*X*L_pred - e_st*q_pcc-Cg_ns*qg_pred+Qd;
    %                   L_pred - (Pij_pred.^2+Qij_pred.^2)./U_pred(from_bus,:)];
    %     d_z_cor  = -jac_z\g_pred;
    %     L_comp   = L_pred + d_z_cor(3*Nbus+1:3*Nbus+Nbranch);
    %     PQ_loss  = ([branch_r';branch_x']*L_comp);
    % 
    %     compensated_points(i,:) = (pcc_grid(:,i) + PQ_loss)';
    % end
end

