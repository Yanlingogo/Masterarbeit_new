function [compensated_points] = compensation_collection(base_point, endpoints, power_dispatch, resolution)
    [x_power, ~] = size(power_dispatch);
    p_dispatch = power_dispatch(1:x_power/2,:);
    q_dispatch = power_dispatch(x_power/2+1:end,:);

    % load the data from base point
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

    V_0      = z_0(1:Nbus);
    R = diag(branch_r);
    X = diag(branch_x);
    Z2       = R.^2 + X.^2;
    % compensation
    compensated_points = zeros(resolution+1,2);
    pcc_grid = vector_segment(endpoints, resolution+1);
    for i = 0: resolution
        p_step_iter = p_dispatch(:,1) + (i/resolution)*(p_dispatch(:,2)-p_dispatch(:,1));
        q_step_iter = q_dispatch(:,1) + (i/resolution)*(q_dispatch(:,2)-q_dispatch(:,1));

        power_step = [p_step_iter;q_step_iter];
        d_z_pred = - jac_z\jac_u*(power_step-u_0);
        z_pred   = z_0 + d_z_pred;
        V_pred   = z_pred(1:Nbus);
        Pij_pred = z_pred(Nbus+1:Nbus+Nbranch);
        Qij_pred = z_pred(Nbus+Nbranch+1:Nbus+2*Nbranch);
        p_pcc    = z_pred(Nbus+2*Nbranch+1); 
        q_pcc    = z_pred(Nbus+2*Nbranch+2);
        L_pred   = z_pred(3*Nbus+1:end);
        e_st     = sparse(id_slack,1,1,Nbus,1);
        g_pred = [e_st'*V_pred-1;
                  C*V_pred - 2*(R*Pij_pred+X*Qij_pred) + Z2*L_pred;
                  C'*Pij_pred + Ct'*R*L_pred - e_st*p_pcc-Cg_ns*p_step_iter+Pd;
                  C'*Qij_pred + Ct'*X*L_pred - e_st*q_pcc-Cg_ns*q_step_iter+Qd;
                  L_pred - (Pij_pred.^2+Qij_pred.^2)./V_pred(from_bus,:)];
        d_z_cor = -jac_z\g_pred;
        z_comp = z_0 + d_z_pred + d_z_cor;
        V_comp = z_comp(1:Nbus,:);

        % fix the voltage, where violation happens
        if any(V_comp < 0.81)
            % index_V = find(V_comp<0.81);
            [ ~,index_V] = min(V_comp);
            remaining_index_U = setdiff(1:size(jac_z,2),index_V);
            jac_z_remaining = jac_z(:,remaining_index_U);
            jac_z_violation = jac_z(:,index_V);

            % first validation: equalities are still satisfied?
            d_z_pred_remaining = d_z_pred(remaining_index_U);
            d_z_pred_violation = d_z_pred(index_V);
            validation_1 = jac_z_remaining*d_z_pred_remaining + jac_u*(power_step-u_0)...
                +jac_z_violation*d_z_pred_violation;

            tol = 1e-8;
            if all(abs(validation_1) <= tol)
                disp('Decomposition of Jacobian results are correct')
            else
                error('Decomposition result error')
            end

            % Second estimation
            index_Q = size(power_step,1)-size(index_V)+1 : size(power_step,1);
            remaining_index_Q = setdiff(1:size(jac_u,2),index_Q);
            jac_u_remaining = jac_u(:,remaining_index_Q);
            jac_u_Q         = jac_u(:,index_Q);

            jac_combi = [jac_z_remaining jac_u_Q];
            delta_V = 0.81 - V_0(index_V);
            delta_u = power_step(remaining_index_Q)-u_0(remaining_index_Q); 
            % new_estimation = -jac_combi\(jac_u_remaining*delta_u+jac_z_violation*delta_V);
            new_estimation = -pinv(full(jac_combi))*(jac_u_remaining*delta_u+jac_z_violation*delta_V);


            d_z_pred_new = [new_estimation(1:index_V-1);delta_V;new_estimation(index_V:end-size(index_V))];

            %% correction again?
            
            % no correction
            z_comp   = z_0 + d_z_pred_new;
            V_comp   = z_comp(1:Nbus);
            % second validation: voltage fixed?
            if all(V_comp(index_V) <= 0.81+tol) && all(V_comp(index_V) >= 0.81-tol)
                disp('The voltage minimum is met')
            else
                error('Voltage minimum failed to be fixed')
            end
            Pij_comp = z_comp(Nbus+1:Nbus+Nbranch);
            Qij_comp = z_comp(Nbus+Nbranch+1:Nbus+2*Nbranch);
            L_comp     = (Pij_comp.^2 + Qij_comp.^2)./V_comp(from_bus,:);
            PQ_loss    = ([branch_r';branch_x']*L_comp);
            compensated_points(i+1,:) = pcc_grid(:,i+1) + PQ_loss;
            
        end
        V_comp   = z_comp(1:Nbus);
        Pij_comp = z_comp(Nbus+1:Nbus+Nbranch);
        Qij_comp = z_comp(Nbus+Nbranch+1:Nbus+2*Nbranch);
        L_comp     = (Pij_comp.^2 + Qij_comp.^2)./V_comp(from_bus,:);
        PQ_loss    = ([branch_r';branch_x']*L_comp);
    
        compensated_points(i+1,:) = pcc_grid(:,i+1) + PQ_loss;

end