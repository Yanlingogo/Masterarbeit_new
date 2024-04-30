function [compensated_points] = compensation_Interpolation(base_point, endpoints, power_dispatch, resolution)
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
    Cf       = C + Ct;
    Cg_ns    = base_point.Cg_ns;
    Pd       = base_point.Pd;
    Qd       = base_point.Qd;
    u_0      = base_point.u_0;
    P_pf     = base_point.P_pf;
    Q_pf     = base_point.Q_pf;
    Mp       = base_point.Mp;
    Mq       = base_point.Mq;
    H_l      = base_point.H_l;

    Nbranch  = numel(branch_r);
    Nbus     = Nbranch+1;
    R = diag(branch_r);
    X = diag(branch_x);
    Z2       = R.^2 + X.^2;

    compensated_points = zeros(resolution,2);
    pcc_grid = vector_segment(endpoints, resolution);
    Jz_inv = inv(jac_z);
    %% setting:
    % participation factor or linear interpolation: 1 for participation, 2
    % for linear interpolation
    DER_setting = 1;
    % Whether to use correction
    %% compensation iteration
    for i = 1: resolution
        if DER_setting == 1
            p_step_iter = P_pf * (pcc_grid(1,i) - pcc_0(1));
            q_step_iter = Q_pf * (pcc_grid(2,i) - pcc_0(2));
        elseif DER_setting == 2
            p_step_iter = p_dispatch(:,1) + ((i-1)/(resolution-1))*(p_dispatch(:,2)-p_dispatch(:,1));
            q_step_iter = q_dispatch(:,1) + ((i-1)/(resolution-1))*(q_dispatch(:,2)-q_dispatch(:,1));
        end

        power_step = [p_step_iter;q_step_iter];
        d_z_pred = - jac_z\jac_u*(power_step-u_0);
        z_pred   = z_0 + d_z_pred;
        
        U_pred   = z_pred(1:Nbus);
        Pij_pred = z_pred(Nbus+1:Nbus+Nbranch);
        Qij_pred = z_pred(Nbus+Nbranch+1:Nbus+2*Nbranch);
        p_pcc    = z_pred(Nbus+2*Nbranch+1); 
        q_pcc    = z_pred(Nbus+2*Nbranch+2);
        L_pred   = z_pred(3*Nbus+1:end);
        e_st     = sparse(id_slack,1,1,Nbus,1);
        g_pred = [e_st'*U_pred-1;
                  C*U_pred - 2*(R*Pij_pred+X*Qij_pred) + Z2*L_pred;
                  C'*Pij_pred + Ct'*R*L_pred - e_st*p_pcc-Cg_ns*p_step_iter+Pd;
                  C'*Qij_pred + Ct'*X*L_pred - e_st*q_pcc-Cg_ns*q_step_iter+Qd;
                  L_pred .* Cf*U_pred - Pij_pred.^2  -  Qij_pred .^2];
                  % L_pred - (Pij_pred.^2+Qij_pred.^2)./U_pred(from_bus,:)];
        d_z_cor = -jac_z\g_pred;
        z_cor = z_0 + d_z_pred+ d_z_cor;
        U_cor = z_cor(1:Nbus);

        if any(U_cor < 0.81)
            % find the minimum voltage
            [U_cor_min, U_viol_index] = min(U_cor);

            % solve the equation for scaling factor
            syms scale
            d_z_pred_scale = -Jz_inv*jac_u*power_step*scale;
            z_pred_scale = z_0 + d_z_pred_scale;

            U_pred_scale   = z_pred_scale(1:Nbus);
            Pij_pred_scale = z_pred_scale(Nbus+1:Nbus+Nbranch);
            Qij_pred_scale = z_pred_scale(Nbus+Nbranch+1:Nbus+2*Nbranch);
            p_pcc_scale    = z_pred_scale(Nbus+2*Nbranch+1); 
            q_pcc_scale    = z_pred_scale(Nbus+2*Nbranch+2);
            L_pred_scale   = z_pred_scale(3*Nbus+1:end);
            g_pred_scale = [e_st'*U_pred_scale-1;
                      C*U_pred_scale - 2*(R*Pij_pred_scale+X*Qij_pred_scale) + Z2*L_pred_scale;
                      C'*Pij_pred_scale + Ct'*R*L_pred_scale - e_st*p_pcc_scale-Cg_ns*scale*p_step_iter+Pd;
                      C'*Qij_pred_scale + Ct'*X*L_pred_scale - e_st*q_pcc_scale-Cg_ns*scale*q_step_iter+Qd;
                      L_pred_scale .* Cf*U_pred_scale - Pij_pred_scale.^2  -  Qij_pred_scale.^2];
                      % L_pred_scale - (Pij_pred_scale.^2+Qij_pred_scale.^2)./U_pred_scale(from_bus,:)];
            d_z_cor_scale = -Jz_inv(U_viol_index,:)*g_pred_scale;
            z_cor_scale = z_0(U_viol_index) + d_z_pred_scale(U_viol_index)+ d_z_cor_scale;

            equation = z_cor_scale == 0.81;
            sol_scale = vpasolve(equation, scale, [0 1]);

            if isempty(sol_scale)
                error('There are no positive solutions');
            else
                % Filtering out positive solutions
                sol_scale = sol_scale(sol_scale > 0);    
                % Checking for positive solutions
                if isempty(sol_scale)
                    error('No positive solution was found');
                end
            end

            % validaton of scaling
            power_step = sol_scale*power_step;
            p_step_iter = power_step(1:size(power_step,1)/2);
            q_step_iter = power_step(size(power_step,1)/2+1:end);
            d_z_pred = - jac_z\jac_u*(power_step-u_0);
            z_pred   = z_0 + d_z_pred;

            U_pred   = z_pred(1:Nbus);
            Pij_pred = z_pred(Nbus+1:Nbus+Nbranch);
            Qij_pred = z_pred(Nbus+Nbranch+1:Nbus+2*Nbranch);
            p_pcc    = z_pred(Nbus+2*Nbranch+1); 
            q_pcc    = z_pred(Nbus+2*Nbranch+2);
            L_pred   = z_pred(3*Nbus+1:end);
            e_st     = sparse(id_slack,1,1,Nbus,1);
            g_pred = [e_st'*U_pred-1;
                      C*U_pred - 2*(R*Pij_pred+X*Qij_pred) + Z2*L_pred;
                      C'*Pij_pred + Ct'*R*L_pred - e_st*p_pcc-Cg_ns*p_step_iter+Pd;
                      C'*Qij_pred + Ct'*X*L_pred - e_st*q_pcc-Cg_ns*q_step_iter+Qd;
                      L_pred .* Cf*U_pred - Pij_pred.^2  -  Qij_pred.^2];
                      % L_pred - (Pij_pred.^2+Qij_pred.^2)./U_pred(from_bus,:)];
            d_z_cor = -jac_z\g_pred;
            z_cor = z_0 + d_z_pred+ d_z_cor;
            U_cor = z_cor(1:Nbus);

            tol = 1e-6;
            if any(U_cor < 0.81-tol)
                error('Scaling error, violation still exists')
            else
                disp('Scaling results, voltage bounds hold')
            end
        end

        Pij_cor = z_cor(Nbus+1:Nbus+Nbranch);
        Qij_cor = z_cor(Nbus+Nbranch+1:Nbus+2*Nbranch);
        L_cor     = (Pij_cor.^2 + Qij_cor.^2)./U_cor(from_bus,:);
        PQ_loss    = ([branch_r';branch_x']*L_cor);
        % pcc_comp   = z_cor(Nbus+2*Nbranch+1:Nbus+2*Nbranch+2);
    
        compensated_points(i,:) = pcc_grid(:,i) + PQ_loss;
    end
end

