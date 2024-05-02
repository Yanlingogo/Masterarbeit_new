function [line_points] = compensation_line(vertexes,dispatch,base_point, voltage_check,comp_type,resolution,method)
    line_points = cell(size(vertexes,1),1);
    numlines = size(vertexes,1);
    for i = 1: numlines
        % endpoints of lines
        endpoints = [vertexes(i,:)' vertexes(mod(i,numlines)+1,:)'];
        % DER dispatch at endpoints
        power_dispatch = [dispatch{i} dispatch{mod(i,numlines)+1}];

        line_points{i} = compensation(base_point, endpoints, power_dispatch,...
            voltage_check,comp_type,resolution,method);
    end
end
% compensation for each line
function [compensated_points] = compensation(base_point, endpoints, power_dispatch,voltage_check,comp_type,resolution,method)
    [x_power, ~] = size(power_dispatch);
    p_dispatch = power_dispatch(1:x_power/2,:);
    q_dispatch = power_dispatch(x_power/2+1:end,:);

    pcc_0    = base_point.pcc_0;
    z_0      = base_point.z_0;
    jac_z    = base_point.jac_z; Jz_inv = inv(jac_z);
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
    L_0      = base_point.L_0;

    Nbranch  = numel(branch_r);
    Nbus     = Nbranch+1;
    R = diag(branch_r);
    X = diag(branch_x);
    Z2       = R.^2 + X.^2;

    compensated_points = zeros(resolution,2);
    pcc_grid = vector_segment(endpoints, resolution);
    %% setting:
    % participation factor or linear interpolation: 1 for participation, 2
    % for linear interpolation
    DER_setting = 2;
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
        Vmin_nodal = 0.81;
        if any(U_pred < Vmin_nodal) && voltage_check == 1
            % find the minimum voltage
            [U_pred_min, U_viol_index] = min(U_pred);
            delta_U = U_pred_min - Vmin_nodal;
            U_viol_index

            % 1. equal weight scaling
            % intermediate = -jac_z\jac_u*power_step;
            % scale_U = 1 - (delta_U/intermediate(U_viol_index));
            
            % 2. Fixed scaling
            % scale_U = 0.5;
            
            % 3. variable weight scaling
            % weight = abs(power_step)/norm(power_step);
            % intermediate = -Jz_inv*jac_u*diag(power_step);
            % intermediate = intermediate(U_viol_index,:);
            % scale_U = -(delta_U-sum(intermediate))/(intermediate*weight);
            % scale_U = diag(weight*scale_U);
            % 
            % if any(scale_U > 1)
            %     warning('scaling error')
            % end

            % 4. estimate the relationship between voltage and pwer
            intermediate = 2*(Mp*Cg_ns(2:end,:)*p_step_iter+Mq*Cg_ns(2:end,:)*q_step_iter);
            scale_U = 1 - delta_U/intermediate(U_viol_index-1);

            % scale_U = 1;

            power_step = scale_U*power_step;
            p_step_iter = power_step(1:size(power_step,1)/2);
            q_step_iter = power_step(size(power_step,1)/2+1:end);
            d_z_pred = - jac_z\jac_u*(power_step-u_0);
            z_pred   = z_0 + d_z_pred;
            U_pred   = z_pred(1:Nbus);

        end
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
        d_z_cor = -jac_z\g_pred;
        % z_comp = z_0 + d_z_pred;

        z_comp = z_0 + d_z_pred+ d_z_cor;
        
        U_comp   = z_comp(1:Nbus);
        % if any(U_comp < 0.81)
        %     error('correction error')
        % end
        Pij_comp = z_comp(Nbus+1:Nbus+Nbranch);
        Qij_comp = z_comp(Nbus+Nbranch+1:Nbus+2*Nbranch);
        L_comp     = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:);
        if method == 2 % current recover
            L_comp = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:) - L_0;
        end
        PQ_loss    = ([branch_r';branch_x']*L_comp);
        pcc_comp   = z_comp(Nbus+2*Nbranch+1:Nbus+2*Nbranch+2);
        if comp_type == 1
            compensated_points(i,:) = pcc_grid(:,i) + PQ_loss;
        elseif comp_type == 2
            compensated_points(i,:) = pcc_comp;
        end

    end
end

function [points] = vector_segment(A,resolution)
    [row, ~] = size(A);
    points = zeros(row, resolution);
    for i = 1:row
        points(i,:) = linspace(A(i,1),A(i,2),resolution);
    end
end

