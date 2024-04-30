function [compensated_points] = compensation_pf(base_point, vertexes, resolution)
    % load parameters at basepoint
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
    P_pf     = base_point.P_pf;
    Q_pf     = base_point.Q_pf;

    Nbranch  = numel(branch_r);
    Nbus     = Nbranch+1;
    R = diag(branch_r);
    X = diag(branch_x);
    Z2       = R.^2 + X.^2;
    Cf = C + Ct;

    % segmentation for vertexes
    compensated_points = zeros(resolution,2);
    t = linspace(0, 1, resolution);
    pcc_segmentation = interp1([0,1],vertexes,t');
    delta_s = pcc_segmentation - pcc_0;

    for i = 1:size(delta_s,1)

        p_step_iter = P_pf.*delta_s(i,1);
        q_step_iter = Q_pf.*delta_s(i,2);

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
        z_comp = z_0 + d_z_pred + d_z_cor;
    
        U_comp = z_comp(1:Nbus,:);
        Pij_comp = z_comp(Nbus+1:Nbus+Nbranch,:);
        Qij_comp = z_comp(Nbus+Nbranch+1:Nbus+2*Nbranch,:);
        L_comp = (Pij_comp.^2 + Qij_comp.^2)./U_comp(from_bus,:);
        PQ_loss = ([branch_r';branch_x']*L_comp)';

        PCC = pcc_segmentation(i,:) + PQ_loss;     
        compensated_points(i,:) = PCC;
    end
end

