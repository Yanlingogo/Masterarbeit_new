function [q] = create_coupling_branch_limit_q(v, ang, id_slack,connected_buses,id_cline, Gf, Bf, Gt, Bt)
    % active and reactive power along branch (k,l)
    U = v .* cos(ang);
    W = v .* sin(ang);
    % from slack bus to all other buses, > 0 means that DSO absorb active
    % power
    qij = W(id_slack) .* (Gf(id_cline,:)*U - Bf(id_cline,:)*W) -...
        U(id_slack) .* (Gf(id_cline,:)*W + Bf(id_cline,:)*U);
    % to connected bus
%     qji = W(connected_buses) .* (Gt(id_cline,:)*U - Bt(id_cline,:)*W) -...
%         U(connected_buses) .* (Gt(id_cline,:)*W + Bt(id_cline,:)*U);
    % detect active line limit - 0 as no line limit
    q = qij;
    
end

