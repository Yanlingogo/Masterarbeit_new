function [p] = create_coupling_branch_limit_p(v, ang, id_slack,connected_buses,id_cline, Gf, Bf, Gt, Bt)
    % active and reactive power along branch (k,l)
    U = v .* cos(ang);
    W = v .* sin(ang);
    % from slack bus
    pij = W(id_slack) .* (Gf(id_cline,:)*W + Bf(id_cline,:)*U) +...
        U(id_slack) .* (Gf(id_cline,:)*U - Bf(id_cline,:)*W);
    % to connected bus
%     pji = W(connected_buses) .* (Gt(id_cline,:)*W + Bt(id_cline,:)*U) +...
%         U(connected_buses) .* (Gt(id_cline,:)*U - Bt(id_cline,:)*W);
    % detect active line limit - 0 as no line limit
    p = pij;
    
end

