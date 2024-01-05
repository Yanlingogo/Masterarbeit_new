function [h] = create_local_branch_limit_rec(v, ang, Gf, Bf, Gt, Bt, Fmax, from_bus,to_bus)
    % active and reactive power along branch (k,l)
    U = v .* cos(ang);
    W = v .* sin(ang);
    % from side
    pij = W(from_bus) .* (Gf*W + Bf*U) + U(from_bus) .* (Gf*U - Bf*W);
    qij = W(from_bus) .* (Gf*U - Bf*W) - U(from_bus) .* (Gf*W + Bf*U);
    % to side
    pji = W(to_bus) .* (Gt*W + Bt*U) + U(to_bus) .* (Gt*U - Bt*W);
    qji = W(to_bus) .* (Gt*U - Bt*W) - U(to_bus) .* (Gt*W + Bt*U);
    % detect active line limit - 0 as no line limit
    idx_limit = find(Fmax);
    Sij = pij(idx_limit).^2 + qij(idx_limit).^2 - Fmax(idx_limit).^2;
    Sji = pji(idx_limit).^2 + qji(idx_limit).^2 - Fmax(idx_limit).^2;
    % inequality constraints for line limit
    h   = vertcat(Sij, Sji);

    
end