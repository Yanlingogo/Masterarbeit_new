function [base] = create_base(v_u,v_l,Delta_v_u,Delta_v_l,Phi_u,Phi_l,Delta_Phi_u,Delta_Phi_l,v0,Phi0,id)

    b1 = v_u(id) - v0(id);
    b2 = v_u(id) - v_l(id);
    b3 = v_u - Delta_v_u - v0;
    b4 = v_l - Delta_v_l - v0;
    b5 = Phi_u - Delta_Phi_u - Phi0;
    b6 = Phi_l - Delta_Phi_l - Phi0;

    % b_tu = theta_u(id);
    % b_tl = theta_l(id);
    % b_tPu = E'*theta_u - Phi_u;
    % b_tPl = E'*theta_l - Phi_l;
    % inequality
    b7 = -Phi_u + Phi_l;
    b8 = - v_u + v_l;
    base = vertcat(b1,b2,b3,b4,b5,b6,b7,b8);
    % 2 + 2*Nbus+2*Nbranch == 0
end

