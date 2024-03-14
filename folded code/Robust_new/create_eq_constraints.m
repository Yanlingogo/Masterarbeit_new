function [eq] = create_eq_constraints(v0,Phi0,id_gen,v_u,v_l,D_v_u,D_v_l,Phi_u,Phi_l,D_Phi_u,D_Phi_l)
    eq1 = v_u(id_gen) - v_l(id_gen);
    eq2 = v_u - v0 - D_v_u;
    eq3 = v_l - v0 - D_v_l;
    eq4 = Phi_u - Phi0 - D_Phi_u;
    eq5 = Phi_l - Phi0 - D_Phi_l;
    eq = vertcat(eq1,eq2,eq3,eq4,eq5);
    %Ngen+2Nbus+2Nbranch
end

