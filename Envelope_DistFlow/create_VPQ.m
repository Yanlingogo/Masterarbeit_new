function [eq_1] = create_VPQ(pg,qg,V_u,V_l,P_u,P_l,Q_u,Q_l,l_u,l_l,Cg,C,DR,DX_p,DX_m,Pd,Qd,V_0,Mp,Mq,H_p,H_m)
    p = Cg*pg - Pd;
    q = Cg*qg - Qd;

    Nbus = size(V_0,1);

    eq_V_u = V_u(2:end) - (V_0(1)*ones(Nbus-1,1) + Mp*p+Mq*q-H_p*l_l-H_m*l_u);
    eq_V_l = V_l(2:end) - (V_0(1)*ones(Nbus-1,1) + Mp*p+Mq*q-H_p*l_u-H_m*l_l);
    eq_P_u = P_u - (C*p-DR*l_l);
    eq_P_l = P_l - (C*p-DR*l_u);
    eq_Q_u = Q_u - (C*q -DX_p*l_l-DX_m*l_u);
    eq_Q_l = Q_l - (C*q -DX_p*l_u-DX_m*l_l);
    
    eq_1 = vertcat(eq_V_u,eq_V_l,eq_P_u,eq_P_l,eq_Q_u,eq_Q_l);
end

