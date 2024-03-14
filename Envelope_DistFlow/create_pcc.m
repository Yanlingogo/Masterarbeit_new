function [ineq_pcc] = create_pcc(p_pcc,q_pcc,P_u,P_l,Q_u,Q_l,l_u,l_l,R,X,Csline)
    pcc_pu = p_pcc - Csline*(-P_u + R*l_u);
    pcc_pl = -p_pcc + Csline*(-P_l+ R*l_l);
    pcc_qu = q_pcc - Csline*(-Q_u + X*l_u);
    pcc_ql = -q_pcc + Csline*(-Q_l + X*l_l);

    ineq_pcc = vertcat(pcc_pu,pcc_pl,pcc_qu,pcc_ql);
end

