function [inj] = create_ineq_inj(g_u_pinj,g_l_pinj,g_u_qinj,g_l_qinj,p_line_fr_u,q_line_fr_u,p_line_to_u,q_line_to_u,pg_opt,qg_opt,pl_opt,ql_opt,gamma_opt, slack_pg, slack_qg, slack_line, margin_opt,Cg,Cl,pg_max,pg_min,qg_max,qg_min,s_line_max)
    g_g1 = -g_u_pinj + Cg*pg_opt - Cl*pl_opt; 
    g_g2 =  g_l_pinj - Cg*pg_opt + Cl*pl_opt; 
    g_g3 = -g_u_qinj + Cg*qg_opt - Cl*ql_opt;
    g_g4 =  g_l_qinj - Cg*qg_opt + Cl*ql_opt;

    g_pgmax = -pg_max + pg_opt + slack_pg;   
    g_pgmin =  pg_min - pg_opt + slack_pg; 
    g_qgmax = -qg_max + qg_opt + slack_qg;
    g_qgmin =  qg_min - qg_opt + slack_qg;

    g_line1 = p_line_fr_u.^2+q_line_fr_u.^2+slack_line-s_line_max.^2;% <= s_line_max.^2
    g_line2 = p_line_to_u.^2+q_line_to_u.^2+slack_line-s_line_max.^2;% <= s_line_max.^2

    g_margin1 = -slack_pg + margin_opt;
    g_margin2 = -slack_qg + margin_opt;
    g_margin3 = -slack_line + margin_opt;
    g_margin4 = -gamma_opt + margin_opt;
 
    inj = vertcat(g_g1,g_g2,g_g3,g_g4,g_pgmax,g_pgmin,g_qgmax,g_qgmin,g_line1,g_line2,...
        g_margin1,g_margin2,g_margin3,g_margin4);
    %4*Nbus+4*Ngen+2*Nbranch+4
end

