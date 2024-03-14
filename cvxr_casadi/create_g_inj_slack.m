function [g_inj] = create_g_inj_slack(g_u_pinj,g_l_pinj,g_u_qinj,g_l_qinj,pg_opt,qg_opt,slack_pg,slack_qg,Cg,Cl,pl0,ql0,pg_max,pg_min,qg_max,qg_min)
   g1 = - g_u_pinj + Cg*pg_opt - Cl*pl0;
   g2 = g_l_pinj - Cg*pg_opt + Cl*pl0;
   g3 = -g_u_qinj + Cg*qg_opt - Cl*ql0;
   g4 = g_l_qinj - Cg*qg_opt + Cl*ql0;
   g5 = pg_opt + slack_pg - pg_max;
   g6 = -pg_opt + slack_pg + pg_min;
   g7 = qg_opt + slack_qg - qg_max;
   g8 = -qg_opt +slack_qg + qg_min;

   g_inj = vertcat(g1,g2,g3,g4,g5,g6,g7,g8);
end

