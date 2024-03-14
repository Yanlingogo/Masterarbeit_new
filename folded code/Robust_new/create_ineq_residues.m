function [residues] = create_ineq_residues(g_u_pinj,g_l_pinj,g_u_qinj,g_l_qinj,g_u_vvcos,g_l_vvcos,g_u_vvsin,g_l_vvsin,g_u_vv,g_l_vv,K_plus,K_minus,Phi_u,Phi_l,v_u,v_l,idx_nslack)
    
    g_K1 = K_plus*[g_u_pinj;g_u_qinj;g_u_vvcos;g_u_vvsin;g_u_vv]+...
        K_minus*[g_l_pinj;g_l_qinj;g_l_vvcos;g_l_vvsin;g_l_vv]-...
        [Phi_u;v_u(idx_nslack);Pg_max;Qg_max;-Phi_l;-v_l(idx_nslack);-Pg_min;-Qg_min];
    
    residues = g_K1;
    %2Nbranch+2Nnslack
end

