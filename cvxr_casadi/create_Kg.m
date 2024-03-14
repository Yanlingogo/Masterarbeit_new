function [K_eq] = create_Kg(g_u_pinj,g_l_pinj,g_u_qinj,g_l_qinj,g_u_vvcos,g_l_vvcos,g_u_vvsin,g_l_vvsin,g_u_vv,g_l_vv,v_u,v_l,Phi_u,Phi_l,K_plus,K_minus)
 K_eq = K_plus*[g_u_pinj; g_u_qinj; g_u_vvcos; g_u_vvsin; g_u_vv]+K_minus*[g_l_pinj; g_l_qinj; g_l_vvcos; g_l_vvsin; g_l_vv]-[Phi_u;v_u;-Phi_l;-v_l];
end

