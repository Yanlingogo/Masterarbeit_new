function [gvv] = create_ineq_gvv(g_u_vv,g_l_vv,v_u,v_l,v0,onoff_pq,Psi_vv0)

    gvv1 = -g_u_vv + v_u.^2 - 2*v0.*onoff_pq.*v_u;
    gvv2 = -g_u_vv + v_l.^2 - 2*v0.*onoff_pq.*v_l;

    gvv3 = g_l_vv - Psi_vv0 - 2*v0.*(v_u-v0) + 2*v0.*onoff_pq.*v_u;
    gvv4 = g_l_vv - Psi_vv0 - 2*v0.*(v_l-v0) + 2*v0.*onoff_pq.*v_l;
    
    gvv = vertcat(gvv1,gvv2,gvv3,gvv4);
    %4*Nbus
end

