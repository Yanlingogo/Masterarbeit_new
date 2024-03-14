function [ineq_current] = create_current(V_u,V_l,P_u,P_l,Q_u,Q_l,l_u,l_l,V_0,P_0,Q_0,l_0,J_p,J_m,Hes,from_bus)
    % theta_1 = [P_u - P_0; Q_u - Q_0; V_u(to_bus) - V_0(to_bus)];
    % theta_2 = [P_u - P_0; Q_u - Q_0; V_l(to_bus) - V_0(to_bus)];
    % theta_3 = [P_u - P_0; Q_l - Q_0; V_u(to_bus) - V_0(to_bus)];
    % theta_4 = [P_u - P_0; Q_l - Q_0; V_l(to_bus) - V_0(to_bus)];
    % theta_5 = [P_l - P_0; Q_u - Q_0; V_u(to_bus) - V_0(to_bus)];
    % theta_6 = [P_l - P_0; Q_u - Q_0; V_l(to_bus) - V_0(to_bus)];
    % theta_7 = [P_l - P_0; Q_l - Q_0; V_u(to_bus) - V_0(to_bus)];
    % theta_8 = [P_l - P_0; Q_l - Q_0; V_l(to_bus) - V_0(to_bus)];

    theta_1 = [P_u - P_0; Q_u - Q_0; V_u(from_bus) - V_0(from_bus)];
    theta_2 = [P_u - P_0; Q_u - Q_0; V_l(from_bus) - V_0(from_bus)];
    theta_3 = [P_u - P_0; Q_l - Q_0; V_u(from_bus) - V_0(from_bus)];
    theta_4 = [P_u - P_0; Q_l - Q_0; V_l(from_bus) - V_0(from_bus)];
    theta_5 = [P_l - P_0; Q_u - Q_0; V_u(from_bus) - V_0(from_bus)];
    theta_6 = [P_l - P_0; Q_u - Q_0; V_l(from_bus) - V_0(from_bus)];
    theta_7 = [P_l - P_0; Q_l - Q_0; V_u(from_bus) - V_0(from_bus)];
    theta_8 = [P_l - P_0; Q_l - Q_0; V_l(from_bus) - V_0(from_bus)];
    
    ineq_l_u_1 = -l_u + abs(2*(J_p'*theta_1 + J_m'*theta_8))+l_0;
    
    ineq_l_u_2 = -l_u + theta_1'*Hes*theta_1+l_0;
    ineq_l_u_3 = -l_u + theta_2'*Hes*theta_2+l_0;
    ineq_l_u_4 = -l_u + theta_3'*Hes*theta_3+l_0;
    ineq_l_u_5 = -l_u + theta_4'*Hes*theta_4+l_0;
    ineq_l_u_6 = -l_u + theta_5'*Hes*theta_5+l_0;
    ineq_l_u_7 = -l_u + theta_6'*Hes*theta_6+l_0;
    ineq_l_u_8 = -l_u + theta_7'*Hes*theta_7+l_0;
    ineq_l_u_9 = -l_u + theta_8'*Hes*theta_8+l_0;

    ineq_l_l_1 = l_l - (J_p'*theta_8 + J_m'*theta_1+l_0);

    ineq_current = vertcat(ineq_l_u_1,ineq_l_u_2,ineq_l_u_3,ineq_l_u_4,...
        ineq_l_u_5,ineq_l_u_6,ineq_l_u_7,ineq_l_u_8,ineq_l_u_9,ineq_l_l_1);
end

