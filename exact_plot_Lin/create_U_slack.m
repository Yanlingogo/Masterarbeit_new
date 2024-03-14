function [U_slack] = create_U_slack(U,slack_U,Umax,Umin,id_nslack)
    ineq1 =  U(id_nslack) + slack_U - Umax(id_nslack);
    ineq2 = -U(id_nslack) + slack_U + Umin(id_nslack);

    U_slack = vertcat(ineq1,ineq2);

    %2Nbus-2
end

