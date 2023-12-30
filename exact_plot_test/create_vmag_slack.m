function [v_slack] = create_vmag_slack(v,slack_v,vmax,vmin,id_nslack)
    ineq1 =  v(id_nslack) + slack_v - vmax(id_nslack);
    ineq2 = -v(id_nslack) + slack_v + vmin(id_nslack);

    v_slack = vertcat(ineq1,ineq2);

    %2Nbus-2
end

