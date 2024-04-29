function [gen_slack] = create_gen_slack(Pg,Qg,slack_P,slack_Q,Pmax,Pmin,Qmax,Qmin)
    P1 = Pg + slack_P - Pmax;
    P2 = -Pg + slack_P + Pmin;
    Q1 = Qg + slack_Q - Qmax;
    Q2 = -Qg + slack_Q + Qmin;

    gen_slack = vertcat(P1,P2,Q1,Q2);
    %4Ngen
end

