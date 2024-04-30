function [DER_output] = create_DER_participation_factor(Pg,Qg,Pgmax,Qgmax,gen_nslack,base_point)
    % P_p = -Pgmax(gen_nslack)/sum(Pgmax(gen_nslack));
    % Q_p = -Qgmax(gen_nslack)/sum(Qgmax(gen_nslack));
    % 
    % s_0 = base_point.pcc_0;
    % gen_slack = 1;
    % 
    % DER_ouput_1 = Pg(gen_nslack) - P_p.*(Pg(gen_slack) - s_0(1));
    % DER_ouput_2 = Qg(gen_nslack) - Q_p.*(Qg(gen_slack) - s_0(2));

    % DER_output = vertcat(DER_ouput_1, DER_ouput_2);

    DER_output_p = diff(Pg(gen_nslack));
    DER_output_q = diff(Qg(gen_nslack));

    DER_output = vertcat(DER_output_p, DER_output_q);
end


