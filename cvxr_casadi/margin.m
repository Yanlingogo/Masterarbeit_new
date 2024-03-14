function [max_margin] = margin(slack_pg,slack_qg,margin)

    margin1 = margin - slack_pg;
    margin2 = margin - slack_qg;

    max_margin = vertcat(margin1,margin2);
end

