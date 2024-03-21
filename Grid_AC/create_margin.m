function [margin_cons] = create_margin(slack_v,slack_phi,slack_P,slack_Q,slack_line,margin)
    cons1 = margin - slack_v;
    cons2 = margin - slack_phi;
    cons3 = margin - slack_P;
    cons4 = margin - slack_Q;
    cons5 = margin - slack_line;

    margin_cons = vertcat(cons1,cons2,cons3,cons4,cons5);
end

