function [margin_cons] = create_margin(slack_U,slack_P,slack_Q,margin)
    cons1 = margin - slack_U;
    cons2 = margin - slack_P;
    cons3 = margin - slack_Q;

    margin_cons = vertcat(cons1,cons2,cons3);
end

