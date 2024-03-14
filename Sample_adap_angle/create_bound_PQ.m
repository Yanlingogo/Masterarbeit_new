function [bound] = create_bound_PQ(P, Q, alpha, id_gen_nslack,id_gen_slack)
    bound1 = -1*alpha*Q(id_gen_nslack) - P(id_gen_nslack);
    bound2 = alpha*Q(id_gen_nslack) - P(id_gen_nslack);
    %bound3 = (alpha*Q(id_gen_slack)).^2 - P(id_gen_slack).^2;
    bound = vertcat(bound1,bound2); 
end

