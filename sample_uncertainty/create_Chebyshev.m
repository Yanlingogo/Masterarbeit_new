function [ineq] = create_Chebyshev(pg,eps,pg_f,e_sum,mu_wind,A_wind,sig_w)

    b_p = -pg + pg_f;
    ineq1 = -sqrt((1-eps)./eps).*(b_p - A_wind*mu_wind) + sqrt(sig_w);
    eq = sum(eps) - e_sum; % == 0

    ineq = vertcat(ineq1,eq);
end

