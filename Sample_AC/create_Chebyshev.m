function [ineq] = create_Chebyshev(pg,pd_error,pg_f,K_e,mu_wind,mu_load,A_wind,A_load,cov_wind,cov_load)

    b_p = -pg + pg_f;
    ineq1 = -K_e*(b_p -mu_wind'*A_wind) + sqrt(A_wind'*cov_wind*A_wind);
    b_load =  -pd_error;
    ineq2 = -K_e*(b_load - mu_load'*A_load) + sqrt(A_load'*cov_load*A_load);

    ineq = vertcat(ineq1,ineq2);
end

