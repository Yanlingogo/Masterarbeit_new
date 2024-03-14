function [CVaR] = create_CVaR(pg,b,p,t,z,s,pg0,ep,mu_u,sig_u)
    ineq1 = b + (1/ep)*(p+s);
    ineq2 = t^2 + sig_u - 4*z*s;
    ineq3 = -p + (pg-pg0) - b - t + mu_u + z;
    ineq4 = -z;
    ineq5 = -p;

    CVaR = vertcat(ineq1,ineq2,ineq3,ineq4,ineq5);
end

