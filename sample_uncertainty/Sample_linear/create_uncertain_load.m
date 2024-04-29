function [uncertain] = create_uncertain_load(Pg,alpha,A,e_l,e_u,Sigma,epsl,Pg_max,Pg_min,id_gen_nslack)
    eq_un = sum(alpha);
    un1 = alpha*A*e_u+sqrt((1-epsl)/epsl)*sqrt(alpha.^2*(A*Sigma*A'))+Pg(id_gen_nslack)-Pg_max(id_gen_nslack);
    un2 = -alpha*A*e_l+sqrt((1-epsl)/epsl)*sqrt(alpha.^2*(A*Sigma*A'))-Pg(id_gen_nslack)+Pg_min(id_gen_nslack);
    un3 = alpha*A*e_u-sqrt((1-epsl)/epsl)*sqrt(alpha.^2*(A*Sigma*A'))+Pg(id_gen_nslack)-Pg_max(id_gen_nslack);
    un4 = -alpha*A*e_l-sqrt((1-epsl)/epsl)*sqrt(alpha.^2*(A*Sigma*A'))-Pg(id_gen_nslack)+Pg_min(id_gen_nslack);
    uncertain = vertcat(eq_un,un1,un2,un3,un4);
end

