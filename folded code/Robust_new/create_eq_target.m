function [target] = create_eq_target(pg_opt,qg_opt,pl_opt,ql_opt,v_u,mpc,idx_load,id_gen)

    target1 = pg_opt - mpc.gen(:,2)/mpc.baseMVA; % 2 PG
    target2 = qg_opt - mpc.gen(:,3)/mpc.baseMVA; % 3 QG
    target3 = pl_opt - mpc.bus(idx_load,3)/mpc.baseMVA; % 3 PD
    target4 = ql_opt - mpc.bus(idx_load,4)/mpc.baseMVA; % 4 QD
    target5 = v_u(id_gen) - mpc.bus(id_gen,8); % 8 VM

    target = vertcat(target1,target2,target3,target4,target5);
    % 3Ngen+2Nload
end

