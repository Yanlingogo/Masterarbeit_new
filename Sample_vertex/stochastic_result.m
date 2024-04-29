mpc = ext2int(loadcase('case33_modified'));
% stochastic output (WT,PV)
mpc_stochastic = mpc;
hold all
%% generate stochastic result (normal distribution)
Numdata = 1000;
stochastic_AF = cell(Numdata,1);
for i = 1:Numdata
    % generate data for each loads(WT, PV), PD = 3
    mpc_stochastic.bus(:,3) = stochastic_loads(mpc);
    stochastic_AF{i} = fun_sample_vertex_ipopt(mpc_stochastic);
end
