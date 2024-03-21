function [u1_plot,u2_plot,mesh_all,mesh_feasibility,u_base]=get_plot_data(mpc,plot_gen,mesh_axis,resolution)
limit_mode='vpq';

%% Compute data for plotting
num_bus=size(mpc.bus,1);
num_gen=size(mpc.gen,1);
num_line=size(mpc.branch,1);
gen_status=mpc.gen(:,8);
pg_base=diag(gen_status)*mpc.gen(:,2)/mpc.baseMVA;
pl_base=mpc.bus(:,3)/mpc.baseMVA;
ql_base=mpc.bus(:,4)/mpc.baseMVA;

v_min=mpc.bus(:,13); U_min = v_min.^2;
v_max=mpc.bus(:,12); U_max = v_max.^2;
pg_max=diag(gen_status)*mpc.gen(:,9)/mpc.baseMVA;
pg_min=diag(gen_status)*mpc.gen(:,10)/mpc.baseMVA;
qg_max=diag(gen_status)*mpc.gen(:,4)/mpc.baseMVA;
qg_min=diag(gen_status)*mpc.gen(:,5)/mpc.baseMVA;
sline_max=mpc.branch(:,6)/mpc.baseMVA; sline_max(sline_max==0)=1e10;
Etheta_max=mpc.branch(:,13)*pi/180;
Etheta_min=mpc.branch(:,12)*pi/180;

name2idx=sparse(1,mpc.bus(:,1),1:num_bus);
idx_fr=name2idx(mpc.branch(:,1))';
idx_to=name2idx(mpc.branch(:,2))';
[Ybus, Yf, Yt]=makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
idx_sgen = find(mpc.gen(:,1)==plot_gen);%plot 对应节点的编号

Cg = sparse(name2idx(mpc.gen(:,1)), (1:num_gen), ones(num_gen,1), num_bus, num_gen);
p_base=Cg*pg_base-pl_base;

u_base=p_base(plot_gen);

[u1_plot,u2_plot]=meshgrid(linspace(mesh_axis(1),mesh_axis(2),resolution),linspace(mesh_axis(3),mesh_axis(4),resolution));

mesh_all=cell(1,7); mesh_all(1,:)={0};
solve_mesh=ones(resolution);
U_max_mesh=-ones(resolution,resolution,num_bus);
U_min_mesh=-ones(resolution,resolution,num_bus);
Pg_max_mesh=-ones(resolution,resolution,num_gen);
Pg_min_mesh=-ones(resolution,resolution,num_gen);
Qg_max_mesh=-ones(resolution,resolution,num_gen);
Qg_min_mesh=-ones(resolution,resolution,num_gen);

for i=1:resolution
    for j=1:resolution
%         mpc_run=mpc;
%         mpc_run.bus(plot_bus(1),3) = u1_plot(i, j)*mpc.baseMVA;
%         mpc_run.bus(plot_bus(2),3) = u2_plot(i, j)*mpc.baseMVA;
%         pf1 = mpc.bus(plot_bus(1),4) / mpc.bus(plot_bus(1),3);
%         pf2 = mpc.bus(plot_bus(2),4) / mpc.bus(plot_bus(2),3);
%         mpc_run.bus(plot_bus(1),4) = pf1 * u1_plot(i, j)*mpc.baseMVA;
%         mpc_run.bus(plot_bus(2),4) = pf2 * u2_plot(i, j)*mpc.baseMVA;

        mpc_run = mpc;
        mpc_run.gen(idx_gen,2) = u1_plot(i,j)*mpc.baseMVA;
        mpc_run.gen(idx_gen,3) = u2_plot(i,j)*mpc.baseMVA;

%         mpc_run.bus(plot_bus,3)=(Cg(plot_bus(1),:)*pg_base-u1_plot(i,j))*mpc.baseMVA;
%         mpc_run.bus(plot_bus,3)=u_plot*mpc.baseMVA;
        
        [results, margin]=margin_search(mpc_run);

        U_cur = results.U;
        Pg_cur = results.pg;
        Qg_cur = results.qg;
                
        solve_mesh(i,j)=margin;

        tol = 0;
        U_max_mesh(i,j,:)=U_max'+tol-U_cur';
        U_min_mesh(i,j,:)=U_cur'-U_min'+tol;
        Pg_max_mesh(i,j,:)=pg_max'+tol-Pg_cur';
        Pg_min_mesh(i,j,:)=Pg_cur'-pg_min'+tol;
        Qg_max_mesh(i,j,:)=qg_max'+tol-Qg_cur';
        Qg_min_mesh(i,j,:)=Qg_cur'-qg_min'+tol;
        
    end
end

mesh_feasibility=solve_mesh; mesh_all{1}=-solve_mesh; % minus sign to avoid warning: mesh_all{1} encloses unsolvable region
if sum(limit_mode=='v'); mesh_feasibility=min(min(min(U_max_mesh,[],3),min(U_min_mesh,[],3)),mesh_feasibility); mesh_all{2}=U_max_mesh; mesh_all{3}=U_min_mesh; end
if sum(limit_mode=='p'); mesh_feasibility=min(min(min(Pg_max_mesh,[],3),min(Pg_min_mesh,[],3)),mesh_feasibility);mesh_all{4}=Pg_max_mesh; mesh_all{5}=Pg_min_mesh; end
if sum(limit_mode=='q'); mesh_feasibility=min(min(min(Qg_max_mesh,[],3),min(Qg_min_mesh,[],3)),mesh_feasibility); mesh_all{6}=Qg_max_mesh; mesh_all{7}=Qg_min_mesh; end

end


