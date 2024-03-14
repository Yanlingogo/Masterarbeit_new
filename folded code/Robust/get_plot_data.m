function [u1_plot,u2_plot,mesh_all,mesh_feasibility,u_base]=plot_FeasibleRegion(mpc,plot_bus,mesh_axis,resolution)
limit_mode='vpqla';

%% Compute data for plotting
num_bus=size(mpc.bus,1);
num_gen=size(mpc.gen,1);
num_line=size(mpc.branch,1);
gen_status=mpc.gen(:,8);
pg_base=diag(gen_status)*mpc.gen(:,2)/mpc.baseMVA;
pl_base=mpc.bus(:,3)/mpc.baseMVA;
ql_base=mpc.bus(:,4)/mpc.baseMVA;

v_min=mpc.bus(:,13);
v_max=mpc.bus(:,12);
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

Cg = sparse(name2idx(mpc.gen(:,1)), (1:num_gen), ones(num_gen,1), num_bus, num_gen);
p_base=Cg*pg_base-pl_base;

u_base=p_base(plot_bus);

[u1_plot,u2_plot]=meshgrid(linspace(mesh_axis(1),mesh_axis(2),resolution),linspace(mesh_axis(3),mesh_axis(4),resolution));

mesh_all=cell(1,9); mesh_all(1,:)={0};
solve_mesh=ones(resolution);
V_max_mesh=-ones(resolution,resolution,num_bus);
V_min_mesh=-ones(resolution,resolution,num_bus);
Pg_max_mesh=-ones(resolution,resolution,num_gen);
Pg_min_mesh=-ones(resolution,resolution,num_gen);
Qg_max_mesh=-ones(resolution,resolution,num_gen);
Qg_min_mesh=-ones(resolution,resolution,num_gen);
Sline_mesh=-ones(resolution,resolution,num_line);
Etheta_max_mesh=-ones(resolution,resolution,num_line);
Etheta_min_mesh=-ones(resolution,resolution,num_line);

for i=1:resolution
    for j=1:resolution
        mpc_run=mpc;
        mpc_run.bus(plot_bus(1),3) = u1_plot(i, j)*mpc.baseMVA;
        mpc_run.bus(plot_bus(2),3) = u2_plot(i, j)*mpc.baseMVA;
        pf1 = mpc.bus(plot_bus(1),4) / mpc.bus(plot_bus(1),3);
        pf2 = mpc.bus(plot_bus(2),4) / mpc.bus(plot_bus(2),3);
        mpc_run.bus(plot_bus(1),4) = pf1 * u1_plot(i, j)*mpc.baseMVA;
        mpc_run.bus(plot_bus(2),4) = pf2 * u2_plot(i, j)*mpc.baseMVA;
        
        
%         mpc_run.bus(plot_bus,3)=(Cg(plot_bus(1),:)*pg_base-u1_plot(i,j))*mpc.baseMVA;
%         mpc_run.bus(plot_bus,3)=u_plot*mpc.baseMVA;
        
        mpc_run=runpf_cvxr(mpc_run);
        vmag_cur=mpc_run.bus(:,8);
        Pg_cur=mpc_run.gen(:,2)/mpc_run.baseMVA+(mpc_run.gen(:,end))*mpc_run.delta;
        Qg_cur=mpc_run.gen(:,3)/mpc_run.baseMVA;
        v_cplx=mpc_run.bus(:,8).*exp(1i*mpc_run.bus(:,9));
        
        Sline_fr=v_cplx(idx_fr).*conj(Yf*v_cplx); Sline_to=v_cplx(idx_to).*conj(Yt*v_cplx);
        Sline_cur=max(abs(Sline_fr),abs(Sline_to));
        Etheta_cur=(mpc_run.bus(idx_fr,9)-mpc_run.bus(idx_to,9));
        
        if mpc_run.success==0
            solve_mesh(i,j)=-1;
        else %大于0即表示满足约束
            tol = 1e-5;
            V_max_mesh(i,j,:)=v_max'+tol-vmag_cur';
            V_min_mesh(i,j,:)=vmag_cur'-v_min'+tol;
            Pg_max_mesh(i,j,:)=pg_max'+tol-Pg_cur';
            Pg_min_mesh(i,j,:)=Pg_cur'-pg_min'+tol;
            Qg_max_mesh(i,j,:)=qg_max'+tol-Qg_cur';
            Qg_min_mesh(i,j,:)=Qg_cur'-qg_min'+tol;
            Sline_mesh(i,j,:)=sline_max'+tol-Sline_cur';
            Etheta_max_mesh(i,j,:)=Etheta_max'+tol-Etheta_cur';
            Etheta_min_mesh(i,j,:)=Etheta_cur'-Etheta_min'+tol;
        end
    end
end

mesh_feasibility=solve_mesh; mesh_all{1}=-solve_mesh; % minus sign to avoid warning: mesh_all{1} encloses unsolvable region
if sum(limit_mode=='v'); mesh_feasibility=min(min(min(V_max_mesh,[],3),min(V_min_mesh,[],3)),mesh_feasibility); mesh_all{2}=V_max_mesh; mesh_all{3}=V_min_mesh; end
if sum(limit_mode=='p'); mesh_feasibility=min(min(min(Pg_max_mesh,[],3),min(Pg_min_mesh,[],3)),mesh_feasibility);mesh_all{4}=Pg_max_mesh; mesh_all{5}=Pg_min_mesh; end
if sum(limit_mode=='q'); mesh_feasibility=min(min(min(Qg_max_mesh,[],3),min(Qg_min_mesh,[],3)),mesh_feasibility); mesh_all{6}=Qg_max_mesh; mesh_all{7}=Qg_min_mesh; end
if sum(limit_mode=='l'); mesh_feasibility=min(min(Sline_mesh,[],3),mesh_feasibility); mesh_all{8}=Sline_mesh; end
if sum(limit_mode=='a'); mesh_feasibility=min(min(min(Etheta_max_mesh,[],3),min(Etheta_min_mesh,[],3)),mesh_feasibility); mesh_all{9}=Etheta_max_mesh; mesh_all{10}=Etheta_min_mesh; end


