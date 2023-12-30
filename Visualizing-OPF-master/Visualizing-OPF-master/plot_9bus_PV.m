clc;

%% Get system data
mpc=test;
mpc.gen(:,end) = [0.5;0.3;0.2];
% mpc=runopf(mpc,mpoption('verbose',0,'out.all',0));
% mpc.gen(:,end+1)=[1;0;0];

%% get plot data
idx_slack = find(mpc.bus(:,2)==3);
mesh_intensity=20; 
plot_gen=find(mpc.gen(:,1)==idx_slack); 
mesh_axis=[0 2 -3 3];
[u1_plot,u2_plot,mesh_all,mesh_feasibility,u_base]=get_plot_data_PV(mpc,plot_gen,mesh_axis,mesh_intensity);

%% Plot figure
fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 450, 350]);
clim([0 400]); colormap([0.7461 0.832 0.9062 ;0.875*0.9 0.9492*0.9 0.6875*0.9; 0.875 0.9492 0.6875; 1 1 1]);

contourf(u1_plot,u2_plot,mesh_feasibility,[1e-4 1e-4])

pcolor_meshall=[0.25 0.25 0.25; ones(2,1)*[0 0.447 0.7410]; ones(4,1)*[0.929 0.694 0.125]; 0.494 0.184 0.556; ones(2,1)*[0.635 0.078 0.184]];
ptype_meshall={'-','-','--','-','--','-','--','-','--','-'};
for i=1:size(mesh_all,2)
    if size(mesh_all{i},1)~=1
        for j=1:size(mesh_all{i},3)
            if max(max(mesh_all{i}(:,:,j)))>0
                contour(u1_plot,u2_plot,mesh_all{i}(:,:,j),[1e-3 1e-3],ptype_meshall{i},'color',pcolor_meshall(i,:),'LineWidth',2);
            end
        end
    end
end

x_label=xlabel(['$p_{d,' num2str(plot_gen) '}$ (p.u.)']); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
y_label=ylabel(['$q_{d,' num2str(plot_gen) '}$ (p.u.)']); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
set(gca,'fontsize',15,'FontName','Times New Roman')
