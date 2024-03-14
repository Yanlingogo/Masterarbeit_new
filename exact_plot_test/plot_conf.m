clear;
clc;

mesh_all = load('mesh_all_rad.mat');
mesh_all = mesh_all.mesh_all;
u1_plot  = load('u1_plot_rad.mat');
u1_plot  = u1_plot.u1_plot;
u2_plot  = load('u2_plot_rad.mat');
u2_plot  = u2_plot.u2_plot;
load('mesh_feasibility_rad.mat');

fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 850, 650]);
clim([0 400]); colormap([0.7461 0.832 0.9062 ;0.875*0.9 0.9492*0.9 0.6875*0.9; 0.875 0.9492 0.6875; 1 1 1]);
tol = 0;
% contourf(u1_plot,u2_plot,mesh_feasibility,[-tol tol])
pcolor_meshall=[0.25 0.25 0.25; ones(2,1)*[0 0.447 0.7410]; ones(2,1)*[0.929 0.694 0.125]; ones(2,1)*[120, 171, 50]/255 ;0.494 0.184 0.556; ones(2,1)*[98, 64, 150]/255];
ptype_meshall={'-','-','--','-','--','-','--','-','--','-'};
for i=1:size(mesh_all,2)-2
    if size(mesh_all{i},1)~=1
        for j=1:size(mesh_all{i},3)
            if max(max(mesh_all{i}(:,:,j)))>=0 
                contour(u1_plot,u2_plot,mesh_all{i}(:,:,j),[-tol tol],ptype_meshall{i},'color',pcolor_meshall(i,:),'LineWidth',1.5);
            end
        end
    end
end

%% almost exact
Points_exact = load('AC_conf.mat');
Points_exact = Points_exact.sortedPoints;
lightbl = [193, 212, 230]/255;
h2 = fill(Points_exact(:,1), Points_exact(:,2), lightbl);
set(h2, 'facealpha', 0.9, 'EdgeColor', 'none');
%% almost exact
% Points_exact = load('rad.mat');
% Points_exact = Points_exact.sortedPoints;
% lightzs = [181, 180, 214]/255;
% h2 = fill(Points_exact(:,1), Points_exact(:,2), lightzs);
% set(h2, 'facealpha', 0.9, 'EdgeColor', 'none');
%% LinDistFlow
Points_exact = load('LinRad.mat');
Points_exact = Points_exact.sortedPoints;
lightorange = [250, 188, 113]/255;
h3 = fill(Points_exact(:,1), Points_exact(:,2), lightorange);
set(h3, 'facealpha', 0.6, 'EdgeColor', 'none');
%% observe point
% red = [240, 69, 49]/255;
% h3=plot(-0.85,1.54,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', red, 'MarkerFaceColor', red);
% p2 = [84, 84, 84]/255;
% h4=plot(-0.981411034156963,1.86212451429350,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', p2, 'MarkerFaceColor', p2);
%% range
 xlim([-1.5 1.6]);
 ylim([-2 2]);
%% legend
unique_colors = [0 0.447 0.7410;
                 0.929 0.694 0.125;
                 [120, 171, 50]/255;
                 [122, 50, 139]/255;
                 lightzs;
                 lightbl;
                 ];
hold on;
lineWidths = [2*ones(4,1); 5*ones(2,1);];
h_fake = arrayfun(@(x) plot(NaN, NaN, 'Color', unique_colors(x,:), 'LineWidth', lineWidths(x)), 1:size(unique_colors, 1));
h_sum = [h_fake, h3, h4];
legend(h_sum, {'V bounds', 'P_{gen} bounds','Q_{gen} bounds', 'S_{line} bounds', 'Radial Network','Mesh Network','A','B'},'NumColumns', 4,'FontSize',12);

%% legend
x_label = xlabel('$p^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');

y_label = ylabel('$q^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');

set(gca,'fontsize',20,'FontName','Times New Roman')
exportgraphics(gca, 'legendOnly.pdf', 'ContentType', 'vector');