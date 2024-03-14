clear;
clc;

mesh_all = load('mesh_all_rad.mat');
mesh_all = mesh_all.mesh_all;
u1_plot  = load('u1_plot_rad.mat');
u1_plot  = u1_plot.u1_plot;
u2_plot  = load('u2_plot_rad.mat');
u2_plot  = u2_plot.u2_plot;
load('mesh_feasibility_rad.mat')


fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 850, 650]);
hold all
clim([0 400]); colormap([0.7461 0.832 0.9062 ;0.875*0.9 0.9492*0.9 0.6875*0.9; 0.875 0.9492 0.6875; 1 1 1]);
tol = 0;
% contourf(u1_plot,u2_plot,mesh_feasibility,[-tol tol])
pcolor_meshall=[0.25 0.25 0.25; ones(2,1)*[0 0.447 0.7410]; ones(2,1)*[0.929 0.694 0.125]; ones(2,1)*[120, 171, 50]/255 ;0.494 0.184 0.556; ones(2,1)*[98, 64, 150]/255];
ptype_meshall={'-','-','--','-','--','-','--','-','--','-'};
for i=2:size(mesh_all,2)-2
    if size(mesh_all{i},1)~=1
        for j=1:size(mesh_all{i},3)
            if ismember(i,[2,3,4,5,6,7]) && j == 1
                continue
            elseif max(max(mesh_all{i}(:,:,j)))>=0
                contour(u1_plot,u2_plot,mesh_all{i}(:,:,j),[-tol tol],ptype_meshall{i},'color',pcolor_meshall(i,:),'LineWidth',1.5);
            end
        end
    end
end
deepgreen = [63, 141, 93]/255;
deepyellow = [231, 146, 53]/255;
yline(2, 'Color', deepgreen,'LineWidth',2.5)
yline(-2, 'Color', deepgreen,'LineWidth',2.5,'LineStyle','--')
xline(1.6, 'Color', deepyellow,'LineWidth',2.5)
xline(-1.5, 'Color', deepyellow,'LineWidth',2.5,'LineStyle','--')
%% almost exact
Points_exact = load('AC_conf.mat');
Points_exact = Points_exact.sortedPoints;
lightbl = [181, 180, 214]/255;
h2 = fill(Points_exact(:,1), Points_exact(:,2), lightbl);
set(h2, 'facealpha', 0.9, 'EdgeColor', 'none');
%% LinDistFlow
% Points = xlsread('filename.xlsx');
% lightRedColor = [1, 0.6, 0.6]; 
% plot(Points(:,1), Points(:,2), 'o-', 'Color', lightRedColor, 'MarkerSize', 6, 'MarkerFaceColor', lightRedColor);
% hold on; 
% h1 = fill(Points(:,1), Points(:,2), lightRedColor);
% set(h1, 'facealpha', 0.6, 'EdgeColor', 'none');
%% Sample AC model
% load('sample_AC.mat')
% sample_AC = sortedPoints;
% lightorange = [1 0.5 0];
% plot(sample_AC(:,1), sample_AC(:,2), 'o-', 'Color', lightorange, 'MarkerSize', 6, 'MarkerFaceColor', lightorange);
% hold on; 
% h3 = fill(sample_AC(:,1), sample_AC(:,2), lightorange);
% set(h3, 'facealpha', 0.9, 'EdgeColor', 'none');
%% AC Model
% Points_exact = load('exact_200.mat');
% Points_exact = Points_exact.sortedPoints;
% lightBlue = [0.4, 0.6, 1]; 
% h2 = fill(Points_exact(:,1), Points_exact(:,2), lightBlue);
% set(h2, 'facealpha', 0.6, 'EdgeColor', 'none');
%% add explanation
xlim([-2 2]);
ylim([-2.2 2.2]);
set(gca, 'FontSize', 14,'FontName', 'Times New Roman');

x_label = xlabel('$p^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');

y_label = ylabel('$q^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');

% text(2.0, 0.5, 'Exact Feasible Region', 'Color', 'k', 'FontSize', 12,'FontName','Times New Roman');
% text(0.5, 0.5, 'hyperbox approximation', 'Color', 'k', 'FontSize', 12,'FontName','Times New Roman');
% add explanation for the bounds
unique_colors = [0 0.447 0.7410;
                 0.929 0.694 0.125;
                 0.494 0.184 0.556;
                 0.635 0.078 0.184;
                 0.7461 0.832 0.9062];
% hold on;
% lineWidths = [2*ones(4,1); 5*ones(1,1);];
% h_fake = arrayfun(@(x) plot(NaN, NaN, 'Color', unique_colors(x,:), 'LineWidth', lineWidths(x)), 1:size(unique_colors, 1));
% legend(h_fake, {'V bounds', 'P_{gen}/Q_{gen} bounds', 'S_{line} bounds', '\phi_{line} bounds','AC Model(imprecise)'});