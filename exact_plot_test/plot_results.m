clear;
clc;

mesh_all = load('mesh_all.mat');
mesh_all = mesh_all.mesh_all;
u1_plot  = load('u1_plot.mat');
u1_plot  = u1_plot.u1_plot;
u2_plot  = load('u2_plot.mat');
u2_plot  = u2_plot.u2_plot;


fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 850, 650]);
clim([0 400]); colormap([0.7461 0.832 0.9062 ;0.875*0.9 0.9492*0.9 0.6875*0.9; 0.875 0.9492 0.6875; 1 1 1]);
tol = 0;
pcolor_meshall=[0.25 0.25 0.25; ones(2,1)*[0 0.447 0.7410]; ones(4,1)*[0.929 0.694 0.125]; 0.494 0.184 0.556; ones(2,1)*[0.635 0.078 0.184]];
ptype_meshall={'-','-','--','-','--','-','--','-','--','-'};
for i=2:size(mesh_all,2)
    if size(mesh_all{i},1)~=1
        for j=1:size(mesh_all{i},3)
            if max(max(mesh_all{i}(:,:,j)))>=0 
                contour(u1_plot,u2_plot,mesh_all{i}(:,:,j),[-tol tol],ptype_meshall{i},'color',pcolor_meshall(i,:),'LineWidth',2);
            end
        end
    end
end
%% almost exact
% load('mesh_feasibility.mat')
% lightBlue = [0.4, 0.6, 1];
% contourf(u1_plot,u2_plot,mesh_feasibility,[-tol tol],'Color',lightBlue);
%% LinDistFlow
% Points = xlsread('filename.xlsx');
% lightRedColor = [1, 0.6, 0.6]; 
% plot(Points(:,1), Points(:,2), 'o-', 'Color', lightRedColor, 'MarkerSize', 6, 'MarkerFaceColor', lightRedColor);
% hold on; 
% h1 = fill(Points(:,1), Points(:,2), lightRedColor);
% set(h1, 'facealpha', 0.6, 'EdgeColor', 'none');
%% Sample AC model
load('sample_AC.mat')
sample_AC = sortedPoints;
lightorange = [1 0.5 0];
plot(sample_AC(:,1), sample_AC(:,2), 'o-', 'Color', lightorange, 'MarkerSize', 6, 'MarkerFaceColor', lightorange);
hold on; 
h3 = fill(sample_AC(:,1), sample_AC(:,2), lightorange);
set(h3, 'facealpha', 0.9, 'EdgeColor', 'none');
%% AC Model
Points_exact = load('exact_200.mat');
Points_exact = Points_exact.sortedPoints;
lightBlue = [0.4, 0.6, 1]; 
h2 = fill(Points_exact(:,1), Points_exact(:,2), lightBlue);
set(h2, 'facealpha', 0.6, 'EdgeColor', 'none');
%% add explanation
x_label=xlabel(['$P_{pcc,' num2str(1) '}$ (p.u.)']); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
y_label=ylabel(['$Q_{pcc,' num2str(1) '}$ (p.u.)']); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
title(['Feasible Region of Bus ' num2str(1) ' '],'FontSize',15,'FontName','Times New Roman');
set(gca,'fontsize',15,'FontName','Times New Roman')
% text(2.0, 0.5, 'Exact Feasible Region', 'Color', 'k', 'FontSize', 12,'FontName','Times New Roman');
% text(0.5, 0.5, 'hyperbox approximation', 'Color', 'k', 'FontSize', 12,'FontName','Times New Roman');
% add explanation for the bounds
unique_colors = [0 0.447 0.7410;
                 0.929 0.694 0.125;
                 0.494 0.184 0.556;
                 0.635 0.078 0.184;
                 lightorange;
                 lightBlue];
hold on;
lineWidths = [2*ones(4,1); 6*ones(2,1);];
h_fake = arrayfun(@(x) plot(NaN, NaN, 'Color', unique_colors(x,:), 'LineWidth', lineWidths(x)), 1:size(unique_colors, 1));
legend(h_fake, {'V bounds', 'P_{gen}/Q_{gen} bounds', 'S_{line} bounds', '\phi_{line} bounds','AC Mpdel(sampled)', 'AC Model'});