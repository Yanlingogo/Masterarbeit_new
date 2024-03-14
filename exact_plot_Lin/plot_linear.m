load("mesh_all.mat");
load("mesh_feasibility.mat")
load("u1_plot.mat")
load("u2_plot.mat")
load("intersections.mat")
fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 850, 650]);
hold all
clim([0 400]); colormap([0.7461 0.832 0.9062 ;0.875*0.9 0.9492*0.9 0.6875*0.9; 0.875 0.9492 0.6875; 1 1 1]);
tol = 0;
% contourf(u1_plot,u2_plot,mesh_feasibility,[-tol tol])
pcolor_meshall=[0.25 0.25 0.25; ones(2,1)*[0 0.447 0.7410]; ones(2,1)*[0.929 0.694 0.125]; ones(2,1)*[120, 171, 50]/255 ;0.494 0.184 0.556; ones(2,1)*[98, 64, 150]/255];
ptype_meshall={'-','-','--','-','--','-','--','-','--','-'};
for i=2:size(mesh_all,2)
    if size(mesh_all{i},1)~=1
        for j=1:size(mesh_all{i},3)
            if ismember(i,[4,5,6,7]) && j == 1
                continue
            elseif max(max(mesh_all{i}(:,:,j)))>=0
                contour(u1_plot,u2_plot,mesh_all{i}(:,:,j),[-tol tol],ptype_meshall{i},'color',pcolor_meshall(i,:),'LineWidth',1.5);
            end
        end
    end
end

lightorange = [250, 188, 113]/255;
fill(valid_intersections(:,1), valid_intersections(:,2), lightorange, 'FaceAlpha',0.7, 'EdgeColor', 'none');

deepgreen = [114, 156, 97]/255;
deepyellow = [226, 201, 93]/255;
yline(2, 'Color', deepgreen,'LineWidth',2.5)
yline(-2, 'Color', deepgreen,'LineWidth',2.5,'LineStyle','--')
xline(1.6, 'Color', deepyellow,'LineWidth',2.5)
xline(-1.5, 'Color', deepyellow,'LineWidth',2.5,'LineStyle','--')

xlim([-2 2]);
ylim([-2.2 2.2]);
set(gca, 'FontSize', 14,'FontName', 'Times New Roman');
x_label = xlabel('$p^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');

y_label = ylabel('$q^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
% exportgraphics(gca, 'grid_lin.pdf', 'ContentType', 'vector');