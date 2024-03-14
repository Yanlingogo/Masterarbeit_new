fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 650, 550])


%% exact plot
load("AC_mesh.mat")
Points_exact = sortedPoints;
lightbl = [181, 180, 214]/255;
h1 = fill(Points_exact(:,1), Points_exact(:,2), lightbl);
set(h1, 'facealpha', 0.9, 'EdgeColor', 'none');

%% sample points
% power range
deepgreen = [125, 170, 66]/255;
deepyellow = [250, 214, 61]/255;
x1 = [-2,2];
y_upper = [2 2];
y_lower = [-2 -2];
y1 = [-2.2,2.2];
x_lower = [-1.5 -1.5];
x_upper = [1.6 1.6];
plot(x1,y_upper, 'LineStyle','-','Color',deepgreen,'LineWidth',2);
plot(x1,y_lower, 'LineStyle','--','Color',deepgreen,'LineWidth',2);
plot(x_upper,y1, 'LineStyle','-','Color',deepyellow,'LineWidth',2);
plot(x_lower,y1, 'LineStyle','--','Color',deepyellow,'LineWidth',2);
% fill 
load('mesh045.mat');
Points_exact = Points;
lightgreen = [222, 235, 176]/255;
h2 = fill(Points_exact(:,1), Points_exact(:,2), lightgreen);
set(h2, 'facealpha', 0.7, 'EdgeColor', 'none');

% initial points
load('initial_mesh.mat')
Points_exact = Points;
h3 = plot(Points_exact(:,1),Points_exact(:,2),'rx', 'LineWidth',1.5);

% threshold =0.5
load('mesh045.mat')
A_str = join(string(Points), ",");
B_str = join(string(Points_exact), ",");
[C_str, ia] = setdiff(A_str, B_str, 'rows');
Points =  str2double(split(C_str, ","));
green = [123, 200, 91]/255;
h4 = scatter(Points(:,1), Points(:,2), 25, 'filled',...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', green,'MarkerEdgeColor', 'none');

%% explaination 


xlim([-2 2]);
ylim([-2.2 2.2]);
set(gca, 'FontSize', 14,'FontName', 'Times New Roman');
x_label = xlabel('$p^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');

y_label = ylabel('$q^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');

legend([h3, h4, h2, h1], {'Initial Points', 'Sampled Points','Filled Region','Exact AF'}, 'Interpreter', 'latex', 'FontName', 'Times New Roman','FontSize',12);
%% store the output
% exportgraphics(gca, 'legendOnly.pdf', 'ContentType', 'vector');
% 导出到pdf
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% filename = 'adap_rad05'; % 设定导出文件名
% print(gcf,filename,'-dpdf','-r0')
% close(gcf)
