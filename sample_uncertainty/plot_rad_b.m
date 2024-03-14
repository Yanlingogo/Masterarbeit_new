fig=figure; box on; grid on; hold on; set(fig, 'Position', [100, 100, 650, 550])

%% bounds on slack generator
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

%% bounds of combination
% fill the range
Points_exact = load('AC_conf.mat');
Points_exact = Points_exact.sortedPoints;
lightzs = [181, 180, 214]/255;
h2 = fill(Points_exact(:,1), Points_exact(:,2), lightzs);
set(h2, 'facealpha', 0.9, 'EdgeColor', 'none');
% draw the boundaries
load("rad_b1.mat");
load("rad_b2.mat");
load("rad_b3.mat");
load("rad_b4.mat");
load("rad_b5.mat");
load("rad_b8.mat");
boundaries = {b1,b2,b3,b4,b5,b8};

colorliste = [[237, 177, 32]/255;
              [0, 114, 189]/255;
              [0, 114, 189]/255;
              [0, 114, 189]/255;
              [65, 199, 101]/255;
              [0, 114, 189]/255;];
lineliste = {'-','-','-','-','-','--'};

for i = 1:size(boundaries,2)
    x = boundaries{i}(:,1);
    y = boundaries{i}(:,2);
    if i == 1
        plot(x,y,lineliste{i},'LineWidth',1.5,'Color',colorliste(i,:));
    else
        xi = linspace(min(x), max(x), 100);
        yi = pchip(x, y, xi);
        plot(xi, yi,lineliste{i},'LineWidth',1.5,'Color',colorliste(i,:));
    end
end


%% mark the intersections
p1 = [-1.0228,0.9471];
p2 = [0.0559,-0.9004];
p3 = [0.6901,-1.4102];
p4 = [1.2746,-1.757];
p5 = [1.6, -1.7099];
p6 = [1.6, -0.9856];
p7 = [-0.2626, 2];
p8 = [-0.6967, 2];
inter = [p1;p2;p3;p4;p5;p6;p7;p8];
deepviolet = [157, 105, 208]/255;
plot(inter(:,1), inter(:,2), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', deepviolet);

%% overlaid the results of LinDistFlow
% load("LinDistFlow.mat");
% lightorange = [250, 188, 113]/255;
% h3 = fill(valid_intersections(:,1), valid_intersections(:,2), lightorange, 'FaceAlpha',0.7, 'EdgeColor', 'none');
% set(h3, 'facealpha', 0.6, 'EdgeColor', 'none');

%% region considering uncertainty
load("rad_uncert.mat")
Points_un = sortedPoints;
lightred = [240, 128, 122]/255;
h3 = fill(Points_un(:,1), Points_un(:,2), lightred);
set(h3, 'facealpha', 0.6, 'EdgeColor', 'none');
%% set the range
xlim([-2 2]);
ylim([-2.2 2.2]);
set(gca, 'FontSize', 14,'FontName', 'Times New Roman');

x_label = xlabel('$p^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
y_label = ylabel('$q^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');