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
Points_exact = load('AC_mesh.mat');
Points_exact = Points_exact.sortedPoints;
lightbl = [193, 212, 230]/255;
h2 = fill(Points_exact(:,1), Points_exact(:,2), lightbl);
set(h2, 'facealpha', 0.9, 'EdgeColor', 'none');
% draw the boundaries
load("mesh_b1.mat"); 
load("mesh_b2.mat"); 
load("mesh_b3.mat"); 
load("mesh_b4.mat"); 
load("mesh_b5.mat");
load("mesh_b6u.mat");
load("mesh_b6l.mat");
load("mesh_b8.mat"); 
load("mesh_b9.mat");
boundaries = {b1,b2,b3,b4,b5,b6_u,b6_l,b8,b9};

colorliste = [[237, 177, 32]/255; 
              [0, 114, 189]/255;  
              [0, 114, 189]/255;  
              [0, 114, 189]/255; 
              [246, 152, 50]/255;  
              [246, 152, 50]/255;  
              [246, 152, 50]/255;
              [0, 114, 189]/255; 
              [0, 114, 189]/255;];
lineliste = {'-','-','-','-','-','-','-','--','--'};

for i = 1:size(boundaries,2)
    x = boundaries{i}(:,1);
    y = boundaries{i}(:,2);
    % if i == 1
        plot(x,y,lineliste{i},'LineWidth',1.5,'Color',colorliste(i,:));
    % else
        % xi = linspace(min(x), max(x), 100);
        % yi = pchip(x, y, xi);
        % plot(xi, yi,lineliste{i},'LineWidth',1.5,'Color',colorliste(i,:));
    % end
end


%% mark the intersections
p1 = [-0.9555,1.9817];
p2 = [-0.9165,2];
p3 = [-1.1958,0.0787];
p4 = [-0.3206,-1.1132];
p5 = [0.7057, -1.9712];
p6 = [0.8687, -1.9535];
p7 = [1.0994, -1.9055];
p8 = [1.6, -1.5098];
p9 = [1.6, -0.0927];
p10 = [1.1865,0.3975];
p11 = [0.047,2];
inter = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11];
deepbl2 = [79, 108, 220]/255;
plot(inter(:,1), inter(:,2), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', deepbl2);
%% Radial Network for comparison
Points_exact = load('AC_conf.mat');
Points_exact = Points_exact.sortedPoints;
lightzs = [181, 180, 214]/255;
h2 = fill(Points_exact(:,1), Points_exact(:,2), lightzs);
set(h2, 'facealpha', 0.9, 'EdgeColor', 'none');
%% region considering uncertainty
% load("uncertainty_mesh.mat")
% Points_un = sortedPoints;
% lightred = [240, 128, 122]/255;
% h3 = fill(Points_un(:,1), Points_un(:,2), lightred);
% set(h3, 'facealpha', 0.6, 'EdgeColor', 'none');
%% compare the power flow
% red = [240, 69, 49]/255;
% h3=plot(-0.85,1.54,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', red, 'MarkerFaceColor', red);
% p2 = [84, 84, 84]/255;
% h4=plot(-0.981411034156963,1.86212451429350,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', p2, 'MarkerFaceColor', p2);
%% set the range
xlim([-2 2]);
ylim([-2.2 2.2]);
set(gca, 'FontSize', 14,'FontName', 'Times New Roman');

x_label = xlabel('$p^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
y_label = ylabel('$q^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');