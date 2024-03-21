function [] = plot_bus(exact_mpc,plot_bus,plot_rng,resolution)
%% Get system data
mpc=exact_mpc;
% mpc=runopf(mpc,mpoption('verbose',0,'out.all',0));
% mpc.gen(:,end+1)=[1;0;0];

%% get plot data
mesh_intensity=resolution; 
mesh_axis=plot_rng;
[u1_plot,u2_plot,mesh_all,mesh_feasibility,u_base]=get_plot_data(mpc,plot_bus,mesh_axis,mesh_intensity);

%% Plot figure
fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 850, 650]);
clim([0 400]); colormap([0.7461 0.832 0.9062 ;0.875*0.9 0.9492*0.9 0.6875*0.9; 0.875 0.9492 0.6875; 1 1 1]);
tol = 0;
contourf(u1_plot,u2_plot,mesh_feasibility,[-tol tol])
% 蓝色-电压，黄色-PQ, 紫色-线路功率，深红-线路辐角
pcolor_meshall=[0.25 0.25 0.25; ones(2,1)*[0 0.447 0.7410]; ones(4,1)*[0.929 0.694 0.125]; 0.494 0.184 0.556; ones(2,1)*[0.635 0.078 0.184]];
ptype_meshall={'-','-','--','-','--','-','--','-','--','-'};
for i=1:size(mesh_all,2)
    if size(mesh_all{i},1)~=1
        for j=1:size(mesh_all{i},3)
            if max(max(mesh_all{i}(:,:,j)))>=0 
                contour(u1_plot,u2_plot,mesh_all{i}(:,:,j),[-tol tol],ptype_meshall{i},'color',pcolor_meshall(i,:),'LineWidth',2);
            end
        end
    end
end

% for i = 1:size(mesh_all, 2)
%     if size(mesh_all{i}, 1) ~= 1
%         for j = 1:size(mesh_all{i}, 3)
%             if max(max(mesh_all{i}(:,:,j))) >= 0 % 满足约束才能画图
%                 hsv = rgb2hsv(pcolor_meshall(i, :));
%                 hsv(3) = min(hsv(3) * (1.05^j), 1); % 确保亮度不超过 1
%                 new_color = hsv2rgb(hsv);
%                 contour(u1_plot, u2_plot, mesh_all{i}(:,:,j), [-tol tol], ptype_meshall{i}, 'Color', new_color, 'LineWidth', 2);
%             end
%         end
%     end
% end

[x_p,y_q] =find(mesh_feasibility>=0);
u1_x = u1_plot(1,:);u2_y = u2_plot(:,1)';
P = u1_x(x_p); Q = u2_y(y_q);
C = contourc(u1_x, u2_y, mesh_feasibility, [0 0]);
idx = 1;
while idx < size(C, 2)
    level = C(1, idx);  % 等高线的高度值
    nPoints = C(2, idx);  % 等高线上的点数
    points = C(:, idx+1 : idx+nPoints);  % 等高线的坐标点
    % 处理 points...
    idx = idx + 1 + nPoints;  % 移动到下一个等高线段
end
% points = unique(points', 'rows', 'stable'); % delete the same points
% Ei = max_ei(points);

%% draw the max rectangle
matrix = (mesh_feasibility>=-tol);
[minC, maxC, minR, maxR, maxArea] = max_rectangle(matrix);
P = u1_x([minR,maxR]); Q = u2_y([minC,maxC]);
vertex = [P(1),P(1),P(2),P(2);
          Q(1),Q(2),Q(2),Q(1)];
rec = fill(vertex(1,:),vertex(2,:),[1, 0.6, 0.6]);


%plot(Ei(1,:), Ei(2,:), 'g', 'LineWidth', 2);

% 
x_label=xlabel(['$p_{pcc,' num2str(plot_bus) '}$ (p.u.)']); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
y_label=ylabel(['$q_{pcc,' num2str(plot_bus) '}$ (p.u.)']); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
title(['Feasible Region of Bus ' num2str(plot_bus) ' '],'FontSize',15,'FontName','Times New Roman');
set(gca,'fontsize',15,'FontName','Times New Roman')
text(2.0, 0.5, 'Exact Feasible Region', 'Color', 'k', 'FontSize', 12,'FontName','Times New Roman');
text(0.5, 0.5, 'hyperbox approximation', 'Color', 'k', 'FontSize', 12,'FontName','Times New Roman')
% add explanation for the bounds
unique_colors = [0 0.447 0.7410;
                 0.929 0.694 0.125;
                 0.494 0.184 0.556;
                 0.635 0.078 0.184;
                 ];
hold on;
h_fake = arrayfun(@(x) plot(NaN, NaN, 'Color', unique_colors(x,:), 'LineWidth', 2), 1:size(unique_colors, 1));
legend(h_fake, {'V bounds', 'P_{gen}/Q_{gen} bounds', 'S_{line} bounds', '\phi_{line} bounds'});