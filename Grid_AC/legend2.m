figure;
red = [240, 69, 49]/255;
h3=plot(-0.85,1.54,'o', 'MarkerSize', 7, 'MarkerEdgeColor', red, 'MarkerFaceColor', red);
% h3=plot(NaN,NaN,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', red, 'MarkerFaceColor', red);
p2 = [84, 84, 84]/255;
h4=plot(-0.98,1.86,'o', 'MarkerSize', 7, 'MarkerEdgeColor', p2, 'MarkerFaceColor', p2);
% h4=plot(NaN,NaN,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', p2, 'MarkerFaceColor', p2);
legend([h3, h4], {'Red Line', 'Blue Line'}, 'Location', 'best');
axis off; % 关闭坐标轴的显示
