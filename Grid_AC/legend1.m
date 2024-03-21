
figure;
red = [240, 69, 49]/255;
% h3=plot(-0.85,1.54,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', red, 'MarkerFaceColor', red);
h3=plot(NaN,NaN,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', red, 'MarkerFaceColor', red);
p2 = [84, 84, 84]/255;
% h4=plot(-0.981411034156963,1.86212451429350,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', p2, 'MarkerFaceColor', p2);
h4=plot(NaN,NaN,'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', p2, 'MarkerFaceColor', p2);
unique_colors = [0 0.447 0.7410;
                 0.929 0.694 0.125;
                 [120, 171, 50]/255;
                 [122, 50, 139]/255;
                 lightorange
                 lightzs;
                 lightbl;
                 ];
hold on;
lineWidths = [2*ones(4,1); 5*ones(3,1);];
h_fake = arrayfun(@(x) plot(NaN, NaN, 'Color', unique_colors(x,:), 'LineWidth', lineWidths(x)), 1:size(unique_colors, 1));
h_sum = [h_fake, h3, h4];
legend(h_sum, {'V bounds', 'P_{gen} bounds','Q_{gen} bounds', 'S_{line} bounds','LinDistFlow', 'Radial Network','Mesh Network','A','B'},'NumColumns', 9,'FontSize',12);
axis off
exportgraphics(gca, 'legendOnly.pdf', 'ContentType', 'vector');