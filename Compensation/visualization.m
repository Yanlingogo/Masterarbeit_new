function visualization(AC_points,LinDistFlow_points,line_points)
    % set the position of the figure
    fig=figure; box on; grid on; hold on; set(fig, 'Position', [100, 100, 650, 550])
    xlim([-1.5 2]);
    ylim([-1.5 2]);
    set(gca, 'FontSize', 14,'FontName', 'Times New Roman');
    
    % AC model
    lightzs = [181, 180, 214]/255;
    violett = [158, 49, 252]/255;
    h1 = fill(AC_points(:,1), AC_points(:,2), lightzs);
    set(h1, 'facealpha', 0.9, 'EdgeColor', violett,'LineWidth',2);
    hold all
    % LinDistFlow 
    lightorange = [250, 188, 113]/255;
    orange = [255, 99, 0]/255;
    h2 = fill(LinDistFlow_points(:,1), LinDistFlow_points(:,2), lightorange, 'FaceAlpha',0.7);
    set(h2, 'facealpha', 0.6, 'EdgeColor', orange,'LineWidth', 2);
    h3 = scatter(LinDistFlow_points(:,1), LinDistFlow_points(:,2), 25, 'filled',...
         'MarkerFaceColor', orange,'MarkerEdgeColor', 'none');
    % compensated result
    for i = 1:length(line_points)
        data = line_points{i}; % 获取当前 cell 元素的数据
        
        % 检查数据尺寸
        if size(data, 2) ~= 2 || size(data, 1) < 2
            error('Each cell element must contain an Nx2 matrix.');
        end
        
        x = data(:, 1); % x 坐标
        y = data(:, 2); % y 坐标
        
        plot(x, y, 'LineWidth', 2); % 绘制线条
        % scatter(x,y);
    end
    x_label = xlabel('$p^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
    set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
    y_label = ylabel('$q^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
    set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
end



