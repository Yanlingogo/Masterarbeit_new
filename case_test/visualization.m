function  visualization(result_projection, result_grid, result_optimization, visualization_options)
    fig=figure; box on; grid on; hold on; set(fig, 'Position', [100, 100, 650, 550]);
    % load the color parameters
    lightorange = [250, 188, 113]/255;
    orange = [255, 99, 0]/255;
    light_violett = [181, 180, 214]/255;
    violett = [158, 49, 252]/255;
    light_green = [27, 214, 108]/255;
    green = [36, 138, 36]/255;

    fill_optimization = [222, 235, 176;
                         181, 180, 214]/255;
    line_optimization = [36, 138, 36;
                         123, 200, 91]/255;
    % Projection method
    if sum(visualization_options=='P')
        h1 = fill(result_projection(:,1), result_projection(:,2), lightorange);
        set(h1, 'facealpha', 0.7, 'EdgeColor', orange,'LineWidth',2);
    end
    % Grid method
    h2 = gobjects(size(result_grid, 1), 1);
    if sum(visualization_options=='G')
        for i = size(result_grid,1)
            fill_gradient = light_violett*(i/size(result_grid,1));
            line_gradient = violett*(i/size(result_grid,1));
            h2(i)  = fill(result_grid{i}(:,1), result_grid{i}(:,2),fill_gradient);
            set(h2(i), 'facealpha', 0.7, 'EdgeColor', line_gradient,'LineWidth', 2);
        end
    end
    % Optimization method
    h3 = gobjects(size(result_optimization, 2), 1);
    if sum(visualization_options=='O')
        for i = 1: size(result_optimization,2)
            h3(i) = fill(result_optimization{i}(:,1), result_optimization{i}(:,2),fill_optimization(i,:));
            set(h3(i), 'facealpha', 0.7, 'EdgeColor', line_optimization(i,:),'LineWidth', 2);
        end        
    end
    validHandles = [h1(1), h2(1), h3(1)];
    validHandles = validHandles(ishandle(validHandles)); % Filter out invalid handles

    if ~isempty(validHandles) 
        legend(validHandles, {'Projection', 'Grid', 'Optimization'});
    else
        warning('没有有效的图形句柄用于创建图例。');
    end


    x_label = xlabel('$p^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
    set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
    
    y_label = ylabel('$q^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
    set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');

end



