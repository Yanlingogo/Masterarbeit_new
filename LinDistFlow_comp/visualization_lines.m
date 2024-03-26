function visualization_lines(line_data)
    for i = 1:length(line_data)
        data = line_data{i}; % 获取当前 cell 元素的数据
        
        % 检查数据尺寸
        if size(data, 2) ~= 2 || size(data, 1) < 2
            error('Each cell element must contain an Nx2 matrix.');
        end
        
        x = data(:, 1); % x 坐标
        y = data(:, 2); % y 坐标
        
        plot(x, y, 'LineWidth', 2); % 绘制线条
    end
end

