function [gridPointsInside] = grid_region(vert)
    
    % 步骤1: 找到最小矩形的边界
    minX = min(vert(:,1))+1e-5;
    maxX = max(vert(:,1))-1e-5;
    minY = min(vert(:,2))+1e-5;
    maxY = max(vert(:,2))-1e-5;
    
    % 步骤2: 创建网格
    N = 50;
    xGrid = linspace(minX, maxX, N);
    yGrid = linspace(minY, maxY, N);

    % 生成网格
    [XGrid, YGrid] = meshgrid(xGrid, yGrid);
    xGrid = XGrid(:);
    yGrid = YGrid(:);
    
    % 步骤3: 去除区域外的点
    [in, ~] = inpolygon(xGrid, yGrid, vert(:,1), vert(:,2));
    xGridInside = xGrid(in);
    yGridInside = yGrid(in);
    
    gridPointsInside = [xGridInside, yGridInside];
    
end

