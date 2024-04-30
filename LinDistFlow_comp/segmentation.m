% function Points = segmentation(vertexes, num_interp)
% 
%     % Add the first vertex to the end of the vertex list to form a closed loop
%     vertexes(end+1, :) = vertexes(1, :);
% 
%     % Initialize the cell array to store the point sets for each segment.
%     Points = cell(size(vertexes, 1)-1, 1);
% 
%     % Iterate through all pairs of points (now including from the last point back to the first)
%     for i = 1:size(vertexes, 1)-1
%         % Take the current point and the next point
%         p1 = vertexes(i, :);
%         p2 = vertexes(i+1, :);
% 
%         % Each dimension is interpolated separately
%         x = linspace(p1(1), p2(1), num_interp + 2);
%         y = linspace(p1(2), p2(2), num_interp + 2);
% 
%         % Add interpolated points (including start and end points) to the cell
%         interp_points = [x' y'];
% 
%         % Save to the corresponding cell
%         Points{i} = interp_points;
%     end
% end
function [gridPointsInside] = segmentation(vert,N)
    
    % 步骤1: 找到最小矩形的边界
    minX = min(vert(:,1))+1e-5;
    maxX = max(vert(:,1))-1e-5;
    minY = min(vert(:,2))+1e-5;
    maxY = max(vert(:,2))-1e-5;
    
    % 步骤2: 创建网格
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
    