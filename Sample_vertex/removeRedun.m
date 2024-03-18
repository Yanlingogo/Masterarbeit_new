function [filteredPoints] = removeRedun(points, tolerance)
    % 计算凸包
    k = convhull(points(:,1), points(:,2));
    hullPoints = points(k, :);
    
    % 初始化一个逻辑向量，标记所有点最初都不是冗余的
    nRedundant = true(size(points, 1), 1);
    
    % 对于每个点，检查它是否在凸包的边上（考虑误差）
    for i = 1:size(points, 1)
        point = points(i, :);
        otherpoints = unique_normal(point,points,1e-6);
        otherpoints = [otherpoints;otherpoints(1,:)];
        for j = 1:size(otherpoints,1)
            dist = vectorizedDistanceToLines(otherpoints,point);
            if any(dist <= 1e-6)
                idx = find(dist <= 1e-6);
                p1 = otherpoints(idx,:);
                p2 = otherpoints(idx+1,:);
                nRedundant(i) = ~isPointBetween(point,p1,p2);
            end
        end
    end
    
    % 凸包上的点不考虑为冗余
    % isRedundant(k) = false;
    
    % 过滤掉标记为冗余的点
    filteredPoints = points(nRedundant, :);
end

function distances = vectorizedDistanceToLines(points, singlePoint)
    % 确保points是一个闭环
    if ~isequal(points(1,:), points(end,:))
        points(end+1,:) = points(1,:);
    end
    
    % 计算所有线段的方向向量
    directionVectors = diff(points, 1, 1);
    
    % 计算点到每个线段起点的向量
    pointVectors = singlePoint - points(1:end-1, :);
    
    % 利用叉积公式计算垂直距离
    crossProd = cross([directionVectors zeros(size(directionVectors, 1), 1)], [pointVectors zeros(size(pointVectors, 1), 1)], 2);
    distances = abs(crossProd(:,3)) ./ vecnorm(directionVectors, 2, 2);
end

function isBetween = isPointBetween(point, lineStart, lineEnd)
    % 检查点是否位于线段的两个端点之间
    dot1 = dot(point - lineStart, lineEnd - lineStart);
    dot2 = dot(point - lineEnd, lineStart - lineEnd);
    isBetween = dot1 > 0 && dot2 > 0;
end
