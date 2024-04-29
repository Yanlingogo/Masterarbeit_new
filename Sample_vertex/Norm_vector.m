function normals = Norm_vector(points)
    % 确保闭合的线段序列
    if points(1,:) ~= points(end,:)
        points = [points; points(1,:)]; % 添加第一个点到末尾以闭合区域
    end
    
    % 初始化法向量存储矩阵
    normals = zeros(size(points,1)-1, 2);
    
    for i = 1:size(points,1)-1
        % 当前线段的方向向量
        d = points(i+1, :) - points(i, :);
        
        % normal vector to the outer side
        n = [d(2), -d(1)];
        
        % affine hull of facet
        n = n / abs((n*points(i,:)'));
        
        normals(i, :) = n;
    end
end

