function [area_all, error_factor] = area_error(result_projection,result_optimization)
    %% second cell of result_optimization is the benchmark
    coords1 = result_projection;
    coords2 = result_optimization{1};
    benchmark = result_optimization{2};
    x1 = coords1(:, 1);
    y1 = coords1(:, 2);
    x2 = coords2(:, 1);
    y2 = coords2(:, 2);
    x_benchmark = benchmark(:,1);
    y_benchmark = benchmark(:,2);
    
    % 使用polyshape创建多边形对象
    ps1 = polyshape(x1, y1);
    ps2 = polyshape(x2, y2);
    ps_benchmark = polyshape(x_benchmark, y_benchmark);
    
    % 计算各个多边形的面积
    area1 = area(ps1);
    area2 = area(ps2);
    area_benchmark = area(ps_benchmark);
    
    % 计算两个多边形与基准多边形的交集区域
    ps_projection_intersection = intersect(ps1, ps_benchmark);
    ps_optimization_intersection = intersect(ps2, ps_benchmark);
    
    % 计算交集区域的面积
    areaIntersect_projection = area(ps_projection_intersection);
    areaIntersect_optimization = area(ps_optimization_intersection);
    
    % 计算不相交部分的面积
    areaNonIntersect_projection = (area1 + area_benchmark) - 2 * areaIntersect_projection;
    areaNonIntersect_optimization = (area2 + area_benchmark) - 2 * areaIntersect_optimization;

    area_all = [area1 area2 area_benchmark];
    error_factor = [areaNonIntersect_projection areaNonIntersect_optimization]/area_benchmark*100;
end

