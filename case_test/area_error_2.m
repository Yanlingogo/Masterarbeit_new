function [area_all, error_factor] = area_error(result_projection, result_optimization)
    %% second cell of result_optimization is the benchmark
    coords1 = result_projection;
    coords2 = result_optimization{1};
    benchmark = result_optimization{2};
    
    % 初始化面积值
    area1 = NaN;
    area2 = NaN;
    area_benchmark = NaN;
    
    % 检查输入数据是否为空，并计算非空数据的面积
    if ~isempty(coords1)
        ps1 = polyshape(coords1(:, 1), coords1(:, 2));
        area1 = area(ps1);
    end
    if ~isempty(coords2)
        ps2 = polyshape(coords2(:, 1), coords2(:, 2));
        area2 = area(ps2);
    end
    if ~isempty(benchmark)
        ps_benchmark = polyshape(benchmark(:, 1), benchmark(:, 2));
        area_benchmark = area(ps_benchmark);
    end
    
    % 初始化交集和不相交部分的面积
    areaIntersect_projection = NaN;
    areaIntersect_optimization = NaN;
    areaNonIntersect_projection = NaN;
    areaNonIntersect_optimization = NaN;
    
    % 如果基准数据非空，计算与其相交的面积
    if ~isnan(area_benchmark)
        if ~isnan(area1)
            ps_projection_intersection = intersect(ps1, ps_benchmark);
            areaIntersect_projection = area(ps_projection_intersection);
            areaNonIntersect_projection = (area1 + area_benchmark) - 2 * areaIntersect_projection;
        end
        if ~isnan(area2)
            ps_optimization_intersection = intersect(ps2, ps_benchmark);
            areaIntersect_optimization = area(ps_optimization_intersection);
            areaNonIntersect_optimization = (area2 + area_benchmark) - 2 * areaIntersect_optimization;
        end
    end

    % 组装输出
    area_all = [area1, area2, area_benchmark];
    % 计算误差因子，处理 NaN 的情况
    error_factor_projection = 0;
    error_factor_optimization = 0;
    if ~isnan(areaNonIntersect_projection) && ~isnan(area_benchmark) && area_benchmark ~= 0
        error_factor_projection = (areaNonIntersect_projection / area_benchmark) * 100;
    end
    if ~isnan(areaNonIntersect_optimization) && ~isnan(area_benchmark) && area_benchmark ~= 0
        error_factor_optimization = (areaNonIntersect_optimization / area_benchmark) * 100;
    end
    error_factor = [error_factor_projection, error_factor_optimization];
end
