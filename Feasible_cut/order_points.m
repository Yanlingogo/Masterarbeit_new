function [vert] = order_points(vert)
    centroid = mean(vert, 1); % surely convex
    angles = atan2(vert(:,2) - centroid(2), vert(:,1) - centroid(1));
    [~, order] = sort(angles);
    vert = vert(order, :);
end

