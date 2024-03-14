% 假设 coords1 和 coords2 分别存储两组坐标点
load("mesh_R7.mat"); coords1 = sortedPoints;
load("AC_mesh4.1086.mat"); coords2 = sortedPoints;
x1 = coords1(:, 1);
y1 = coords1(:, 2);
x2 = coords2(:, 1);
y2 = coords2(:, 2);

% 计算两个多边形的面积
area1 = polyarea(x1, y1);
area2 = polyarea(x2, y2);

% 计算两个多边形的交集区域（需要Mapping Toolbox）
[xi, yi] = polybool('intersection', x1, y1, x2, y2);

% 计算交集区域的面积
areaIntersect = polyarea(xi, yi);

% 计算不相交部分的面积
areaNonIntersect = (area1 + area2) - 2 * areaIntersect;

% 显示结果
disp(['The non-intersecting area is: ', num2str(areaNonIntersect)]);
