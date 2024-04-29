function middlePoint = findMidpoint(A, B, C)
    % 将输入点转换为向量
    AB = B - A; % 向量AB
    AC = C - A; % 向量AC
    BA = A - B; % 向量BA
    BC = C - B; % 向量BC

    % 计算标量积
    dotProduct_AB_AC = dot(AB, AC);
    dotProduct_BA_BC = dot(BA, BC);

    % 初始化中间点为空
    middlePoint = [];

    % 检查C是否在A和B之间
    if dotProduct_AB_AC >= 0 && dotProduct_BA_BC >= 0
        middlePoint = C;
    % 检查A是否在B和C之间
    elseif dot(AC, -BC) >= 0 && dot(-AC, AB) >= 0
        middlePoint = A;
    % 检查B是否在A和C之间
    elseif dot(BC, BA) >= 0 && dot(-BA, -AC) >= 0
        middlePoint = B;
    end

    % 如果找到了中间点，显示它
    if ~isempty(middlePoint)
        disp('Middle point is:');
        disp(middlePoint);
    else
        disp('No middle point found.');
    end
end

