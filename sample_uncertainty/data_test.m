numComponents = 4; % 选择高斯分布的数量
options = statset('MaxIter',2000); % 可选，设置最大迭代次数以确保收敛

% 拟合高斯混合模型
gmmModel = fitgmdist(e_wind, numComponents, 'Options', options);

% 查看AIC
aic = gmmModel.AIC;

% 查看BIC
bic = gmmModel.BIC;

figure; % 创建新图形
histogram(e_wind, 'Normalization', 'probability','BinWidth',0.1); % 绘制并归一化以匹配PDF
hold on; % 保持图像，以便在同一图中添加更多图层

% 生成评估GMM PDF的点
xRange = linspace(min(e_wind), max(e_wind), 1000); % 生成一个覆盖误差数据范围的点序列
gmmPdf = pdf(gmmModel, xRange'); % 计算GMM的PDF

% 绘制GMM的PDF
plot(xRange, gmmPdf, 'r--', 'LineWidth', 2); % 绘制GMM的PDF
legend('原始数据PDF', 'GMM拟合PDF');
xlabel('误差值');
ylabel('概率密度');
title('误差数据与GMM拟合结果的比较');
hold off; % 解除保持状态

means = gmmModel.mu;
covariances = gmmModel.Sigma;

e_wind1 = (wind1(:,1)-wind1(:,2))./wind1(:,2);
for box region
load("e_wind1.mat");
load("e_wind2.mat");
error = [e_wind1,e_wind2]';
mu = mean(error,2);
cov_wind = cov(error');

[V,D] = eig(cov_wind);
D_inv_sqrt = D.^(-1/2);
Sigma_inv_sqrt = V * D_inv_sqrt * V';
theta = Sigma_inv_sqrt * (error - mu);

theta_b = 10;

lambda = sdpvar()