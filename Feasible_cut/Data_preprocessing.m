clc;
clear;
T = 6;

WS = 15*rand(2,T);
WS_stat = cell(1,2);
WS_stat{1} = zeros(500,T);
WS_stat{2} = zeros(500,T);
PV_power = 1*rand(2,T);
Lp = zeros(32,T);
Lq = zeros(32,T);
for i = 1:T
    Lp(:,i) = 0.8*rand*rand(32,1);
    Lq(:,i) = 1*rand*rand(32,1);
end
%% distribution function for Wind speed
b = 1.5;
a = WS*(b/(b-1))^(1/b);
for j = 1:2
 for i = 1:T
    WS_stat{j}(:,i) = wblrnd(a(j,i), b, 500, 1); 
 end
end
% parameter for Wind turbine
vi = 4;
vr = 10;
v0 = 25;
Pr = 1; % p.u. 
WT_power = func_speed_power(WS,vi,vr,v0,Pr);
WT_stat =cell(1,2);
for i = 1:2
    WT_stat{i} = func_speed_power(WS_stat{i},vi,vr,v0,Pr);
end
%% distribution function for PV
sigma_PV = 0.1*PV_power;
PV_stat = cell(1,2);
PV_stat{1} = zeros(500,T);
PV_stat{2} = zeros(500,T);
for j = 1:2
 for i = 1:T
    PV_stat{j}(:,i) = PV_power(j,i)+sigma_PV(j,i)*randn(500, 1); 
 end
end
%% distribution function for Loads
sigma_loadp = 0.2*Lp;
sigma_loadq = 0.2*Lq;
Load_stat = cell(2,32);
for j = 1:32
    for i = 1:T
        Load_stat{1,j}(:,i) = Lp(j,i)+sigma_loadp(j,i)*randn(500,1);
        Load_stat{2,j}(:,i) = Lq(j,i)+sigma_loadq(j,i)*randn(500,1);
    end
end
%% covariance matrix at each time periode

for i =1:T
    cov_all{1,i} = cov(WT_stat{1}(:,i)',WT_stat{2}(:,i)');
    cov_all{2,i} = cov(PV_stat{1}(:,i)',WT_stat{2}(:,i)');
    Loadp = [];
    Loadq = [];
    for j = 1:32
        Loadp = [Loadp Load_stat{1,j}(:,i)];
        Loadq = [Loadq Load_stat{2,j}(:,i)];
    end
    cov_all{3,i} = cov(Loadp);
    cov_all{4,i} = cov(Loadq);

end

%% Cholesky decomposition
cov_de = cell(4,T);
for i = 1:4 
    for j = 1:T
        if all(eig(cov_all{i,j})>0)
            cov_de{i,j} = chol(cov_all{i,j},'upper');
        else
            fprintf('cov_all{%d,%d} is not positive definite\n',i,j)
        end
    end
end

%% convert to p.u.
baseMVA = 10;
WT_power = WT_power/baseMVA;
PV_power = PV_power/baseMVA;
Lp = Lp/baseMVA;
Lq = Lq/baseMVA;
scaled_cov = cellfun(@(x) x/baseMVA, cov_de, 'UniformOutput', false);
clearvars -except Lp Lq PV_power WT_power scaled_cov