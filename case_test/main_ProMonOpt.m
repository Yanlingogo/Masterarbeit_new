clc
clear

%% load the case
mpc = loadcase('case141.m');
mpc = ext2int(mpc);
Nbus = size(mpc.bus,1);
Ngen = size(mpc.gen,1);
Nbranch = size(mpc.branch,1);
if Nbus ~= Nbranch+1
    error('The grid is non-radial');
end
if Ngen == 1
    warning('No flexibility currently, DERs will automaticlly be added')
    mpc = create_case_modification(mpc);
end

%% solver
% projection
[result_projection, time_projection]= vertex_ipopt(mpc);
% monte carlo
[result_monte, time_monte] = Monte_solver(mpc);
% optimization
Pd = sum(mpc.bus(:,3))/mpc.baseMVA;
Qd = sum(mpc.bus(:,4))/mpc.baseMVA;
% resolution = [sqrt(Pd^2+Qd^2)/10,sqrt(Pd^2+Qd^2)/100];
resolution = [5,30];
result_optimization = cell(size(resolution,2));
time_optimization = zeros(size(resolution,2),1);
for i = 1: size(resolution,2)
    [result_optimization{i}, time_optimization(i)] = sample_AC(mpc,resolution(i));
end
%% Visualization
hold all
plot(result_projection(:,1),result_projection(:,2));
plot(result_monte(:,1),result_monte(:,2));
plot(result_optimization{1}(:,1),result_optimization{1}(:,2));
%% performance evaluation
% LinDistFlow model
[Lost_Lin, risk_Lin, area_Lin, area_ref] = compensation_error(result_projection,result_optimization{2});
% Monte carlo
[Lost_Monte, risk_Monte, area_Monte, ~] = compensation_error(result_monte,result_optimization{2});
% Optimization-based (angle-based)
[Lost_Angle, risk_Angle, area_Angle,~] = compensation_error(result_optimization{1},result_optimization{2});

result = [area_Lin,area_Monte,area_Angle,area_ref,Lost_Lin,Lost_Monte,Lost_Angle,...
    risk_Lin,risk_Monte,risk_Angle,time_projection,time_monte,time_optimization(1)];