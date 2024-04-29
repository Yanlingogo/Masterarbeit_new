clc
clear
close all
%% load the case
mpc = loadcase('case33_org.m');
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
%% solver Optimization-based Method 

% adaptive angle-based
Pd = sum(mpc.bus(:,3))/mpc.baseMVA;
Qd = sum(mpc.bus(:,4))/mpc.baseMVA;
resolution = sqrt(Pd^2+Qd^2)/10;
[result_optimization, time_AC] = solver_optimization(mpc,resolution);
% normal direction
[vertexes, dispatch, time_LinD] = vertex_ipopt(mpc);
% Jacobian at z_0
z_0 = mean(vertexes);
% z_0 = [0 0];
data_basepoint = jacobian_mpc(mpc,z_0);
% compensation
Resoultion  = 50;
line_points = cell(size(vertexes,1),1);
line_points_bounds = cell(size(vertexes,1),1);
numlines = size(vertexes,1);
for i = 1: numlines
    power_dispatch = [dispatch{i} dispatch{mod(i,numlines)+1}];
    endpoints = [vertexes(i,:)' vertexes(mod(i,numlines)+1,:)'];
    line_points{i} = compensation_line_cor(data_basepoint, endpoints, power_dispatch, Resoultion);
end


%% visualization
fig=figure; box on; grid on; hold on; set(fig, 'Position', [100, 100, 650, 550])
% AC model
lightzs = [181, 180, 214]/255;
violett = [158, 49, 252]/255;
h1 = fill(result_optimization(:,1), result_optimization(:,2), lightzs);
set(h1, 'facealpha', 0.9, 'EdgeColor', violett,'LineWidth',2);
hold all
% LinDistFlow 
lightorange = [250, 188, 113]/255;
orange = [255, 99, 0]/255;
h2 = fill(vertexes(:,1), vertexes(:,2), lightorange, 'FaceAlpha',0.7);
set(h2, 'facealpha', 0.6, 'EdgeColor', orange,'LineWidth', 2);
% compensation result

visualization_lines(line_points)


[lost_Lin, risk_Lin, area_Lin, area_ref] = compensation_error(vertexes,result_optimization);
% for compensation
line_points_all = vertcat(line_points{:,:});
[lost_Comp, risk_Comp, area_Comp, ~] = compensation_error(line_points_all,result_optimization);

result = [area_Lin,area_Comp,area_ref,lost_Lin,lost_Comp,risk_Lin,risk_Comp];