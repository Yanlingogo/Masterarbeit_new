%% solver configuration
% choose the solver P:projection, G:Grid, O:Optimziation
choose_solver = 'PO';
resolution_G = 100;
resolution_O = [0.5 0.05];
% Visualization options 
visualization_options = 'PO';
%% load the case
mpc = loadcase('case1888rte.m');
mpc = ext2int(mpc);
Nbus = size(mpc.bus,1);
Ngen = size(mpc.gen,1);
Nbranch = size(mpc.branch,1);
if Nbus ~= Nbranch+1
    error('The grid is non-radial');
end
if Ngen == 1
    error('No flexibility')
end
%% modify the case file
% mpc = create_case_modification(mpc);

%% solver
% Projection-based Method(linear, radial)
result_projection = [];
time_projection = [];
if sum(choose_solver=='P')
    [result_projection, time_projection]= solver_projection(mpc);
end
% check the center of the AF
plot(result_projection(:,1),result_projection(:,2));
% Grid-Method
result_grid = cell(size(resolution_G,2),10);
time_grid = zeros(size(resolution_G,2));
if sum(choose_solver=='G')
    [result_grid{i,:}, time_grid(i)]= solver_grid(mpc,resolution_G(i));
end
% Optimization-based Method
result_optimization = cell(size(resolution_O));
time_optimization = zeros(size(resolution_O,2),1);
for i = 1: size(resolution_O,2)
    [result_optimization{i}, time_optimization(i)] = solver_optimization(mpc,resolution_O(i));
end
%% Visualization
visualization(result_projection, result_grid, result_optimization, visualization_options);

%% Area/error factor calculation
[area, error_factor] = area_error(result_projection,result_optimization);

%print the result 
time = [time_projection;time_optimization];
numColumns = size(error_factor,1);
for i = 1:numColumns
    % 打印每列对应的面积和误差
    fprintf('method%d: area: %.2f, error factor: %.2f\n', i, area(i), error_factor(i));
end
