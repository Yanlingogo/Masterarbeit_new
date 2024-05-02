clc;
clear;
%% Load case
% m.file 
mpc = loadcase('case33_org.m');
% case modification
Nbus = size(mpc.bus,1);
Ngen = size(mpc.gen,1);
Nbranch = size(mpc.branch,1);
if Nbus ~= Nbranch+1
    error('The grid is non-radial');
end
if Ngen == 1
    warning('No flexibility currently, DERs will automaticlly be added')
    % set the number of DER, and ratio of sum(power)/sum(load)
    Num_DER = 4;
    ratio= 1; %(0.25, 0.5, 1, 2, 4)
    % function 
    mpc = case_modification(mpc,Num_DER,ratio);
end

%% AC model: flexibility aggregation
mpc_AC = mpc;
% more sample point, higher precision
Sampling_resolution = 30;
[AC_points, AC_time]= AC_flexibility_aggregation(mpc_AC,Sampling_resolution);
%% DC model: flexibility aggregation

%% LinDistFlow model: flexibility aggregation, radial grid only!
mpc_LinDistFlow = mpc;
% sample points for vertexes
[LinDistFlow_points, Lin_dispatch,LinDistFlow_time] = LinDistFlow_flexibility_aggregation(mpc_LinDistFlow);

%% Modification for LinDistFlow model
% Jacobian at base point 
data_basepoint = jacobian_mpc(mpc_LinDistFlow);
% model: 2-fixed current, 3-linearized current
model = 2;
[vertexes_esti, dispatch_esti, Compensation_time] = LinDistFlow_flexibility_aggregation(mpc_LinDistFlow,model,data_basepoint);

%% Compensation methods
% choose compensation method
% 1: compensation for LinDistFlow model
% 2: compensation for modified DistFlow model
method = 2;
if method == 1
    vertexes = LinDistFlow_points;
    dispatch = Lin_dispatch;
else
    vertexes = vertexes_esti;
    dispatch = dispatch_esti;
end

% voltage_check: 0-off, 1-on
% compensation_type: 1-quadratic determined current, 2-corrected PCC
% resolution: segmentation of a boundary line
voltage_check = 1;
comp_type     = 2;
resolution    = 30;
line_points = compensation_line(vertexes,dispatch,data_basepoint,...
    voltage_check,comp_type,resolution,method);

%% Visualization 
visualization(AC_points,LinDistFlow_points,line_points)