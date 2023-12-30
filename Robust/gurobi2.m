import casadi.*

% Define your optimization problem using CasADi symbols
x = MX.sym('x', 2);
f = x' * x;  % Objective function
g = x;      % Constraints
nlp = struct('x', x, 'f', f, 'g', g);

% Set solver options
opts = struct();
opts.verbose = true;
opts.print_time = true;

% Create a QP solver instance
solver = qpsol('solver', 'qpoases', nlp, opts);

% Solve the problem
lbg = [0; 0];
ubg = [1; 1];
sol = solver('lbg', lbg, 'ubg', ubg);

disp(sol)
