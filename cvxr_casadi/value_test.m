%% test for operation point
% min(g_pu)
g_pu = zeros(9,1);
for i = 1:9
    obj_g = @(x) x(entries{9}(i));
    constraint = g(x_SX);
    objective = obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_pu(i) = full(sol.f);
end
% max(g_pl)
g_pl = zeros(9,1);
for i = 1:9
    obj_g = @(x) x(entries{10}(i));
    constraint = g(x_SX);
    objective = -obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_pl(i) = -full(sol.f);
end
% min(g_qu)
g_qu = zeros(9,1);
for i = 1:9
    obj_g = @(x) x(entries{11}(i));
    constraint = g(x_SX);
    objective = obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_qu(i) = full(sol.f);
end
% max(g_ql)
g_ql = zeros(9,1);
for i = 1:9
    obj_g = @(x) x(entries{12}(i));
    constraint = g(x_SX);
    objective = -obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_ql(i) = -full(sol.f);
end
% min(g_vvcos_u)
g_cu = zeros(nl,1);
for i = 1:nl
    obj_g = @(x) x(entries{13}(i));
    constraint = g(x_SX);
    objective = obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_cu(i) = full(sol.f);
end
% max(g_vvcos_l)
g_cl = zeros(nl,1);
for i = 1:nl
    obj_g = @(x) x(entries{14}(i));
    constraint = g(x_SX);
    objective = -obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_cl(i) = -full(sol.f);
end
% min(g_vvsin_u)
g_su = zeros(nl,1);
for i = 1:nl
    obj_g = @(x) x(entries{15}(i));
    constraint = g(x_SX);
    objective = obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_su(i) = full(sol.f);
end
% max(g_vvsin_l)
g_sl = zeros(nl,1);
for i = 1:nl
    obj_g = @(x) x(entries{16}(i));
    constraint = g(x_SX);
    objective = -obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_sl(i) = -full(sol.f);
end
% min(g_vv_u)
g_vvu = zeros(nb,1);
for i = 1:nb
    obj_g = @(x) x(entries{17}(i));
    constraint = g(x_SX);
    objective = obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_vvu(i) = full(sol.f);
end
% max(g_vv_l)
g_vvl = zeros(nb,1);
for i = 1:nb
    obj_g = @(x) x(entries{18}(i));
    constraint = g(x_SX);
    objective = -obj_g(x_SX);
    nlp = struct('x',x_SX,'f',objective,'g',constraint);
    S   = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    g_vvl(i) = -full(sol.f);
end