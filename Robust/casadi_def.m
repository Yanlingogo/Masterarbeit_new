import casadi.*
%% set of all variables
v_u = MX.sym('v_u', Nbus, 1);
v_l = MX.sym('v_l', Nbus, 1);
Delta_v_u = MX.sym('Delta_v_u', Nbus, 1);
Delta_v_l = MX.sym('Delta_v_l', Nbus, 1);
Phi_u = MX.sym('Phi_u', Nbranch, 1);
Phi_l = MX.sym('Phi_l', Nbranch, 1);
Delta_Phi_u = MX.sym('Delta_Phi_u', Nbranch, 1);
Delta_Phi_l = MX.sym('Delta_Phi_l', Nbranch, 1);

delta_u = MX.sym('delta_u');
delta_l = MX.sym('delta_l');
Ddelta_u = MX.sym('Ddelta_u');
Ddelta_l = MX.sym('Ddelta_l');

g_u_pinj = MX.sym('g_u_pinj', Nbus, 1);
g_l_pinj = MX.sym('g_l_pinj', Nbus, 1);
g_u_qinj = MX.sym('g_u_qinj', Nbus, 1);
g_l_qinj = MX.sym('g_l_qinj', Nbus, 1);
g_u_vvsin = MX.sym('g_u_vvsin', Nbranch, 1);
g_l_vvsin = MX.sym('g_l_vvsin', Nbranch, 1);
g_u_vvcos = MX.sym('g_u_vvcos', Nbranch, 1);
g_l_vvcos = MX.sym('g_l_vvcos', Nbranch, 1);
g_u_vv = MX.sym('g_u_vv', Nbus, 1);
g_l_vv = MX.sym('g_l_vv', Nbus, 1);

Psi_u_vvcos = MX.sym('Psi_u_vvcos', Nbranch, 1);
Psi_l_vvcos = MX.sym('Psi_l_vvcos', Nbranch, 1);
Psi_u_vvsin = MX.sym('Psi_u_vvsin', Nbranch, 1);
Psi_l_vvsin = MX.sym('Psi_l_vvsin', Nbranch, 1);
Psi_u_vv = MX.sym('Psi_u_vv', Nbus, 1);
Psi_l_vv = MX.sym('Psi_l_vv', Nbus, 1);

Delta_vv_u_uu = MX.sym('Delta_vv_u_uu', Nbranch, 1);
Delta_vv_u_ll = MX.sym('Delta_vv_u_ll', Nbranch, 1);
Delta_vv_u_ul = MX.sym('Delta_vv_u_ul', Nbranch, 1);
Delta_vv_u_lu = MX.sym('Delta_vv_u_lu', Nbranch, 1);
Delta_vv_l_uu = MX.sym('Delta_vv_l_uu', Nbranch, 1);
Delta_vv_l_ll = MX.sym('Delta_vv_l_ll', Nbranch, 1);
Delta_vv_l_ul = MX.sym('Delta_vv_l_ul', Nbranch, 1);
Delta_vv_l_lu = MX.sym('Delta_vv_l_lu', Nbranch, 1);

Delta_cos_u_u = MX.sym('Delta_cos_u_u', Nbranch, 1);
Delta_cos_u_l = MX.sym('Delta_cos_u_l', Nbranch, 1);
Delta_cos_l_u = MX.sym('Delta_cos_l_u', Nbranch, 1);
Delta_cos_l_l = MX.sym('Delta_cos_l_l', Nbranch, 1);
Delta_sin_u_u = MX.sym('Delta_sin_u_u', Nbranch, 1);
Delta_sin_u_l = MX.sym('Delta_sin_u_l', Nbranch, 1);
Delta_sin_l_u = MX.sym('Delta_sin_l_u', Nbranch, 1);
Delta_sin_l_l = MX.sym('Delta_sin_l_l', Nbranch, 1);

p_line_fr_u = MX.sym('p_line_fr_u', Nbranch, 1);
q_line_fr_u = MX.sym('q_line_fr_u', Nbranch, 1);
p_line_to_u = MX.sym('p_line_to_u', Nbranch, 1);
q_line_to_u = MX.sym('q_line_to_u', Nbranch, 1);

pg_opt = MX.sym('pg_opt', Ngen, 1);
pl_opt = MX.sym('pl_opt', Npq, 1);
ql_opt = MX.sym('ql_opt', Npq, 1);
pg_opt_u = MX.sym('pg_opt_u', Ngen, 1);
alpha_opt = MX.sym('alpha_opt', Ngen, 1);
Dalpha_opt = MX.sym('Dalpha_opt', Ngen, 1);

cost_opt = MX.sym('cost_opt');
gamma_opt = MX.sym('gamma_opt');
slack_pg = MX.sym('slack_pg');
slack_qg = MX.sym('slack_qg');
slack_line = MX.sym('slack_line');
margin_opt = MX.sym('margin_opt');
%% Describe the problem
% Define initial value
x0 = vertcat(v0,v0,zeros(2*Nbus,1),Phi0,Phi0,zeros(2*Nbranch,1),delta0,delta0,zeros(2+8*Nbus+28*Nbranch,1),...
    pg0,ppq0,qpq0,pg0,alpha0,zeros(Ngen,1),mpc.cost,zeros(5,1));
% Define the lower bound of variables
lbx = vertcat(vmin,vmin,vmin-v0,vmin-v0,Phi_min,Phi_min,Phi_min-Phi0,Phi_min-Phi0,...
    -inf(10+8*Nbus+28*Nbranch+4*Ngen+2*Npq,1));
% Define the upper bound of variables
ubx = vertcat(vmax,vmax,vmax-v0,vmax-v0,Phi_max,Phi_max,Phi_max-Phi0,Phi_max-Phi0,...
    inf(10+8*Nbus+28*Nbranch+4*Ngen+2*Npq,1));
% Define gfun
% Equality constraints
g_v1 = v_u(id_gen)-v_l(id_gen);
g_v2 = v_u-Delta_v_u;
g_v3 = v_l-Delta_v_l;
g_Phi1 = Phi_u-Delta_Phi_u;
g_Phi2 = Phi_l-Delta_Phi_l;
g_delta1 = delta_u-Ddelta_u;
g_delta2 = delta_l-Ddelta_l;
g_alpha = alpha_opt-Dalpha_opt;
g_alpha2 = alpha_opt; % ==alpha0
% inequality constraints 
g_vv1 = Delta_vv_u_uu-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)-1/4*(Delta_v_u(idx_fr)+Delta_v_u(idx_to)).^2;
g_vv2 = Delta_vv_u_ll-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)-1/4*(Delta_v_l(idx_fr)+Delta_v_l(idx_to)).^2;
g_vv3 = Delta_vv_u_ul-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)-1/4*(Delta_v_u(idx_fr)+Delta_v_l(idx_to)).^2;
g_vv4 = Delta_vv_u_lu-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)-1/4*(Delta_v_l(idx_fr)+Delta_v_u(idx_to)).^2;
% >= / <= 
g_vv5 = Delta_vv_l_uu-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)+1/4*(Delta_v_u(idx_fr)+Delta_v_u(idx_to)).^2;
g_vv6 = Delta_vv_l_ll-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)+1/4*(Delta_v_l(idx_fr)+Delta_v_l(idx_to)).^2;
g_vv7 = Delta_vv_l_ul-Delta_v_u(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_l(idx_to)+1/4*(Delta_v_u(idx_fr)+Delta_v_l(idx_to)).^2;
g_vv8 = Delta_vv_l_lu-Delta_v_l(idx_fr).*v0(idx_to)-v0(idx_fr).*Delta_v_u(idx_to)+1/4*(Delta_v_l(idx_fr)+Delta_v_u(idx_to)).^2;
% >= / <=
g_cos1 = Delta_cos_u_u+sin(Phi0).*Delta_Phi_u-1/2*(Delta_Phi_u).^2;
g_cos2 = Delta_cos_u_l+sin(Phi0).*Delta_Phi_l-1/2*(Delta_Phi_l).^2;
g_cos3 = Delta_cos_l_u+sin(Phi0).*Delta_Phi_u+1/2*(Delta_Phi_u).^2;
g_cos4 = Delta_cos_l_l+sin(Phi0).*Delta_Phi_l+1/2*(Delta_Phi_l).^2;
% >= / <=
g_sin1 = Delta_sin_u_u-cos(Phi0).*Delta_Phi_u-1/2*(Delta_Phi_u).^2;
g_sin2 = Delta_sin_u_l-cos(Phi0).*Delta_Phi_l-1/2*(Delta_Phi_l).^2;
g_sin3 = Delta_sin_l_u-cos(Phi0).*Delta_Phi_u+1/2*(Delta_Phi_u).^2;
g_sin4 = Delta_sin_l_l-cos(Phi0).*Delta_Phi_l+1/2*(Delta_Phi_l).^2;
% >= ψ_vvcos0
g_gcos11 = g_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_uu+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
g_gcos12 = g_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_uu+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
g_gcos13 = g_u_vvcos - cos(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_ul+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
g_gcos14 = g_u_vvcos - cos(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_ul+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
g_gcos15 = g_u_vvcos - cos(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_lu+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
g_gcos16 = g_u_vvcos - cos(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_lu+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
g_gcos17 = g_u_vvcos - cos(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_ll+Delta_cos_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
g_gcos18 = g_u_vvcos - cos(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_ll+Delta_cos_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
% <= ψ_vvcos0
g_gcos21 = g_l_vvcos - cos(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_uu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
g_gcos22 = g_l_vvcos - cos(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_uu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
g_gcos23 = g_l_vvcos - cos(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_uu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
g_gcos24 = g_l_vvcos - cos(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_uu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
g_gcos25 = g_l_vvcos - cos(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_uu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
g_gcos26 = g_l_vvcos - cos(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_uu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
g_gcos27 = g_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_uu-Delta_cos_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_u;
g_gcos28 = g_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_uu-Delta_cos_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*cos(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*cos(Phi0) - v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi_l;
% >= Psi_vvsin0
g_gsin11 = g_u_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_uu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin12 = g_u_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_uu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin13 = g_u_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_ul+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin14 = g_u_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_ul+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin15 = g_u_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_lu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin16 = g_u_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_lu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin17 = g_u_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_ll+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin18 = g_u_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_u_ll+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
% >= Psi_vvsin0
g_gsin21 = g_u_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_uu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin22 = g_u_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_uu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin23 = g_u_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_ul+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin24 = g_u_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_ul+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin25 = g_u_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_lu+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin26 = g_u_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_lu+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin27 = g_u_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_ll+Delta_sin_u_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin28 = g_u_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_u_l - 1/4*(Delta_vv_l_ll+Delta_sin_u_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
% <= Psi_vsin0
g_gsin31 = g_l_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_uu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin32 = g_l_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_uu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin33 = g_l_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_ul-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin34 = g_l_vvsin - sin(Phi0).*Delta_vv_l_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_ul-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin35 = g_l_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_lu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin36 = g_l_vvsin - sin(Phi0).*Delta_vv_l_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_lu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin37 = g_l_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_l_ll-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin38 = g_l_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_ll-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
% <= Psi_vsin0
g_gsin41 = g_l_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_uu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin42 = g_l_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_uu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin43 = g_l_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_ul-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin44 = g_l_vvsin - sin(Phi0).*Delta_vv_u_ul - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_ul-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_u(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin45 = g_l_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_lu-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin46 = g_l_vvsin - sin(Phi0).*Delta_vv_u_lu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_lu-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_u(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
g_gsin47 = g_l_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_u + 1/4*(Delta_vv_u_ll-Delta_sin_l_u).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_u;
g_gsin48 = g_l_vvsin - sin(Phi0).*Delta_vv_u_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_ll-Delta_sin_l_l).^2 + onoff_pq(idx_fr).*v_l(idx_fr).*v0(idx_to).*sin(Phi0) + v0(idx_fr).*onoff_pq(idx_to).*v_l(idx_to).*sin(Phi0) + v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi_l;
% >= 0
g_gvv1 = g_u_vv - v_u.^2 + 2*v0.*onoff_pq.*v_u;
g_gvv2 = g_u_vv - v_l.^2 + 2*v0.*onoff_pq.*v_l;
% <= Psi_vv0
g_gvv3 = g_l_vv - 2*v0.*(v_u-v0) + 2*v0.*onoff_pq.*v_u;
g_gvv4 = g_l_vv - 2*v0.*(v_l-v0) + 2*v0.*onoff_pq.*v_l;
% >= Psi_vvcos0
g_Psi1 = Psi_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_u - 1/4*(Delta_vv_u_uu+Delta_cos_u_u).^2;
g_Psi2 = Psi_u_vvcos - cos(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_cos_u_l - 1/4*(Delta_vv_u_uu+Delta_cos_u_l).^2;
g_Psi3 = Psi_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_u + 1/4*(Delta_vv_l_ll-Delta_cos_l_u).^2;
g_Psi4 = Psi_l_vvcos - cos(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_cos_l_l + 1/4*(Delta_vv_l_ll-Delta_cos_l_l).^2;
% >= Psi_vvsin
g_Psi5 = Psi_u_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_u_uu+Delta_sin_u_u).^2;
g_Psi6 = Psi_u_vvsin - sin(Phi0).*Delta_vv_l_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_u_u - 1/4*(Delta_vv_l_uu+Delta_sin_u_u).^2;
g_Psi7 = Psi_u_vvsin - sin(Phi0).*Delta_vv_l_ll - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_l_ll-Delta_sin_l_l).^2;
g_Psi8 = Psi_u_vvsin - sin(Phi0).*Delta_vv_u_uu - v0(idx_fr).*v0(idx_to).*Delta_sin_l_l + 1/4*(Delta_vv_u_uu-Delta_sin_l_l).^2;

g_Psi9 = Psi_u_vv - v_u.^2; %>=0
g_Psi10 = Psi_u_vv - v_l.^2;%>=0
% <= Psi_vv0
g_Psi11 = Psi_l_vv - 2*v0.*Delta_v_u;
g_Psi12 = Psi_l_vv - 2*v0.*Delta_v_l;

g_g1 = g_u_pinj - Cg*pg_opt + Cl*pl_opt; %>=0
g_g2 = g_l_pinj - Cg*pg_opt + Cl*pl_opt; %<=0
g_pgmax = pg_opt + alpha0.*delta_u+slack_pg; %<=pg_max
g_pgmin = pg_opt + alpha0.*delta_l-slack_pq; %>=pg_min

zeta = sqrt(sum((Cl(idx_pvs,:)*diag(power_factor)*chol(Sigma_0,'upper')).^2,2));

g_M1 = M_ineq_minus*[zeros(Nbus,1); Cg*qg_max-Cl*ql_opt;Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]+M_ineq_plus*[zeros(Nbus,1);Cg*qg_max-Cl*ql0;Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv]-slack_qg-zeta.*gamma_opt;%>=0
g_M2 = M_ineq_plus*[zeros(Nbus,1); Cg*qg_min-Cl*ql_opt;Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]+M_ineq_minus*[zeros(Nbus,1);Cg*qg_min-Cl*ql0;Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv]+slack_qg+zeta.*gamma_opt;%<=0
g_M3 = M_line_plus*[zeros(2*Nbus,1);Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]+M_line_minus*[zeros(2*Nbus,1);Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv] - [p_line_fr_u;p_line_to_u;q_line_fr_u;q_line_to_u];%<=0
g_M4 = -M_line_minus*[zeros(2*Nbus,1);Psi_u_vvcos;Psi_u_vvsin;Psi_u_vv]-M_line_plus*[zeros(2*Nbus,1);Psi_l_vvcos;Psi_l_vvsin;Psi_l_vv] - [p_line_fr_u;p_line_to_u;q_line_fr_u;q_line_to_u];%<=0
g_line1 = p_line_fr_u.^2+q_line_fr_u.^2+slack_line.^2;% <= s_line_max.^2
g_line2 = p_line_to_u.^2+q_line_to_u.^2+slack_line.^2;% <= s_line_max.^2

xi = sqrt(sum((K(:,1:2*Nbus)*[Cl;Cl*diag(power_factor)]*chol(Sigma_0,'upper')).^2,2));
g_K1 = K_plus*[g_u_pinj; -Cl*ql_opt; g_u_vvcos; g_u_vvsin;g_u_vv]+K_minus*[g_l_pinj; -Cl*ql_opt; g_l_vvcos; g_l_vvsin;g_l_vv]+xi.*gamma_opt-[Phi_u;v_u(idx_pq);delta_u;-Phi_l;-v_l(idx_pq);-delta_l];%<=0
g_margin1 = slack_pg - margin_opt;
g_margin2 = slack_qg - margin_opt;
g_margin3 = slack_line -margin_opt;
g_margin4 = gamma_opt - margin_opt;
g_pgopt = pg_opt+alpha0.*delta_u-pg_opt_u; %<=0
g_cost = gencost(:,COST)'*(pg_opt_u).^2+gencost(:,COST+1)'*pg_opt_u+sum(gencost(:,COST+2));%<=cost_opt
gfun = vertcat()
% Define the lower bound of constraints

% Define the upper bound of constraints

% Define objective function
options.ipopt.tol         = 1.0e-8;
options.ipopt.print_level = 5;
options.print_time        = 5;
options.ipopt.max_iter    = 100;
if option == "margin" && isempty(target_mpc)
    % Constraints
    ffun = gamma_opt;
    nlp = struct('x',x,'f',ffun,'g',gfun);
    S = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    if sol.stats == 0
        %gamma0 = double(gamma_opt);
        %mpc.uncertainty.Sigma0 = Sigma0;
        %mpc.uncertainty.gamma0 = gamma0;
    else
        disp('Error in solving the problem.');
    end
elseif option == "margin" && ~isempty(target_mpc)

elseif option == "obj"
    Constraints = [Constraints, pl_opt == ppq0, ql_opt == qpq0, margin_opt>=0, gamma_opt>=gamma0];
    ffun = cost_opt;
    nlp = struct('x',x,'f',ffun,'g',gfun);
    S = nlpsol('solver','ipopt', nlp,options);
    sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
            'lbx', lbx, 'ubx', ubx);
    if sol.stats == 0
        mpc.gen(:,PG) = value(pg_opt);
        mpc.gen(:,end) = value(alpha_opt(id_gen));
        mpc.gen(:,VG) = value(v_u(id_gen));
    else
        disp('Error in solving the problem.');
    end
end

% Define g