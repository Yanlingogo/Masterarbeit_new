clc;
clear;
%%Index setting
% bus idx
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
% branch idx
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% gen idx
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
% cost idx
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

mpc = loadcase('case18');
mpc = runpf(mpc);
mpc = ext2int(mpc);
% Branch flow: 
y_line = mpc.branch(:,11)./(mpc.branch(:,3)+1i*mpc.branch(:,4));
G_line = diag(real(y_line));
B_line = diag(imag(y_line));
idx_fr = mpc.branch(:,1);
idx_to = mpc.branch(:,2);
Nbus        = size(mpc.bus,1);        %size(,1) get the rows of matrix
Nbranch     = size(mpc.branch,1);
E_fr = sparse(idx_fr, 1:Nbranch, 1, Nbus, Nbranch);
E_to = sparse(idx_to, 1:Nbranch, 1, Nbus, Nbranch);
E = E_fr - E_to;
v = mpc.bus(:,8);
ang = deg2rad(mpc.bus(:,9));
Phi = E'*ang;
p_kl = v(idx_fr).*v(idx_to).*(G_line*cos(Phi)+B_line*sin(Phi));
q_kl = v(idx_fr).*v(idx_to).*(G_line*sin(Phi)-B_line*cos(Phi));
p_line = mpc.branch(:,14)-mpc.branch(:,16);

[Ybus, Yf, Yt] = makeYbus(mpc);
Gbus           = real(Ybus);
Bbus           = imag(Ybus);
Gf             = real(Yf);
Bf             = imag(Yf);
Gt             = real(Yt);
Bt             = imag(Yt);
U = v .* cos(ang);
W = v .* sin(ang);

% from side
pij = W(idx_fr) .* (Gf*W + Bf*U) + U(idx_fr) .* (Gf*U - Bf*W);
qij = W(idx_fr) .* (Gf*U - Bf*W) - U(idx_fr) .* (Gf*W + Bf*U);
% to side
pji = W(idx_to) .* (Gt*W + Bt*U) + U(idx_to) .* (Gt*U - Bt*W);
qji = W(idx_to) .* (Gt*U - Bt*W) - U(idx_to) .* (Gt*W + Bt*U);
p_line2 = (pij-pji)*mpc.baseMVA;


id_slack = find(mpc.bus(:,BUS_TYPE)==3);
connected_buses = unique([mpc.branch(mpc.branch(:, 1) == id_slack, 2);...
    mpc.branch(mpc.branch(:, 2) == id_slack, 1)]);
id_cline = find(mpc.branch(:, 1) == id_slack | mpc.branch(:, 2) == id_slack); % connecting line
pij_slack = W(id_slack) .* (Gf(id_cline,:)*W + Bf(id_cline,:)*U) +...
    U(id_slack) .* (Gf(id_cline,:)*U - Bf(id_cline,:)*W);
% to connected bus
pji_slack = W(connected_buses) .* (Gt(id_cline,:)*W + Bt(id_cline,:)*U) +...
    U(connected_buses) .* (Gt(id_cline,:)*U - Bt(id_cline,:)*W);
% detect active line limit - 0 as no line limit
p = pij_slack - pji_slack;