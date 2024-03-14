%%
clc
clear
close all
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
%% standard OPF as reference result
% mpc        = ext2int(loadcase('case15_modified'));
mpc        = ext2int(loadcase('case33bw'));
% mpc        = ext2int(loadcase('case69'));
% mpc        = ext2int(loadcase('case141'));
idx_remove_branch = ~mpc.branch(:,BR_STATUS);

mpc.branch(idx_remove_branch,:)=[];


%% reference solution by matpower
% remove all shunt element
mpc.branch(:,BR_B) = 0;
% remove all line limit
mpc.branch(:,RATE_A) = 0;
mpc.branch(:,RATE_B) = 0;
mpc.branch(:,RATE_C) = 0;
% reference solution
mpc = runopf(mpc);
idx_remove_branch = ~mpc.branch(:,BR_STATUS);

mpc.branch(idx_remove_branch,:)=[];

%% SOCP problem formulation
% Extract necessary data from Matpower casefile
Nbus       = size(mpc.bus,1);
Nbranch    = size(mpc.branch,1);
Ngen       = size(mpc.gen,1);
% number of all decision variables
Nx         = Nbus + 3*Nbranch + 2*Ngen;
baseMVA    = mpc.baseMVA;              % ratio
% lower & upper bounds for voltiage magnitude (square)
U0         = mpc.bus(:,VM).^2;
Umin       = mpc.bus(:,VMIN).^2;
Umax       = mpc.bus(:,VMAX).^2;
L0         = ones(Nbranch,1);
Lmax       = 1e3 *  ones(Nbranch,1);
Lmin       = sparse(Nbranch,1);
% reference bus
ref_idx    = find(mpc.bus(:,BUS_TYPE)==REF);
% lower & upper bounds for generators
gen_idx    = mpc.gen(:,GEN_BUS);             % gen bus
gen_cost   = mpc.gencost(:,5:end);     % objective coefficients
Pg0        = mpc.gen(:,PG);
Qg0        = mpc.gen(:,QG);
Pgmin      = mpc.gen(:,PMIN)/baseMVA;
Qgmin      = mpc.gen(:,QMIN)/baseMVA;
Pgmax      = mpc.gen(:,PMAX)/baseMVA;
Qgmax      = mpc.gen(:,QMAX)/baseMVA;
Cg         = sparse(1:Ngen,gen_idx,ones(Ngen,1),Ngen,Nbus);
% demand (load)
Pd         = mpc.bus(:,PD)/baseMVA;   
Qd         = mpc.bus(:,QD)/baseMVA;   
% idx of from/to bus of a branch
from_bus   = mpc.branch(:, F_BUS);                           
to_bus     = mpc.branch(:, T_BUS);
Cf         = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
Ct         = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
% branch info
branch_r   = mpc.branch(:,BR_R);
branch_x   = mpc.branch(:,BR_X);
R   = diag(mpc.branch(:,BR_R));
X   = diag(mpc.branch(:,BR_X));
Z2  = R.^2 + X.^2; 
% parameters for relation between U, P, Q
B              = (Cf + Ct)';
A              = [zeros(Nbranch,1) eye(Nbranch)]*B - eye(Nbranch);
D              = (eye(Nbranch) - A)\eye(size(A));
DR             = (eye(Nbranch) - A)\(A*R);
DX             = (eye(Nbranch) - A)\(A*X);
DX_p = max(DX,zeros(size(DX))).*(DX>0);
DX_m = min(DX,zeros(size(DX))).*(DX<0);
Mp = D'*R*D;
Mq = D'*X*D;

H   = D'*(2*(R*DR+X*DX) + Z2);
% line limit
Qijmax = 100*ones(Nbranch,1);
Qij0   = 1*ones(Nbranch,1);
Pijmax = 100*ones(Nbranch,1);
Pij0   = 1*ones(Nbranch,1);
% conic constraint 
B = cell(Nbranch,1);
d = cell(Nbranch,1);
for i = 1:Nbranch
    bb = sparse(3,Nx);
    bb(1,from_bus(i)) = -1; % Uij
    bb(1,Nbus+i) = 1;       % Lij
    bb(2,(Nbus+Nbranch)+i) = 2; % 2Pij
    bb(3,(Nbus+Nbranch*2)+i) = 2; %2Qij
    B{i} = bb;
    dd = sparse(Nx,1);
    dd(from_bus(i)) = 1;
    dd(Nbus+i)=1;
    d{i} = dd;
end
% dependent generator - PCC
idx_dg = 1;
% independent generator
if Ngen>1
    idx_ig = 2:Ngen;
else
    idx_ig = [];
end

% A(1:Nbranch, 1:(Nbus+Nbranch)) = horzcat(-Cf,speye(Nbranch));
% A((Nbranch+1):(Nbranch*2),(Nbus+Nbranch+1):(Nbus+Nbranch*2))=2*speye(Nbranch));
% A((Nbranch*2+1):(Nbranch*3),(Nbus+Nbranch+1):(Nbus+Nbranch*2))=2*speye(Nbranch));

%% referenece solutions

u = (mpc.bus(:,8)).^2;
pg = mpc.gen(:,PG)/mpc.baseMVA;
qg = mpc.gen(:,QG)/mpc.baseMVA;
pij = mpc.branch(:,PF)/mpc.baseMVA;
qij = mpc.branch(:,QF)/mpc.baseMVA;
lij = (pij.^2+qij.^2)./(Cf*u);
pij_to = mpc.branch(:,PT)/mpc.baseMVA;
x_ref_DisFlow   = vertcat(u,pij,qij,pg(idx_dg),qg(idx_dg),lij,pg(idx_ig),qg(idx_ig));


% x_ref_LinDisFlow   = vertcat(u,pij,qij,pg,qg);

% error test
% current_error = norm(pij_to + (pij - branch_r.*lij),inf)
% power_balance = norm(Pd - Cg'*pg + Cf'*pij - Ct'*(pij - branch_r.*lij),inf)
% voltage_error = norm(u(from_bus)-u(to_bus) -2*(branch_r .* pij + branch_x .* qij) ...
%     + (branch_r.^2+branch_x.^2) .* lij,inf)
% conic_error = norm(pij.^2 + qij .^2 - lij .* Cf*u,inf)
% % conic_error = zeros(Nbranch,1);
% % for i = 1:Nbranch
% %     conic_error(i) = norm(B{i}*x_ref_DisFlow,2) - d{i}'*x_ref_DisFlow;
% % end
% % conic_error = norm(conic_error,2)
%% distflow model
% decision variables z = (x,u)
import casadi.*
U          = SX.sym('U',Nbus,1);
Lij        = SX.sym('L',Nbranch,1);
Pij        = SX.sym('Pij',Nbranch,1);
Qij        = SX.sym('Qij',Nbranch,1);
Pg         = SX.sym('Pg',Ngen,1);
Qg         = SX.sym('Qg',Ngen,1);
x          = vertcat(U,Pij,Qij,Pg(idx_dg),Qg(idx_dg),Lij);
u          = vertcat(Pg(idx_ig),Qg(idx_ig));
z       = [x;u]; 
% if Ngen>1
%     u      = vertcat(Pg(2:end),Qg(2:end));
%     eta    = [x;u];
% else
%     
% end
% equality constraints
volt_eq    = U(from_bus)-U(to_bus) -2*(branch_r .* Pij + branch_x .* Qij) ...
    + (branch_r.^2+branch_x.^2) .* Lij;
% volt_eq = U(2:end) - (2*Mq*(Cg(2:end,2:end)*Qg(2:end)-Qd(2:end))+2*Mp*(Cg(2:end,2:end)*Pg(2:end)-Pd(2:end)))-1 + H*Lij;

% power balance
pf_p_eq    = Pd - Cg'*Pg + Cf'*Pij - Ct'*(Pij - branch_r.*Lij);
pf_q_eq    = Qd - Cg'*Qg + Cf'*Qij - Ct'*(Qij - branch_x.*Lij);
% conic constraint
% conic_ineq = Pij.^2 + Qij .^2 - Lij .* Cf*U;
conic_ineq = Cf*U - (Pij.^2 + Qij .^2)./Lij ;
% ref bus
ref_eq     = U(ref_idx) - mpc.bus(ref_idx, VM);
% objective
f          = baseMVA^2*Pg'*diag(gen_cost(:,1))*Pg...
              + baseMVA*Pg'*gen_cost(:,2);
% interface with casadi (IPOPT)
z0         = vertcat(U0,Pij0,Qij0,Pg0(idx_dg),Qg0(idx_dg),L0,Pg0(idx_ig),Qg0(idx_ig));
Umin(1) = 0.8;
Umax(1) = 1.2;
lbz        = vertcat(Umin,-Pijmax,-Qijmax,Pgmin(idx_dg),Qgmin(idx_dg),Lmin,Pgmin(idx_ig),Qgmin(idx_ig));
ubz        = vertcat(Umax,Pijmax,Qijmax,Pgmax(idx_dg),Qgmax(idx_dg),Lmax,Pgmax(idx_ig),Qgmax(idx_ig));
ffun       = f;
gfun       = vertcat(ref_eq, volt_eq, pf_p_eq, pf_q_eq, conic_ineq);
% lbg        = vertcat(zeros(2*Nbus+Nbranch+1,1),-inf*ones(Nbranch,1));  % SOCP
lbg        = vertcat(zeros(2*Nbus+2*Nbranch+1,1));  % DistFlow
ubg        = zeros(2*Nbus+2*Nbranch+1,1);
% solve SOC OPF problem
nlp = struct('x',z,'f',ffun,'g',gfun);
options.ipopt.tol         = 1.0e-8;
options.ipopt.print_level = 5;
options.print_time        = 1;
options.ipopt.max_iter    = 100;
S = nlpsol('solver','ipopt', nlp,options);
sol = S('x0', z0,'lbg', lbg,'ubg', ubg,...
        'lbx', lbz, 'ubx', ubz);
xopt = full(sol.x);
fval = full(sol.f);
error_distflow = norm(xopt -x_ref_DisFlow)
P    = xopt((1:Nbranch)+Nbus);
Q    = xopt((1:Nbranch)+Nbus+Nbranch);
L    = xopt((1:Nbranch)+Nbus+Nbranch*2+2);

%% reference derivatives - casadi
jac_casadi = jacobian(gfun,x);
derivative = Function('sens',{x},{jac_casadi}); 
idx_x      = 1:Nbus*1 + Nbranch*3+2;
jac_ref = full( derivative(xopt(idx_x)));
Nx    = Nbus*1 + Nbranch*2+2;
D_ref = jac_ref(1:Nx,:);
E_ref = D_ref(:,Nx+1:end);
D_ref = D_ref(:,1:Nx);
F_ref = jac_ref(Nx+1:end,:);
G_ref = F_ref(:,Nx+1:end);
F_ref = F_ref(:,1:Nx);

%% 1st derivative
Cft  = Cf-Ct;
e1   = sparse(ref_idx,1,1,Nbus,1);
D    = full([1, zeros(1,Nbus*2+Nbranch);
        Cft, -2*diag(branch_r),-2*diag(branch_x), zeros(Nbranch,2);
        zeros(Nbus,Nbus),Cft',zeros(Nbus,Nbranch),-e1,zeros(Nbus,1);
        zeros(Nbus,Nbus+Nbranch),Cft',zeros(Nbus,1),-e1]);
RR  = sparse(diag(branch_r));
XX  = sparse(diag(branch_x));
PP  = sparse(diag(P));
QQ  = sparse(diag(Q));
Ga  = sparse(diag(L.^(-1)));

E   = [sparse(1,Nbranch);RR^2+XX^2;Ct' *RR;Ct'*XX];
F   = [Cf, - 2*PP*Ga, -2*QQ*Ga,sparse(Nbranch,2)];
G   = (PP^2+QQ^2) * Ga^2;
D_error = norm(D - D_ref)
E_error = norm(E_ref - E)
F_error = norm(F- F_ref)
G_error = norm(G_ref - G)








