function gg_out = gg(mpc,v,v0,theta,Phi0,p_inj,q_inj)
    %% Index setting
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
    % added idx
    ptc_factor = 26;
    %% 
    Nbranch = size(mpc.branch,1);
    Nbus = size(mpc.bus,1);
    Ngen = size(mpc.gen,1);
    idx_fr = mpc.branch(:,F_BUS);
    idx_to = mpc.branch(:,T_BUS);
    E_fr = sparse(idx_fr, 1:Nbranch, 1, Nbus, Nbranch);
    E_to = sparse(idx_to, 1:Nbranch, 1, Nbus, Nbranch);
    E = E_fr - E_to;
    Phi = E'*theta;
    id_gen = mpc.gen(:,GEN_BUS);
    Cg = sparse(id_gen,1:Ngen,1,Nbus,Ngen);
    onoff_pq = zeros(Nbus,1);
    idx_pq = mpc.bus(:,BUS_TYPE)==1;
    onoff_pq(idx_pq) = 1;
    
    alpha0 = mpc.gen(:,ptc_factor);
    delta0 = mpc.delta;

    g_vvcos = v(idx_fr).*v(idx_to).*cos(Phi) - onoff_pq(idx_fr).*v(idx_fr).*v0(idx_to).*cos(Phi0)...
        -v0(idx_fr).*onoff_pq(idx_to).*v(idx_to).*cos(Phi0) + v0(idx_fr).*v0(idx_to).*sin(Phi0).*Phi;
    g_vvsin = v(idx_fr).*v(idx_to).*sin(Phi) - onoff_pq(idx_fr).*v(idx_fr).*v0(idx_to).*sin(Phi0)...
        -v0(idx_fr).*onoff_pq(idx_to).*v(idx_to).*sin(Phi0) - v0(idx_fr).*v0(idx_to).*cos(Phi0).*Phi;
    g_vv = v.^2 - 2*onoff_pq.*v0.*v;
    gg_out = [p_inj-Cg*alpha0*delta0; q_inj; g_vvcos; g_vvsin; g_vv];
end

