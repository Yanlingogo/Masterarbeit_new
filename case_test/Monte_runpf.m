function [mpc] = Monte_runpf(mpc)
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

    Nbus = size(mpc.bus,1);
    Ngen = size(mpc.gen,1);
    idx_load = find((mpc.bus(:,PD).^2+mpc.bus(:,QD).^2) ~= 0);
    idx_pq   = find(mpc.bus(:,BUS_TYPE)==1); 
    Npq   = sum(mpc.bus(:,BUS_TYPE)==1);
    Nload = numel(idx_load);
    bus = mpc.bus;
    gen = mpc.gen;

    bus_type = bus(:,BUS_TYPE);
    slack_bus = bus_type == 3;
    idx_nslack = find(bus_type ~= 3);
    id_gen = gen(:,GEN_BUS);

    gencost = mpc.gencost(:,COST:end);

    Cg = sparse(id_gen, (1:Ngen), ones(1, Ngen), Nbus, Ngen);
    Cl = sparse(idx_load, (1:Nload), ones(1, Nload), Nbus, Nload);

    pg0 = gen(:,GEN_STATUS).*gen(:,PG)/mpc.baseMVA; 
    pl0 = bus(idx_load,PD)/mpc.baseMVA;
    ql0 = bus(idx_load,QD)/mpc.baseMVA;
    pinj0 = Cg*pg0 -Cl*pl0; qinj0 = -Cl*ql0;
    Sinj = (pinj0 + 1i * qinj0); % in p.u.
    
    [Y,Yf,Yt] = makeYbus(mpc);

    vm = bus(:,VM);
    vm(id_gen) = gen(:,VG);
    va = bus(:,VA); % Angle value in radian system
    va(slack_bus) = 0;

    max_iter = 30;
    nf = zeros(max_iter,1); ndx = zeros(max_iter,1);

    for iter = 1:max_iter
        v_cpx = vm.*cos(va)+1i*(vm.*sin(va));
        e_cpx = v_cpx./(vm);
        S_iter = diag(v_cpx)*conj(Y*v_cpx);
        S_bal = Sinj - S_iter;
        f=[real(S_bal(idx_nslack));imag(S_bal(idx_pq))];
        % J1 = real(diag(v_cpx) * conj(Y .* (ones(Nbus, 1) * transpose(1i *v_cpx))) + diag(1i * v_cpx) * diag(conj(Y * v_cpx)));
        % J2 = real(diag(v_cpx) * conj(Y .* (ones(Nbus, 1) * transpose(v_cpx./vm))) + diag(v_cpx ./ vm) * diag(conj(Y * v_cpx)));
        % J3 = imag(diag(v_cpx) * conj(Y .* (ones(Nbus, 1) * transpose(1i *v_cpx))) + diag(1i * v_cpx) * diag(conj(Y * v_cpx)));
        % J4 = imag(diag(v_cpx) * conj(Y .* (ones(Nbus, 1) * transpose(v_cpx./vm))) + diag(v_cpx ./ vm) * diag(conj(Y * v_cpx)));
        
        Z_1 = diag(v_cpx)*conj(Y)*diag(conj(v_cpx));
        Z_2 = diag(v_cpx)*conj(Y)*diag(conj(e_cpx));
        J1 = imag(Z_1)  - diag(imag(S_iter));
        J2 = real(Z_2)  + diag(real(S_iter)./vm);
        J3 = -real(Z_1) + diag(real(S_iter));
        J4 = imag(Z_2)  + diag(imag(S_iter)./vm);
        
        J  = [J1(idx_nslack,idx_nslack) J2(idx_nslack,idx_pq); ...
              J3(idx_pq,idx_nslack) J4(idx_pq,idx_pq)];
        dx = J\f;

        va(idx_nslack) = va(idx_nslack)+dx(1:Nbus-1);
        vm(idx_pq) = vm(idx_pq)+dx(Nbus:end);

        nf(iter) = norm(f); 
        ndx(iter) = norm(dx);
        if isnan(nf(iter)) || ((nf(iter)<1e-8)&&(ndx(iter)<1e-8))
            mpc.success = 1;
            break; 
        end 
    end

    mpc.bus(:,VM) = vm;
    mpc.bus(:,VA) = rad2deg(va); % va为弧度值


    v_cpx = vm.*cos(va)+1i*vm.*sin(va);
    qg_inj = imag(v_cpx.*conj(Y*v_cpx))+Cl*ql0;
    pg_inj = real(v_cpx.*conj(Y*v_cpx))+Cl*pl0;
    mpc.gen(:,PG) = Cg'*pg_inj*mpc.baseMVA;
    mpc.gen(:,QG) = Cg'*qg_inj*mpc.baseMVA;
   
    
    from_bus = mpc.branch(:, F_BUS); 
    to_bus   = mpc.branch(:, T_BUS);
    Sf = v_cpx(from_bus) .* conj(Yf * v_cpx) * mpc.baseMVA;
    St = v_cpx(to_bus) .* conj(Yt * v_cpx) * mpc.baseMVA;
    mpc.branch(:, [PF, QF, PT, QT]) = [real(Sf) imag(Sf) real(St) imag(St)];
    if nf(end)>1e-8
        warning('Newton-Raphson did not converge!');
        mpc.success = 0;
        mpc.delta = NaN;
    end
end

