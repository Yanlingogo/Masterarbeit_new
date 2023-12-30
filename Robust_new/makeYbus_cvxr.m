function [Y,Yf,Yt,M,M_line] = makeYbus2(mpc,phase_shift)
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

    Nbus        = size(mpc.bus,1);        %size(,1) get the rows of matrix
    Nbranch     = size(mpc.branch,1);
    idx_fr = mpc.branch(:,F_BUS);
    idx_to = mpc.branch(:,T_BUS);
    branch = mpc.branch;
    
    y_line = branch(:,BR_STATUS)./(branch(:,BR_R)+1i*branch(:,BR_X));
    b_c = branch(:,BR_STATUS).*branch(:,BR_B);

    tap_ratio = zeros(Nbranch,1);
    for i = 1: Nbranch
        if branch(i,TAP) == 0
            tap_ratio(i) = 1;
        else
            tap_ratio(i) = branch(i,TAP);
        end
    end
    tap = tap_ratio.*exp(1i*branch(:,SHIFT));
    y_sh = (mpc.bus(:, GS) + 1j * mpc.bus(:, BS)) / mpc.baseMVA;

    E_fr = sparse(idx_fr, 1:Nbranch, 1, Nbus, Nbranch);
    E_to = sparse(idx_to, 1:Nbranch, 1, Nbus, Nbranch);
    E = E_fr - E_to;
    y_tt = y_line+1i*b_c/2;
    y_ff = y_tt./(tap.*conj(tap));
    y_ft = -y_line./conj(tap);
    y_tf = -y_line./tap;
    Yf = sparse([1:Nbranch;1:Nbranch],[idx_fr;idx_to],[y_ff;y_ft],Nbranch,Nbus);
    Yt = sparse([1:Nbranch;1:Nbranch],[idx_fr;idx_to],[y_tf;y_tt],Nbranch,Nbus);

    Y = E_to*diag(y_tt)*E_to' + E_to*diag(y_tf)*E_fr' + E_fr*diag(y_ft)*E_to'+E_fr*diag(y_ff)*E_fr'+diag(y_sh);

    if phase_shift == false
        Y_cos = sparse([idx_fr; idx_to],[(1:Nbranch)'; (1:Nbranch)'],[y_ft;y_tf],Nbus,Nbranch);
        Y_sin = sparse([idx_fr; idx_to],[(1:Nbranch)'; (1:Nbranch)'],[y_ft;-y_tf],Nbus,Nbranch);
        Y_diag = E_to*diag(y_tt)*E_to'+E_fr*diag(y_ff)*E_fr'+diag(y_sh);

        M = [eye(Nbus), zeros(Nbus), -real(Y_cos), -imag(Y_sin), -real(Y_diag);...
            zeros(Nbus), eye(Nbus),  imag(Y_cos), -real(Y_sin),  imag(Y_diag)];
        M_line=[zeros(Nbranch,2*Nbus)  real(diag(y_ft)) imag(diag(y_ft))  real(diag(y_ff)*E_fr');
                zeros(Nbranch,2*Nbus)  real(diag(y_tf)) -imag(diag(y_tf))   real(diag(y_tt)*E_to');
                zeros(Nbranch,2*Nbus) -imag(diag(y_ft)) real(diag(y_ft)) -imag(diag(y_ff)*E_fr');
                zeros(Nbranch,2*Nbus) -imag(diag(y_tf)) -real(diag(y_tf))  -imag(diag(y_tt)*E_to')];
    elseif phase_shift==true
        Y_plus = E_fr*diag(y_ft.*exp.(-1i*Phi0))+E_to*diag(y_tf.*exp.(1i*Phi0));
        Y_minus = E_fr*diag(y_ft.*exp.(-1i*Phi0))-E_to*diag(y_tf.*exp.(1i*Phi0));
        Y_diag = E_to*diag(y_tt)*E_to'+E_fr*diag(y_ff)*E_fr'+diag(y_sh);

        M = [eye(Nbus), zeros(Nbus), -real(Y_plus), -imag(Y_minus), -real(Y_diag);...
        zeros(Nbus), eye(Nbus),  imag(Y_plus), -real(Y_minus),  imag(Y_diag)];
        M_line = [zeros(Nbranch, 2*Nbus),  real(diag(y_ft.*exp(-1i*Phi0))), imag(diag(y_ft.*exp(-1i*Phi0))),  real(diag(y_ff)*E_fr');
              zeros(Nbranch, 2*Nbus),  real(diag(y_tf.*exp(1i*Phi0))),  -imag(diag(y_tf.*exp(1i*Phi0))),   real(diag(y_tt)*E_to');
              zeros(Nbranch, 2*Nbus), -imag(diag(y_ft.*exp(-1i*Phi0))), real(diag(y_ft.*exp(-1i*Phi0))), -imag(diag(y_ff)*E_fr');
              zeros(Nbranch, 2*Nbus), -imag(diag(y_tf.*exp(1i*Phi0))),  -real(diag(y_tf.*exp(1i*Phi0))),  -imag(diag(y_tt)*E_to')];
    end
end

