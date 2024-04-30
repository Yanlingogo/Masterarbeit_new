function [vert,time] = Estimation_vertexes(mpc, base_point)

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
    
    %% parameters 
    tic;
    id_bus      = mpc.bus(:,BUS_I);
    id_gen      = mpc.gen(:,GEN_BUS);
    Nbus        = size(mpc.bus,1);
    Ngen        = numel(id_gen);
    Nbranch     = size(mpc.branch,1);
    
    % idx
    id_slack =  find(mpc.bus(:,BUS_TYPE) == REF);
    id_nslack  = find(mpc.bus(:,BUS_TYPE) ~= REF);
    Nslack = numel(id_slack);
    id_gen_nslack = find(id_gen ~= id_slack);
    id_gen_slack = find(id_gen == id_slack);
    
    baseMVA     = mpc.baseMVA;              % baseMVA
    Umax        = mpc.bus(:,VMAX).^2;                                            
    Umin        = mpc.bus(:,VMIN).^2;
    Phimax      = mpc.branch(:,ANGMAX)/180*pi;
    Phimin      = mpc.branch(:,ANGMIN)/180*pi;
    Pgmin       = mpc.gen(:,PMIN)/baseMVA; Pgmin(id_gen_slack) = -inf;  
    Qgmin       = mpc.gen(:,QMIN)/baseMVA; Qgmin(id_gen_slack) = -inf;
    Pgmax       = mpc.gen(:,PMAX)/baseMVA; Pgmax(id_gen_slack) = inf;
    Qgmax       = mpc.gen(:,QMAX)/baseMVA; Qgmax(id_gen_slack) = inf;
    Fmax        = mpc.branch(:,RATE_A)/baseMVA;
    
    
    Pd          = mpc.bus(:,PD)/baseMVA;  
    Qd          = mpc.bus(:,QD)/baseMVA;
    Cg          = sparse(id_gen,1:Ngen,ones(Ngen,1),Nbus,Ngen);
    
    % branch info
    branch_r   = mpc.branch(:,BR_R);
    branch_x   = mpc.branch(:,BR_X);
    R   = diag(mpc.branch(:,BR_R));
    X   = diag(mpc.branch(:,BR_X));
    Z2  = R.^2 + X.^2; 
    % 
    [Ybus, Yf, Yt] = makeYbus(mpc);
    Gbus           = real(Ybus);
    Bbus           = imag(Ybus);
    Gf             = real(Yf);
    Bf             = imag(Yf);
    Gt             = real(Yt);
    Bt             = imag(Yt);
    from_bus       = mpc.branch(:, F_BUS);                         
    to_bus         = mpc.branch(:, T_BUS);  
    Cf             = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
    Ct             = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);
    C              = Cf - Ct;
    
    B              = (Cf + Ct)';
    A              = [zeros(Nbranch,1) eye(Nbranch)]*B - eye(Nbranch);
    D              = (eye(Nbranch) - A)\eye(size(A));
    DR             = (eye(Nbranch) - A)\(A*R);
    DX             = (eye(Nbranch) - A)\(A*X);
    DX_p = max(DX,zeros(size(DX))).*(DX>0);
    DX_m = min(DX,zeros(size(DX))).*(DX<0);
    
    Mp = D'*R*D;
    Mq = D'*X*D;

    P_p = -Pgmax(id_gen_nslack)/sum(Pgmax(id_gen_nslack));
    Q_p = -Qgmax(id_gen_nslack)/sum(Qgmax(id_gen_nslack));

    %% load data from base_point

    pcc_0    = base_point.pcc_0;
    z_0      = base_point.z_0;
    jac_z    = base_point.jac_z;
    jac_u    = base_point.jac_u;
    branch_r = base_point.branch_r;
    branch_x = base_point.branch_x;
    from_bus = base_point.from_bus;
    id_slack = base_point.id_slack;
    C        = base_point.C;
    Ct       = base_point.Ct;
    Cg_ns    = base_point.Cg_ns;
    Pd       = base_point.Pd;
    Qd       = base_point.Qd;
    u_0      = base_point.u_0;
    %% set the variables 
    import casadi.*
    P_pcc = SX.sym('P_pcc',1,1);
    Q_pcc = SX.sym('Q_pcc',1,1);
    x = vertcat(P_pcc,Q_pcc);
    %% lower & upper bounds
    P_pcc_max = sum(Pd) - sum(Pgmin(id_gen_nslack));
    P_pcc_min = -(sum(Pgmax(id_gen_nslack)) - sum(Pd));
    Q_pcc_max = sum(Qd) - sum(Qgmin(id_gen_nslack));
    Q_pcc_min = -(sum(Qgmax(id_gen_nslack)) - sum(Qd));
    lbx       = vertcat(P_pcc_min,Q_pcc_min);
    ubx       = vertcat(P_pcc_max,Q_pcc_max);
    %% initial state x0
    x0      = pcc_0;
    %% Estimation based on PCC
    delta_pcc = x' - pcc_0;

    pg_estimation = P_p.*delta_pcc(1);
    qg_estimation = Q_p.*delta_pcc(2);

    power_step = [pg_estimation;qg_estimation];

    d_z_pred = - jac_z\jac_u*(power_step-u_0);
    
    z_pred   = z_0 + d_z_pred;
    U_pred   = z_pred(1:Nbus);
    Pij_pred = z_pred(Nbus+1:Nbus+Nbranch);
    Qij_pred = z_pred(Nbus+Nbranch+1:Nbus+2*Nbranch);
    p_pcc    = z_pred(Nbus+2*Nbranch+1); 
    q_pcc    = z_pred(Nbus+2*Nbranch+2);
    L_pred   = z_pred(3*Nbus+1:end);
    e_st     = sparse(id_slack,1,1,Nbus,1);
    g_pred = [e_st'*U_pred-1;
              C*U_pred - 2*(R*Pij_pred+X*Qij_pred) + Z2*L_pred;
              C'*Pij_pred + Ct'*R*L_pred - e_st*p_pcc-Cg_ns*pg_estimation+Pd;
              C'*Qij_pred + Ct'*X*L_pred - e_st*q_pcc-Cg_ns*qg_estimation+Qd;
                L_pred .* Cf*U_pred - Pij_pred.^2  -  Qij_pred .^2];
              % L_pred - (Pij_pred.^2+Qij_pred.^2)./U_pred(from_bus,:)];
    d_z_cor = -jac_z\g_pred;
    z_comp = z_0 + d_z_pred + d_z_cor;
    
    U_comp = z_comp(1:Nbus);

    voltage_lbounds = Umin(id_nslack) - U_comp(id_nslack);
    voltage_ubounds = U_comp(id_nslack) - Umax(id_nslack);

    %% Problem formulation
    gfun = vertcat(voltage_lbounds,voltage_ubounds);
    lbg = -inf(2*(Nbus-1),1);
    ubg = zeros(2*(Nbus-1),1);
    
 
    % objective
    obj_p = P_pcc;
    obj_q = Q_pcc;
    
    %% solver options
    
    % tolerance
    tol        = 1e-6;
    options.ipopt.tol             = tol;
    options.ipopt.constr_viol_tol = tol;
    options.ipopt.compl_inf_tol   = tol;
    options.ipopt.acceptable_tol  = tol;
    options.ipopt.acceptable_constr_viol_tol = tol;
    options.ipopt.print_level = 5;
    % options.ipopt.grad_f = fgrad;
    options.print_time        = 5;
    % options.ipopt.max_iter    = 200;

    %% sampling points
    
    c1c2 = [1 0; -1 0; 0 1; 0 -1];
    vert = zeros(4,2);
    for i = 1:4
        f_samp   = c1c2(i,:)*[obj_p;obj_q]; % p_{k,l}, q_{k,l} of PCC
        objective = f_samp;
        nlp = struct('x',x,'f',objective,'g',gfun);
        S   = nlpsol('solver','ipopt', nlp,options);
        sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                'lbx', lbx, 'ubx', ubx);
        xopt= full(sol.x);
        vert(i,:) = [xopt(1), xopt(2)];
    end
    vert = unique_coor(vert, 1e-4);
    % order the vertices
    [vert, order] = order_points(vert);
    % find the normal vectors
    Normals = Norm_vector(vert);
    
    % refine the vertices
    H = 100;
    vert_c = sum(vert,1)/size(vert,1);
    
    while any(H >= 0.1)
        H = zeros(size(Normals,1),1);
        for i = 1:size(Normals,1)
            f_ext = Normals(i,:)*([obj_p;obj_q]-vert_c');
            objective = f_ext;
            nlp = struct('x',x,'f',-objective,'g',gfun);
            S   = nlpsol('solver','ipopt', nlp,options);
            sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                    'lbx', lbx, 'ubx', ubx);
            xopt= full(sol.x);
            obj_p_opt = xopt(1);
            obj_q_opt = xopt(2);
            
            H(i) = Normals(i,:)*[obj_p_opt;obj_q_opt]-1;
            if  H(i) >= 1e-4 % no movement outwards
                vert(end+1,:) = [obj_p_opt,obj_q_opt];
            end
        end
        % update the central point and normal vector
        vert_c = sum(vert,1)/size(vert,1);
        vert = unique_coor(vert, 1e-4);
        [vert, order] = order_points(vert);
        Normals_new = Norm_vector(vert);
        Normals = unique_normal(Normals, Normals_new,1e-4);
    end
    toc

    [vert, remaining_index] = removeRedun(vert, 1e-6);
    time = toc;
end

%% remove the reduntant vertexes
function [filteredPoints, nRedundant] = removeRedun(points, tolerance)
    % 计算凸包
    k = convhull(points(:,1), points(:,2));
    hullPoints = points(k, :);
    
    % 初始化一个逻辑向量，标记所有点最初都不是冗余的
    nRedundant = true(size(points, 1), 1);
    
    % 对于每个点，检查它是否在凸包的边上（考虑误差）
    for i = 1:size(points, 1)
        point = points(i, :);
        otherpoints = unique_normal(point,points,tolerance);
        otherpoints = [otherpoints;otherpoints(1,:)];
        for j = 1:size(otherpoints,1)
            dist = vectorizedDistanceToLines(otherpoints,point);
            if any(dist <= 1e-6)
                idx = find(dist <= 1e-6);
                if numel(idx) == 1
                    p1 = otherpoints(idx,:);
                    p2 = otherpoints(idx+1,:);
                    nRedundant(i) = ~isPointBetween(point,p1,p2);
                end
            end
        end
    end
    
    % 凸包上的点不考虑为冗余
    % isRedundant(k) = false;
    
    % 过滤掉标记为冗余的点
    filteredPoints = points(nRedundant, :);
end

function distances = vectorizedDistanceToLines(points, singlePoint)
    % 确保points是一个闭环
    if ~isequal(points(1,:), points(end,:))
        points(end+1,:) = points(1,:);
    end
    
    % 计算所有线段的方向向量
    directionVectors = diff(points, 1, 1);
    
    % 计算点到每个线段起点的向量
    pointVectors = singlePoint - points(1:end-1, :);
    
    % 利用叉积公式计算垂直距离
    crossProd = cross([directionVectors zeros(size(directionVectors, 1), 1)], [pointVectors zeros(size(pointVectors, 1), 1)], 2);
    distances = abs(crossProd(:,3)) ./ vecnorm(directionVectors, 2, 2);
end

function isBetween = isPointBetween(point, lineStart, lineEnd)
    % 检查点是否位于线段的两个端点之间
    dot1 = dot(point - lineStart, lineEnd - lineStart);
    dot2 = dot(point - lineEnd, lineStart - lineEnd);
    isBetween = dot1 > 0 && dot2 > 0;
end
