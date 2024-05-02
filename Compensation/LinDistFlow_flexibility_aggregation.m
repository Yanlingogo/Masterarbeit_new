function [vert,dispatch,sol_variables,time] = LinDistFlow_flexibility_aggregation(mpc,model,base_point)

    %%Index setting
    % bus idx
    tic
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
    
    %% check input
    if nargin < 2 
        model = 1;
    end
    %% parameters 
    tic;
    id_bus      = mpc.bus(:,BUS_I);
    id_gen      = mpc.gen(:,GEN_BUS);
    Nbus        = size(mpc.bus,1);
    Ngen        = numel(id_gen);
    Nbranch     = size(mpc.branch,1);
    
    % idx
    id_slack =  find(mpc.bus(:,BUS_TYPE) == REF);
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
    %% set the variables 
    import casadi.*
    U          = SX.sym('U',Nbus,1);
    Pij        = SX.sym('Pij',Nbranch,1);
    Qij        = SX.sym('Qij',Nbranch,1);
    Pg         = SX.sym('Pg',Ngen,1);
    Qg         = SX.sym('Qg',Ngen,1);
    x          = vertcat(U, Pij, Qij, Pg, Qg);
    %% lower & upper bounds
    lbx         = [Umin;-inf(2*Nbranch,1);Pgmin;Qgmin;];       
    ubx         = [Umax; inf(2*Nbranch,1);Pgmax;Qgmax;];
    %% initial state x0
    
    U0          = mpc.bus(:,VM).^2;
    Pg0         = mpc.gen(:,PG)/baseMVA;
    Qg0         = mpc.gen(:,QG)/baseMVA;
    Pij0        = zeros(Nbranch,1);
    Qij0        = zeros(Nbranch,1);
    x0      = vertcat(U0, Pij0, Qij0, Pg0, Qg0);
    %% LinDistFlow
    if (model == 1) ||(model == 2)
        if model == 1
            L_0 = zeros(Nbranch,1);
        elseif model == 2
            L_0 = base_point.L_0;
        end
        % voltage constraints
        volt_eq    = C*U -2*(branch_r .* Pij + branch_x .* Qij) + Z2*L_0;
    
        % power balance
        pf_p_eq    = Cg*Pg - Pd - C'*Pij - Cf'*R*L_0;
        pf_q_eq    = Cg*Qg - Qd - C'*Qij - Cf'*X*L_0;
    
        % Problem formulation
        gfun = vertcat(volt_eq,pf_p_eq,pf_q_eq);
        lbg = zeros(Nbranch+2*Nbus, 1);
        ubg = zeros(Nbranch+2*Nbus, 1);
    end

    %% Linearized DistFlow
    if model == 3
        % date from base point
        jac_z = base_point.jac_z;
        jac_u = base_point.jac_u;
        z_0 = base_point.z_0;
        u_0 = base_point.u_0;
        L_0 = base_point.L_0;
        % Linearization of DistFlow
        L = SX.sym('L',Nbranch,1);
        x = vertcat(x,L);
        x0 = vertcat(x0,L_0);
        lbx = vertcat(lbx,-inf(Nbranch,1));
        ubx = vertcat(ubx,inf(Nbranch,1));

        var_state = vertcat(U,Pij,Qij,Pg(id_gen_slack),Qg(id_gen_slack),L);
        var_deci  = vertcat(Pg(id_gen_nslack),Qg(id_gen_nslack));
        % equalities
        gfun= jac_z*(var_state - z_0) + jac_u*(var_deci-u_0);
        lbg = zeros(Nbranch+3*Nbus,1);
        ubg = zeros(Nbranch+3*Nbus,1);
    end

    
    % objective
    obj_p = Pg(id_gen_slack);
    obj_q = Qg(id_gen_slack);
    
    %% solver options
    
    % tolerance
    tol        = 1e-6;
    options.ipopt.tol             = tol;
    options.ipopt.constr_viol_tol = tol;
    options.ipopt.compl_inf_tol   = tol;
    options.ipopt.acceptable_tol  = tol;
    options.ipopt.acceptable_constr_viol_tol = tol;
    options.ipopt.print_level = 0;
    % options.ipopt.grad_f = fgrad;
    options.print_time        = 0;
    % options.ipopt.max_iter    = 200;

    %% sampling points
    
    c1c2 = [1 0; -1 0; 0 1; 0 -1];
    vert = zeros(4,2);
    dispatch = cell(4,1);
    sol_variables = cell(4,1);
    for i = 1:4
        f_samp   = c1c2(i,:)*[obj_p;obj_q]; % p_{k,l}, q_{k,l} of PCC
        objective = f_samp;
        nlp = struct('x',x,'f',objective,'g',gfun);
        S   = nlpsol('solver','ipopt', nlp,options);
        sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                'lbx', lbx, 'ubx', ubx);
        xopt= full(sol.x);
        obj_p_opt = xopt(Nbus+2*Nbranch+id_gen_slack);
        obj_q_opt = xopt(Nbus+2*Nbranch+Ngen+id_gen_slack);
        vert(i,:) = [obj_p_opt, obj_q_opt];
        dispatch{i} = [xopt(Nbus+2*Nbranch+id_gen_nslack); xopt(Nbus+2*Nbranch+Ngen+id_gen_nslack)];
        sol_variables{i} = xopt;
    end
    tol = 1e-2;
    vert = NFD_unique_coor(vert, tol);
    % order the vertices
    [vert, order] = NFD_order_points(vert);
    dispatch = dispatch(order);
    sol_variables = sol_variables(order);
    % find the normal vectors
    Normals = NFD_Norm_vector(vert);
    
    % refine the vertices
    H = 100;
    vert_c = sum(vert,1)/size(vert,1);
    
    while any(H >= 1 + tol) || any( H <= 1-tol)
        H = zeros(size(Normals,1),1);
        for i = 1:size(Normals,1)
            f_ext = Normals(i,:)*([obj_p;obj_q]-vert_c');
            objective = f_ext;
            nlp = struct('x',x,'f',-objective,'g',gfun);
            S   = nlpsol('solver','ipopt', nlp,options);
            sol = S('x0', x0,'lbg', lbg,'ubg', ubg,...
                    'lbx', lbx, 'ubx', ubx);
            xopt= full(sol.x);
            obj_p_opt = xopt(Nbus+2*Nbranch+id_gen_slack);
            obj_q_opt = xopt(Nbus+2*Nbranch+Ngen+id_gen_slack);
            
            H(i) = abs(Normals(i,:)*([obj_p_opt;obj_q_opt]));
            if  (H(i) >= 1 + tol) || (H(i) <= 1-tol) % no movement outwards
                vert(end+1,:) = [obj_p_opt,obj_q_opt];
                dispatch{end+1} = [xopt(Nbus+2*Nbranch+id_gen_nslack); xopt(Nbus+2*Nbranch+Ngen+id_gen_nslack)];
                sol_variables{end+1} = xopt;
            end
        end
        % update the central point and normal vector
        vert_c = sum(vert,1)/size(vert,1);
        vert = NFD_unique_coor(vert, tol);
        [vert, order] = NFD_order_points(vert);
        dispatch = dispatch(order);
        sol_variables = sol_variables(order);
        Normals_new = NFD_Norm_vector(vert);
        Normals = NFD_unique_normal(Normals, Normals_new,tol);
    end
    toc

    [vert, remaining_index] = removeRedun(vert, 1e-6);
    dispatch = dispatch(remaining_index == 1);
    sol_variables = sol_variables(remaining_index == 1);
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
        otherpoints = NFD_unique_normal(point,points,tolerance);
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
