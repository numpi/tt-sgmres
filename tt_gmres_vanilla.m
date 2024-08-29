function [x, res, info] = tt_gmres_vanilla(A, b, hS, varargin)
% Aggiungere documentazione
p = inputParser;

addOptional(p, 'tol', 1e-6);
addOptional(p, 'ktrunc', 1);
addOptional(p, 'maxit', 200);
addOptional(p, 'iap', 1e-2);
addOptional(p, 'check_residual', false);
addOptional(p, 'kronecker_embedding', false);
addOptional(p, 'preconditioner', []);
addOptional(p, 'max_rank', inf);
% restart threshold
addOptional(p, 'restart_tol', 1e50);

addOptional(p, 'streaming_reorthogonalization', true);

% This can either be: whitening, sketched-ls
addOptional(p, 'whitening', false);

parse(p, varargin{:});


tol_stop = p.Results.tol;
ktrunc = p.Results.ktrunc;
maxit = p.Results.maxit;
check_residual = p.Results.check_residual;
kronecker_embedding = p.Results.kronecker_embedding;
whitening = p.Results.whitening;
preconditioner = p.Results.preconditioner;
restart_tol = p.Results.restart_tol;
max_rank = p.Results.max_rank;
streaming_reorthogonalization = p.Results.streaming_reorthogonalization;

iap = p.Results.iap;

sum_tol = iap*tol_stop;
n = size(b); d = length(n);
info = struct;
info.ranks = [];

if check_residual; info.full_residual = []; end

% If hS is the special field [], we automatically select Khatri-Rao
% sketching based on the maximum number of iterations. 
if isempty(hS)
    m = round(maxit * 2);
    
    if kronecker_embedding
        m_kron = 10; % TODO: find a suitable choice;
        S = TT_kronecker_embedding(ones(1, d) * m_kron, n, d);
        hS = @(v) Kron_sketching(v, S, m);
    else
        S = TT_khatri_rao_embedding(round(m), n, d);
        hS = @(v) KR_sketching(v, S);
    end
    
end

% Convert to A to a function handle, if it's not already in that form
if ~isa(A, 'function_handle')
    A = @(x) A*x;
end

% Factor of which we raise accuracy for the matvec
iap = p.Results.iap;

% Initial solution, and computationa of the starting residual
x = tt_zeros(size(b), ndims(b));
r = round(b - A(x), iap * tol_stop, max_rank);

% Restart counter
restarts = 0;

% Oversampling parameter: we oversample on Z to avoid artificially
% inflating the rank of the representation
%l_over = 20; % maybe d*k (k ~ 3/4 could go)
%[Y, Z, Ycell, Zcell] = STTA_generate_tt_sketches(b.n, rk, [1 ; l_over+rk(2:end-1); 1]);
    
% Outer loop on restarts
while true
    % Cell arrays used to store the sketches of the basis vectors  
    %Omega = {}; Psi = {};

    % We use SAV to store the sketched S*A*V
    if ~whitening; SAV = []; end    

    % The basis is represented by a cell array of TT-vectors; the first vector
    % is b, normalized to have norm 1.
    Sb  = hS(r);
    RSV = norm(Sb);
    
    % Norm of the RHS (which is the residual in the restart loop)
    nrmb = norm(r);
    
    V = { r / nrmb };
    %[Omega{1}, Psi{1}] = STTA_contractions(V{1}, Y, Z, Ycell, Zcell);
    
    % Iteration counter
    j = 1;
    % Sketched basis
    if whitening
        SV  = hS(V{1}) / RSV(1,1);
    end
    
    e1 = eye(maxit+3,1);
    H = zeros(maxit+2,maxit+1);
    
    % Inner Arnoldi iteration
    while true
        if ~isempty(preconditioner)
            w = preconditioner(V{end}, 1e-12);
            w = A(w);
        else
            w = A(V{end});
        end

        % w = round(w, iap * tol_stop);
            
        % Update the sketched SAV matrix: note that it seems fairly
        % relevant to do this _before_ rounding, because otherwise we lose
        % some information on the subspace.
        if ~whitening
            % w = round(w, iap * tol_stop);
            SAV = [SAV, hS(w)];
        end

        % w = round(w, iap * tol_stop, max_rank);        
    
        % Local reorthogonalization -- we do simple Gram-Schimdt to take
        % advantage of the sketching procedure and avoid intermediate
        % roundings. We may consider doing it twice for increased accuracy,
        % but we do not worry too much about loss of orthogonality here.
        for i = 1 : min(ktrunc, length(V))
            H(j-i+1,j) = dot(w, V{end-i+1});

            if ~streaming_reorthogonalization
                w = round(w - H(j-i+1,j) * V{end-i+1}, iap * tol_stop, max_rank);
                % w = w - H(j-i+1,j) * V{end-i+1};
            end
        end

        if streaming_reorthogonalization
            w = streaming_reorthogonalization_step(V, ktrunc);
        end

    
        nrm = norm(w);
    
        H(j+1,j) = nrm;
        V{j+1}   = w / nrm;
        %[Omega{j+1}, Psi{j+1}] = STTA_contractions(V{j+1}, Y, Z, Ycell, Zcell);
    
        % In the whitened version, we keep the matrix R with the 
        % coefficients that would make V * R orthogonal with respect to the
        % sketched inner product.
        if whitening
            sw = hS(V{j+1});    
            for i = 1 : length(V)-1
                RSV(i,j+1) = dot(sw, SV(:, i));  
                sw = sw - RSV(i,j+1) * SV(:,i);
            end

            RSV(j+1,j+1)=norm(sw);
            SV = [SV, sw / RSV(j+1,j+1)];
        end

        if whitening
            V{j} = round(V{j}, iap * tol_stop, max_rank);

            ej=zeros(j,1); ej(j)=1;
            coeff=RSV(1:j,1:j)*H(1:j,1:j)/RSV(1:j,1:j) + RSV(1:j,j+1)*(ej'*H(j+1,j)/RSV(j,j));
            coeff=[coeff; zeros(1,j-1) H(j+1,j)];
            rhs_proj=RSV(1,1)*nrmb*e1(1:j+1);
            y = coeff\rhs_proj;
            res(j) = norm(coeff * y - rhs_proj) / RSV(1,1);
        else
            l = 0; % norm(SAV, 'fro') * 1e-8;
            [Q, R] = qr([ SAV ; l * eye(size(SAV, 2)) ], 0);
            y = R \ (Q(1:size(Sb, 1), :)' * Sb);
            res(j) = norm(SAV * y - Sb) / norm(Sb);
            % norm(y)
            % cond(R)

            %xtmp = build_solution(Psi, Omega, y, preconditioner, x, Y.d, Y.n, Y.r, iap, tol_stop);
            %restrue = norm(A(xtmp) - b) / norm(b);
            %[res(j), restrue, restrue / res(j) ]
        end
       
        % For debugging only: compute the actual residual when requested
        if check_residual
            if whitening; yl = RSV(1:j,1:j) \ y; else; yl = y; end
            
            %xtest = build_solution(Psi, Omega, yl, preconditioner, x, Y.d, Y.n, Y.r, iap, tol_stop);

            xtest = V{1} * y(1);
            for jj = 2 : length(V) - 1
                xtest = round(xtest + V{jj} * y(jj), sum_tol);
            end
            

            %xtest = new_TTsum_Randomize_then_Orthogonalize(V(1:end-1), yl, iap * tol_stop);
            info.full_residual(j) = norm(b - A(xtest))/norm(b);
            fullres = sprintf('%e', info.full_residual(j));
        else
            fullres = 'N/A';
        end

        % the restart is computed whenever things get too ill-conditioned
        if whitening
            cond2check=cond(RSV(1:j,1:j));
        else
            cond2check=cond(SAV);
        end
        fprintf('TT-SGMRES: it = %d, res = %e, rks = %d, full-res = %s, cnd = %e\n', ...
            j, res(j), max(rank(V{end})), fullres, cond2check);

        info.ranks(j) = max(rank(V{end}));
        if res(j) < tol_stop %|| j > maxit || cond2check > restart_tol
            break;
        end
    
        j = j + 1;
    end

    if whitening
        y = RSV(1:j,1:j) \ y;
    end
    
    %x = build_solution(Psi, Omega, y, preconditioner, x, Y.d, Y.n, Y.r, iap, tol_stop);
    
    x = V{1} * y(1);
    for jj = 2 : length(V) - 1
        x = round(x + V{jj} * y(jj), sum_tol);
    end
        
    if res(j) < tol_stop || true
        info.it = j;
        break;
    else
        % Prepare for restarting, we truncate the residual aggressively
        % based on the convergence history.
        r = b - A(x);
        r = round(r, iap * tol_stop);

        tol_stop = tol_stop * norm(b) / norm(r);

        fprintf('TT-SGMRES: | Restarting at iteration j |\n');
        fprintf('           | New tol_stop = %e         |\n', tol_stop);

        restarts = restarts + 1;
    end
    
end

end

function x = build_solution(Psi, Omega, y, preconditioner, x, d, n, r, iap, tol_stop)
xd = STTA_sum_recovery(...
        Psi(1 : length(y)), ...
        Omega(1 : length(y)), y, d, n, r);

    %xd = new_TTsum_Randomize_then_Orthogonalize(V(1:end-1), y, iap * tol_stop);
    if ~isempty(preconditioner)
        xd = preconditioner(xd, iap * tol_stop);
    end

    x = round(x + xd, iap * tol_stop);
    % x = x + xd;
end

function w = streaming_reorthogonalization_step(V, ktrunc)
    il = min(ktrunc, length(V));
    n_addends = il+1;
    N = V{end-il+1}.d;
    r1 =  V{end-il+1}.r;
    rs = zeros(N+1,n_addends);
    rs(:,1) = r1;
    for jj=2:n_addends-1
       rs(:,jj) = V{end-il+jj}.r;
    end 
    rs(:, n_addends) = w.r;

    % heuristic for the target ranks for the randomization procedure
    ro = max(rs, [], 2) * 2 + n_addends;
    ro(1) = 1;
    ro(N + 1) = 1;
     [Y_prod, Z_prod, Ycell_prod, Zcell_prod] = ...
         STTA_generate_tt_sketches(b.n, ro, [1 ; l_over+ro(2:end-1); 1]);
     Omega_prod = cell(1, n_addends);
     Psi_prod = cell(1, n_addends);
     [Omega_prod{1}, Psi_prod{1}] = ...
        STTA_contractions(w, Y_prod, Z_prod, Ycell_prod, Zcell_prod);
     for jj = 2:n_addends
        [Omega_prod{jj}, Psi_prod{jj}] = ...
        STTA_contractions(V{end-n_addends+jj}, Y_prod, Z_prod, Ycell_prod, Zcell_prod);
     end
     w = STTA_sum_recovery(...
         [ Psi_prod(end-n_addends+1:end) ], ...
         [ Omega_prod(end-n_addends+1:end) ], ...
         [1 ; -H(j-il+1:j,j) ], Y_prod.d, Y_prod.n, Y_prod.r);

    w = round(w, iap * tol_stop, max_rank);
end

