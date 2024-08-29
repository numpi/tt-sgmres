% Example Preconditioned 1: PDE


%% Setup
d = 5;
nn = [ 128 256 512 1024 ];

data = zeros(length(nn), 11);

for jjj = 1 : length(nn)

    n = nn(jjj) * ones(1, d);
    setup_pde;
    
    tol = 1e-8;
    eta = 0.1;
    maxit = 200;
    
    %% Unpreconditioned
    tic;
    [x , res, info] = tt_sgmres(A, b, [], ...
        [1 ; 50 * ones(d-1, 1) ; 1], ...
        'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf, ...
        'streaming_reorthogonalization', false);
    t_combined = toc;
    numit_combined = info.it;
    res_combined = norm(A*x - b) / norm(b);
    ranks_combined = info.ranks;

    data(jjj, 1) = t_combined;
    data(jjj, 2) = numit_combined;
    data(jjj, 3) = max(ranks_combined);

    % Write out the ranks
    writematrix([ (1 : length(info.ranks))', info.ranks(:) ], ...
        sprintf('ex_precond_1_ranks_%d_unprecond.dat', n(1)), ...
        'Delimiter', '\t');

    maxit = 20;
    
    %% Preconditioned version (without maxrank)
    RR = cell(1, d); 
    for jj = 1 : d
        RR{jj} = (K * L + w * D)';
    end
    % Optimal: q = 16, rkmax = 16
    q = 8; rkmax = inf;
    P = expsum_preconditioner(q, RR, rkmax);
    tic;
    [x , res, info] = tt_sgmres(A, b, [], ...
        [1 ; 50 * ones(d-1, 1) ; 1], ...
        'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf, ...
        'streaming_reorthogonalization', false, 'preconditioner', P);
    t_prec = toc;
    numit_prec = info.it;
    res_prec = norm(A*x - b) / norm(b);
    ranks_prec = info.ranks;

    data(jjj, 4) = t_prec;
    data(jjj, 5) = numit_prec;
    data(jjj, 6) = max(ranks_prec);

    % Write out the ranks
    writematrix([ (1 : length(info.ranks))', info.ranks(:) ], ...
        sprintf('ex_precond_1_ranks_%d_precond.dat', n(1)), ...
        'Delimiter', '\t');
    
    %% Precondition version (maxrank set)
    % Setup the preconditioner
    maxrank = 30;
    RR = cell(1, d); l = 0;
    for jj = 1 : d
        RR{jj} = (K * L + w * D)';
    end
    % Optimal: q = 16, rkmax = 16
    q = 8; rkmax = inf;
    P = expsum_preconditioner(q, RR, rkmax);
    tic;
    [x , res, info] = tt_sgmres(A, b, [], ...
        [1 ; maxrank * ones(d-1, 1) ; 1], ...
        'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', 30, ...
        'streaming_reorthogonalization', false, 'preconditioner', P);
    t_prec2 = toc;
    numit_prec2 = info.it;
    res_prec2 = norm(A*x - b) / norm(b);
    ranks_prec2 = info.ranks;

    data(jjj, 7) = t_prec2;
    data(jjj, 8) = numit_prec2;
    data(jjj, 9) = max(ranks_prec2);

    % Write out the ranks
    writematrix([ (1 : length(info.ranks))', info.ranks(:) ], ...
        sprintf('ex_precond_1_ranks_%d_precond_maxrank.dat', n(1)), ...
        'Delimiter', '\t');

    %% AMEN
    t_amen = tic;
    x = amen_solve2(A, b, tol);
    t_amen = toc(t_amen);
    res_amen = norm(A*x - b) / norm(b);

    data(jjj, 10) = t_amen;
    data(jjj, 11) = res_amen;

    writematrix([ nn', data ], 'ex_precond_1.dat', 'Delimiter', '\t');
end
