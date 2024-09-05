% Example 4: Markov Chain.


%% Setup
K = 5;

tol = 1e-6;
eta = 0.3;

nn = [ 128 256 512 ];

data = zeros(length(nn), 18);

for jjj = 1 : length(nn)
    M = nn(jjj);
    setup_markov;

    %% Unpreconditioned
    maxit = 200;
    tic;
    [x , res, info] = tt_sgmres(A, b, [], ...
        [1 ; 50 * ones(K-1, 1) ; 1], ...
        'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', 300, ...
        'streaming_reorthogonalization', false);
    t_combined = toc;
    numit_combined = info.it;
    res_combined = norm(A*x - b) / norm(b);
    ranks_combined = info.ranks;

    data(jjj, 1) = t_combined;
    data(jjj, 2) = numit_combined;
    data(jjj, 3) = max(ranks_combined);
    data(jjj, 4) = res_combined;

    % Write out the ranks
    writematrix([ (1 : length(info.ranks))', info.ranks(:) ], ...
        sprintf('ex_precond_2_ranks_%d_unprecond.dat', M), ...
        'Delimiter', '\t');

    maxit = 20;
    
    %% Preconditioned version (without maxrank)
    if M < 1024
    RR = cell(1, K);
    for jj = 1 : K
        RR{jj} = R{jj}; d = sum(RR{jj}, 2); 
        RR{jj} = spdiags(d, 0, size(RR{jj}, 1), size(RR{jj}, 2)) - RR{jj};
    end
    % Optimal: q = 16, rkmax = 16
    q = 16; rkmax = inf;
    P = expsum_preconditioner(q, RR, rkmax);
    tic;
    [x , res, info] = tt_sgmres(A, b, [], ...
        [1 ; 50 * ones(K-1, 1) ; 1], ...
        'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf, ...
        'streaming_reorthogonalization', false, 'preconditioner', P);
    t_prec = toc;
    numit_prec = info.it;
    res_prec = norm(A*x - b) / norm(b);
    ranks_prec = info.ranks;

    data(jjj, 5) = t_prec;
    data(jjj, 6) = numit_prec;
    data(jjj, 7) = max(ranks_prec);
    data(jjj, 8) = res_prec;

    % Write out the ranks
    writematrix([ (1 : length(info.ranks))', info.ranks(:) ], ...
        sprintf('ex_precond_2_ranks_%d_precond.dat', M), ...
        'Delimiter', '\t');
    end
    
    %% Precondition version (maxrank set)
    % Setup the preconditioner
    maxrank = 50;
    RR = cell(1, K); l = 0;
    for jj = 1 : K
        RR{jj} = R{jj}; d = sum(RR{jj}, 2); 
        RR{jj} = spdiags(d, 0, size(RR{jj}, 1), size(RR{jj}, 2)) - RR{jj};
    end
    % Optimal: q = 16, rkmax = 16
    q = 32; rkmax = inf;
    P = expsum_preconditioner(q, RR, rkmax);
    tic;
    [x , res, info] = tt_sgmres(A, b, [], ...
        [1 ; maxrank * ones(K-1, 1) ; 1], ...
        'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', maxrank, ...
        'streaming_reorthogonalization', false, 'preconditioner', P, 'tol_preconditioner', 1e-12);
    t_prec2 = toc;
    numit_prec2 = info.it;
    res_prec2 = norm(A*x - b) / norm(b);
    ranks_prec2 = info.ranks;

    data(jjj, 9) = t_prec2;
    data(jjj, 10) = numit_prec2;
    data(jjj, 11) = max(ranks_prec2);
    data(jjj, 12) = res_prec2;

    % Write out the ranks
    writematrix([ (1 : length(info.ranks))', info.ranks(:) ], ...
        sprintf('ex_precond_2_ranks_%d_precond_maxrank30.dat', M), ...
        'Delimiter', '\t');

    %% Precondition version (maxrank set to 80)
    % Setup the preconditioner
    maxrank = 80;
    RR = cell(1, K); l = 0;
    for jj = 1 : K
        RR{jj} = R{jj}; d = sum(RR{jj}, 2); 
        RR{jj} = spdiags(d, 0, size(RR{jj}, 1), size(RR{jj}, 2)) - RR{jj};
    end
    % Optimal: q = 16, rkmax = 16
    q = 32; rkmax = inf;
    P = expsum_preconditioner(q, RR, rkmax);
    tic;
    [x , res, info] = tt_sgmres(A, b, [], ...
        [1 ; maxrank * ones(K-1, 1) ; 1], ...
        'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', maxrank, ...
        'streaming_reorthogonalization', false, 'preconditioner', P, 'tol_preconditioner', 1e-12);
    t_prec2 = toc;
    numit_prec2 = info.it;
    res_prec2 = norm(A*x - b) / norm(b);
    ranks_prec2 = info.ranks;

    data(jjj, 13) = t_prec2;
    data(jjj, 14) = numit_prec2;
    data(jjj, 15) = max(ranks_prec2);
    data(jjj, 16) = res_prec2;

    % Write out the ranks
    writematrix([ (1 : length(info.ranks))', info.ranks(:) ], ...
        sprintf('ex_precond_2_ranks_%d_precond_maxrank80.dat', M), ...
        'Delimiter', '\t');

    %% AMEN
    t_amen = tic;
    x = amen_solve2(A, b, tol);
    t_amen = toc(t_amen);
    res_amen = norm(A*x - b) / norm(b);

    data(jjj, 17) = t_amen;
    data(jjj, 18) = res_amen;

    writematrix([ nn' data ], 'ex_precond_2.dat', 'Delimiter', '\t');
end
