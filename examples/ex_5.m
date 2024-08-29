% Example 4: Markov Chain.

% Time and rank comparison between TT-SGMRES and TT-GMRES with varying 
% dimension d= 4,5,6 and tolerance tol = 1e-4.

% -------------
% d = 4 
%--------------
M = 24;
K = 4;
setup_markov;

tol = 1e-4;
eta = 0.3;
maxit = 500;

tic;
[x , res, info] = tt_sgmres(A, b, [], ...
    [1 ; 50 * ones(K-1, 1) ; 1], ...
    'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf, ...
    'streaming_reorthogonalization', false);
t_combined(1) = toc;
numit_combined = info.it;
res_combined = norm(A*x - b) / norm(b);
[res_combined, t_combined(1)]
ranks_combined = info.ranks;

tic
[xgmres , resgmres, infogmres] = tt_gmres(A, b, tol, maxit);
t_gmres(1) = toc;
numit_gmres = infogmres.it;
res_gmres = norm(A*xgmres - b) / norm(b);
[res_gmres, t_gmres(1)]
ranks_gmres = infogmres.ranks;

writematrix([ (1:length(ranks_combined))', ranks_combined'], ...
   'ex_5_stt_ranks_dim_4.dat', 'Delimiter', '\t');
writematrix([ (1:length(ranks_gmres))',ranks_gmres'], ...
   'ex_5_tt_ranks_dim_4.dat', 'Delimiter', '\t');
    
figure
plot(1:length(ranks_combined), ranks_combined)
hold on 
plot(1:length(ranks_gmres), ranks_gmres)
legend('TT-SGMRES','TT-GMRES')
title('ranks comparison')
xlabel('iteration', 'FontSize', 14)
ylabel('max tt-rank', 'FontSize', 14)

%--------------- 
% d = 5
%--------------
M = 24;
K = 5;
setup_markov;

tol = 1e-4;
eta = 0.3;
maxit = 500;

tic;
[x , res, info] = tt_sgmres(A, b, [], ...
    [1 ; 50 * ones(K-1, 1) ; 1], ...
    'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf, ...
    'streaming_reorthogonalization', false);
t_combined(2) = toc;
numit_combined = info.it;
res_combined = norm(A*x - b) / norm(b);
[res_combined, t_combined(2)]
ranks_combined = info.ranks;

tic
[xgmres , resgmres, infogmres] = tt_gmres(A, b, tol, maxit);
t_gmres(2) = toc;
numit_gmres = infogmres.it;
res_gmres = norm(A*xgmres - b) / norm(b);
[res_gmres, t_gmres(2)]
ranks_gmres = infogmres.ranks;

writematrix([ (1:length(ranks_combined))', ranks_combined'], ...
   'ex_5_stt_ranks_dim_5.dat', 'Delimiter', '\t');
writematrix([ (1:length(ranks_gmres))',ranks_gmres'], ...
   'ex_5_tt_ranks_dim_5.dat', 'Delimiter', '\t'); 
    
figure
plot(1:length(ranks_combined), ranks_combined)
hold on 
plot(1:length(ranks_gmres), ranks_gmres)
legend('TT-SGMRES','TT-GMRES')
title('ranks comparison')
xlabel('iteration', 'FontSize', 14)
ylabel('max tt-rank', 'FontSize', 14)


%--------------- 
% d = 6
%--------------
M = 24;
K = 6;
setup_markov;

tol = 1e-4;
eta = 0.3;
maxit = 500;

tic;
[x , res, info] = tt_sgmres(A, b, [], ...
    [1 ; 50 * ones(K-1, 1) ; 1], ...
    'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf, ...
    'streaming_reorthogonalization', false);
t_combined(3) = toc;
numit_combined = info.it;
res_combined = norm(A*x - b) / norm(b);
[res_combined, t_combined(3)]
ranks_combined = info.ranks;

tic
[xgmres , resgmres, infogmres] = tt_gmres(A, b, tol, maxit);
t_gmres(3) = toc;
numit_gmres = infogmres.it;
res_gmres = norm(A*xgmres - b) / norm(b);
[res_gmres, t_gmres(3)]
ranks_gmres = infogmres.ranks;

writematrix([ (1:length(ranks_combined))', ranks_combined'], ...
   'ex_5_stt_ranks_dim_6.dat', 'Delimiter', '\t');
writematrix([ (1:length(ranks_gmres))',ranks_gmres'], ...
   'ex_5_tt_ranks_dim_6.dat', 'Delimiter', '\t');
    

figure
plot(1:length(ranks_combined), ranks_combined)
hold on 
plot(1:length(ranks_gmres), ranks_gmres)
legend('TT-SGMRES','TT-GMRES')
title('ranks comparison')
xlabel('iteration', 'FontSize', 14)
ylabel('max tt-rank', 'FontSize', 14)


% save times

writematrix([[4,5,6]', t_combined', t_gmres'], ...
   'ex_5_time.dat', 'Delimiter', '\t');    