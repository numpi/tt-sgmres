% Example 2: PDE convection-diffusion.
% TT ranks evolution over iterations, comparison between TT SGMRES and TT GMRES 
%
% Dimensions and rank for this test
d = 6;
n = 64 * ones(1, d);

setup_pde;
tol = 1e-4;
eta = 0.3;
maxit = 500;

tic;
[x , res, info] = tt_sgmres(A, b, [], ...
    [1 ; 50 * ones(d-1, 1) ; 1], ...
    'streaming_reorthogonalization', false, ...
    'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf);
t_combined = toc;
numit_combined = info.it;
res_combined = norm(A*x - b) / norm(b);
[res_combined, t_combined]
ranks_combined = info.ranks;
%return;

tic
[xgmres , resgmres, infogmres] = tt_gmres(A, b, tol, maxit);
t_gmres = toc;
numit_gmres = infogmres.it;
res_gmres = norm(A*xgmres - b) / norm(b);
[res_gmres, t_gmres]
ranks_gmres = infogmres.ranks;

writematrix([ (1:length(ranks_combined))', ranks_combined'], ...
   'ex_2_stt_ranks.dat', 'Delimiter', '\t');
writematrix([ (1:length(ranks_gmres))',ranks_gmres'], ...
   'ex_2_tt_ranks.dat', 'Delimiter', '\t');
writematrix([t_combined, t_gmres], ...
   'ex_2_time.dat', 'Delimiter', '\t');    

plot(1:length(ranks_combined), ranks_combined)
hold on 
plot(1:length(ranks_gmres), ranks_gmres)
legend('TT-SGMRES','TT-GMRES')
title('ranks comparison')
xlabel('iteration', 'FontSize', 14)
ylabel('max tt-rank', 'FontSize', 14)

 %x = amen_solve2(A, b, tol);