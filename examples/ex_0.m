d = 4;    
n = 34 * ones(1, d);

setup_pde;
tol = 1e-5;
eta = 0.3;
maxit = 500;

tic;
[x , res, info] = tt_gmres_vanilla(A, b, [],...
    'tol', tol*eta, 'maxit', maxit, 'check_residual', true, ...
    'streaming_reorthogonalization', false, ...
    'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf);
t_vanilla = toc;
numit_combined = info.it;
res_vanilla = norm(A*x - b) / norm(b);
[res_vanilla, t_vanilla]
ranks_vanilla = info.ranks;
%return;

writematrix([ (1:info.it)', info.full_residual',res'], ...
   'ex_0_res_vs_sres.dat', 'Delimiter', '\t');

semilogy(1:info.it, [info.full_residual', res'])
hold on 
legend('res', 'sketched res')
title('Sketching reliability')
xlabel('iteration', 'FontSize', 14)
ylabel('max tt-rank', 'FontSize', 14)