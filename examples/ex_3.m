% check how reliable is the sketching by showing both the sketched residual
% and the actual residual.

% Example 3: PDE convection-diffusion.
% TT ranks evolution over iterations, comparison between TT SGMRES and TT GMRES 
%
% Dimensions and rank for this test
for d = 3:2:9
    
    n = 64 * ones(1, d);
    
    setup_pde;
    tol = 1e-6;
    eta = 0.3;
    maxit = 500;
    
    tic;
    [x , res, info] = tt_sgmres(A, b, [], ...
        [1 ; 50 * ones(d-1, 1) ; 1], ...
        'tol', tol*eta, 'maxit', maxit, 'check_residual', true, ...
        'streaming_reorthogonalization', true, ...
        'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf);
    t_combined = toc;
    numit_combined = info.it;
    res_combined = norm(A*x - b) / norm(b);
    [res_combined, t_combined]
    ranks_combined = info.ranks;
    %return;
    
    if d == 3
        writematrix([ (1:info.it)', info.full_residual',res'], ...
           'ex_3_res_vs_sres_dim_3.dat', 'Delimiter', '\t');
    end
    if d == 5
        writematrix([ (1:info.it)', info.full_residual',res'], ...
           'ex_3_res_vs_sres_dim_5.dat', 'Delimiter', '\t');
    end
    if d == 7
        writematrix([ (1:info.it)', info.full_residual',res'], ...
           'ex_3_res_vs_sres_dim_7.dat', 'Delimiter', '\t');
    end
    if d == 9
        writematrix([ (1:info.it)', info.full_residual',res'], ...
           'ex_3_res_vs_sres_dim_9.dat', 'Delimiter', '\t');
    end

end

semilogy(1:info.it, [info.full_residual', res'])
hold on 
legend('res', 'sketched res')
title('Sketching reliability')
xlabel('iteration', 'FontSize', 14)
ylabel('max tt-rank', 'FontSize', 14)

 %x = amen_solve2(A, b, tol);