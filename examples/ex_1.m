% Example 1: PDE convection-diffusion.
% Time comparison between TT SGMRES and TT GMRES 
% over different dimension: 3,4,5,6,7,8,9.
%
maxdim = 9;
maxdim_gmres = 6;
t_combined = zeros(1, maxdim-2);
t_gmres = zeros(1, maxdim_gmres-2);
res_combined = zeros(1, maxdim-2);
res_gmres = zeros(1, maxdim_gmres-2);
% Dimensions and rank for this test
for d = 3:1:maxdim
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
    t_combined(d-2) = toc;
    numit_combined = info.it;
    res_combined(d-2) = norm(A*x - b) / norm(b);
    res_combined(d-2)
    t_combined(d-2)
    %return;
    
    if d <= maxdim_gmres
        tic
        [xgmres , resgmres, infogmres] = tt_gmres(A, b, tol, maxit);
        t_gmres(d-2) = toc;
        numit_gmres = infogmres.it;
        res_gmres(d-2) = norm(A*xgmres - b) / norm(b);
        res_gmres(d-2)
        t_gmres(d-2)
    end
writematrix([(3:length(t_gmres)+2)', t_gmres(1:end)'], ...
   'ex_1_times_tt.dat', 'Delimiter', '\t');
writematrix([ (3:length(t_combined)+2)', t_combined(1:end)'], ...
   'ex_1_times_stt.dat', 'Delimiter', '\t');
    
end    

plot(3:maxdim, t_combined(1:end)')
hold on
plot(3:maxdim_gmres, t_gmres(1:end)')

legend('TT-SGMRES','TT-GMRES')
title('time comparison')
xlabel('d', 'FontSize', 14)
ylabel('t', 'FontSize', 14)

     %x = amen_solve2(A, b, tol);