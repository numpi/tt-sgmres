% Example 1: PDE convection-diffusion.
% Time comparison between TT SGMRES and TT GMRES 
% over different dimension: 3,4,5,6,7,8,9.
%
max_power_of_2 = 8;
max_power_of_2_gmres = 5;
t_combined = zeros(1, max_power_of_2-4);
t_gmres = zeros(1, max_power_of_2_gmres-4);
res_combined = zeros(1, max_power_of_2-4);
res_gmres = zeros(1, max_power_of_2_gmres-4);
d = 5;
% Dimensions and rank for this test
for power = 5:1:max_power_of_2
    n = 2^power * ones(1, d);

    setup_pde;
    tol = 1e-3;
    eta = 0.3;
    maxit = 1000;
    
    tic;
    [x , res, info] = tt_sgmres(A, b, [], ...
        [1 ; 50 * ones(d-1, 1) ; 1], ...
        'streaming_reorthogonalization', false, ...
        'tol', tol*eta, 'maxit', maxit, 'ktrunc', 1, 'iap', 1e-2, 'max_rank', inf);
    t_combined(power-4) = toc;
    numit_combined = info.it;
    res_combined(power-4) = norm(A*x - b) / norm(b);
    res_combined(power-4)
    t_combined(power-4)
    %return;
    
    if power <= max_power_of_2_gmres
        tic
        [xgmres , resgmres, infogmres] = tt_gmres(A, b, tol, maxit);
        t_gmres(power-4) = toc;
        numit_gmres = infogmres.it;
        res_gmres(power-4) = norm(A*xgmres - b) / norm(b);
        res_gmres(power-4)
        t_gmres(power-4)
    end
writematrix([(5:length(t_gmres)+4)', t_gmres(1:end)'], ...
   'ex_4_times_tt.dat', 'Delimiter', '\t');
writematrix([ (5:length(t_combined)+4)', t_combined(1:end)'], ...
   'ex_4_times_stt.dat', 'Delimiter', '\t');
    
end    

maxdim = length(t_combined) + 4;

plot(5:maxdim, t_combined(1:end)')
hold on
plot(5:maxdim_gmres, t_gmres(1:end)')

legend('TT-SGMRES','TT-GMRES')
title('time comparison')
xlabel('d', 'FontSize', 14)
ylabel('t', 'FontSize', 14)

     %x = amen_solve2(A, b, tol);
