function P = expsum_preconditioner(q, R, rkmax)
%EXPSUM_PRECONDITIONER Construct the preconditioner for a Kronecker sum.
%
% This function returns a function handle P such that y = P(x, tol) applies
% the preconditioner y = P \ x; the result of the preconditioner is 
% truncated by a sketching low-rank approximation algorithm with the 
% prescribed relative tolerance.
%
% The preconditioner is computed by an exponential sum approximation to the
% inverse of the Kronecker sum with coefficients R{i}, of length 2*q + 1.

% Coefficients of the exponential sum
eta = pi / sqrt(q);
k = q : -1 : -q;
beta = exp(k * eta);
alpha = eta * beta;

% Dimensionality of the tensor
d = length(R);

% Heuristic for the target ranks for the randomization procedure, 
% if not chosen by the user
if ~exist('rkmax', 'var')
    rkmax = inf;
end

% Precomputed matrix exponentials (if the R{i} are sufficiently small)
E = cell(2*q + 1, d);
for s = 1 : d
    for j = 1 : 2*q + 1
        E{j, s} = expm(-beta(j) * R{s});
    end
end

P = @(X, tol) evaluate_expsum(alpha, beta, E, X, tol, rkmax);

end

function SS = evaluate_expsum(alpha, beta, E, X, tol, rkmax)
    le = length(beta); % Length of the exponential sum

    % Generate sketches for the STTA algorithm
    d = ndims(X);
    n = size(X);
    if rkmax == inf
        rkmax = max(rank(X)) * 2 + d;
    end
    rk = [1 ; ones(d-1, 1) * rkmax ; 1]; l_over = 20;
    [Y, Z, Ycell, Zcell] = STTA_generate_tt_sketches(n, rk, [1 ; l_over+rk(2:end-1); 1]);
    
    Omega = cell(1, le);
    Psi   = cell(1, le);

    for j = 1 : le
        M = X;
        for s = 1 : ndims(X)
            % M = ttm(M, ndims(X) - s + 1, expm(-beta(j) * R{s}));
            M = ttm(M, ndims(X) - s + 1, E{j,s});
        end
             
        [Omega{j}, Psi{j}] = STTA_contractions(M, Y, Z, Ycell, Zcell);
    end

    SS = STTA_sum_recovery(Psi, Omega, alpha, Y.d, Y.n, Y.r);
    SS = round(SS, tol);

    %YY = alpha(1) * Y{1};
    %for j = 2 : length(Y)
    %    YY = YY + alpha(j) * Y{j};
    %end
    %Y = YY;
end

