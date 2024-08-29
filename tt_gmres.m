function [x, res, info] = tt_gmres(A, b, tol_stop, maxit)
%


% Factor of which we raise accuracy for the matvec
iap = 1e-2;
sum_tol = tol_stop*iap;

% The basis is represented by a cell array of TT-vectors; the first vector
% is b, normalized to have norm 1.
nrmb = norm(b);
V = { b / nrmb };

H = zeros(1, 0);
j = 1;

info = struct;
info.ranks = [];

% Arnoldi iteration
while true
    % New matvec
    w = A * V{end}; %TODO: aggiungere ROUND ANCHE QUI

    % Reorthogonalization
    for i = 1 : length(V)
        H(i, j) = dot(w, V{i});
        w = round(w - H(i,j) * V{i}, iap * tol_stop);
    end

    H(j+1, j) = norm(w);
    V{j+1} = w / H(j+1,j);

    % Arnoldi check
%     AA = full(A);
%     VV = zeros(size(AA, 1), length(V)); for i = 1 : length(V); VV(:, i) = full(V{i}); end
%     norm(AA * VV(:, 1:end-1) - VV * H)

    rhs = nrmb * eye(length(V), 1);
    y = H \ rhs;
    res(j) = norm(H*y - rhs) / nrmb;

    % Relaxing the tolerance, and truncating with tol / res(j)
    relaxed_tol = iap * tol_stop / res(j);
    V{end} = round(V{end}, relaxed_tol);

    fprintf('TT-GMRES: it = %d, res = %e, rks = %d, trunc_tol = %e\n', ...
        j, res(j), max(rank(V{end})), relaxed_tol);

    info.ranks(j) = max(rank(V{end}));

    if res(j) < tol_stop || j > maxit
        break;
    end

    j = j + 1;
end

info.res = res;
info.it = j;


x = V{1} * y(1);
for j = 2 : length(V) - 1
    x = round(x + V{j} * y(j), sum_tol);
end
end

