markov_tol = 1e-8;

rng(123)

% Rate at which we trigger the interaction
a = .1;

% Generate the transition matrices of the local systems
R = cell(1, K);

for i = 1 : K
    % Transition probabilities to the next state
    next = 1 + rand(1, M);
    
    % Transition probabilities to the previous state
    prev = 1 + rand(1, M);

    RR = diag(next(1:M-1), 1) + diag(prev(2:M), -1);
    RR(1, M) = prev(1);
    RR(M, 1) = next(M);
    R{i} = RR;
end

% Build the sparse version of R
sR = cell(1, K);
for i = 1 : K
    sR{i} = sparse(R{i});
end

L = cell(1, K-1);

% Construct the links between submodels
for i = 1 : K - 1
    % Create the link from R{i} to R{i+1}
    W = zeros(M); W(M-1, M) = a;
    L{i} = W;
end

% Sparse version of L
sL = cell(1, K);
for i = 1 : K - 1
    sL{i} = sparse(L{i});
end

% Spectral norm upper bound, updated as we sum all contributions
nrmA = 0;

% Build the matrix Q in TT-format, only feasible for small matrices
Q = tt_matrix(tt_zeros(M^2 * ones(1, K)), M * ones(1, K), M * ones(1, K));

for i = 1 : K
    % Generate the matrix I \otimes ... \otimes I \otimes R{i} \otimes I
    % ... \otimes I
    Rl = tkron(tt_eye(M * ones(1,K-i)), tt_matrix(R{i}));
    Rl = tkron(Rl, tt_eye(M * ones(1, i-1)));
    Rl = purify(Rl);
    Q = Q + Rl;
    nrmA = nrmA + sqrt( norm(R{i}, 1) * norm(R{i}, inf) );
end

for i = 1 : K - 1
    Il = tt_eye(M * ones(1, i-1));
    Ir = tt_eye(M * ones(1, K - i - 1));
    % N = tt_matrix(diag(ones(M-1,1), 1));
    N = zeros(M); N(end-1, end) = 1;
    N = tt_matrix(N);

    Wl = tkron(tt_matrix(L{i}), Il);
    Wl = tkron(N, Wl);
    Wl = tkron(Ir, Wl);
    Wl = purify(Wl);

    Q = Q + Wl;
    nrmA = nrmA + 1;
end

% Fix the diagonal entries
e = tt_ones(M * ones(1, K));
d = Q * e;
D = diag(d);
Q = round(Q - D, markov_tol);
nrmA = 2 * nrmA;

% Rank-one correction
re = tt_matrix(tt_ones(M^2 * ones(1, K)), M * ones(1, K), M * ones(1, K)) / M^K;

% Corrected matrix
Qe = round(Q - re, markov_tol);
nrmA = nrmA + 1;

A = Qe';
b = -e;
b = b / norm(b);

Afun = @(x) markovAfun(sR, sL, d, x);

function r = markovAfun(sR, sL, d, x)
    r = 0 * x;
    K = length(sR);
    M = size(sR{1}, 1);
    
    N = sparse(M - 1, M, 1, M, M);    

    for i = 1 : K
        r = r + ttm(x, K-i+1, sR{i});
    end

    for i = 1 : K - 1
        s = ttm(x, K-i+1, N);
        s = ttm(s, K-i, sL{i});
        r = r + s;
    end

    r = r - d .* x;
    r = r - sum(x) / M^K;
end