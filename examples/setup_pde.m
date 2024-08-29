%
% This script defines all matrices needed for the PDE example in the paper
% ...
%
% You need to define the following variables before calling the script:
%
%  - n: [n1, ..., nd] Size of the tensor

if ~exist('n', 'var')
    error('Please define the vector n before calling setup_pde');
end

if any(n ~= n(1))
    error('Only equal dimensions in all modes are supported for now');
end

% Dimensions and rank for this test
d = length(n);

% Convection-Diffusion equation on [-1, 1]^d
h = 2 / (n(1) + 1); K = 1e-2; w = 5e-2;
L =  (1 / h^2) * (2 * eye(n(1)) - diag(1 * ones(n(1)-1, 1), 1) - diag(1 * ones(n(1)-1, 1), -1));
D =  (1 / h) * (diag(1 * ones(n(1)-1, 1), 1) - eye(n(1)));
AA = tt_matrix(K * L + w * D);

% Build the TT-matrix A, which in our case is just the d-dimensional
% Laplace operator. 
% h = 1 / (n1 - 1); AA = 2 * eye(n1) - diag(0.9 * ones(n1-1, 1), 1) - ...
%     diag(1.1 * ones(n1-1, 1), -1);
% AA = tt_matrix(AA / h^2);

for j = 1 : d
    if j == 1; Ax = AA; else; Ax = tt_eye(n(1)); end

    for i = 2 : d
        if i == j
            Ax = tkron(Ax, AA);
        else
            Ax = tkron(Ax, tt_eye(n(1)));
        end
    end

    if j == 1
        A = Ax;
    else
        A = round(A + Ax, 1e-8);
    end
end

% Generate the RHS as exp(-10*(x1.^2 + ... + xd.^2))
x = linspace(-1, 1, n(1)+2); x = x(2:end-1);
f = exp(-10 * x.^2);
b = tt_tensor(f);
for j = 1 : d-1
    b = tkron(b, tt_tensor(f));
end
