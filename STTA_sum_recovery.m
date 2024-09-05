function T = STTA_sum_recovery(Psi, Omega, y, d, N, ranks)

Nx = d;
Ix = N;

C = cell(Nx, 1);

n_addend = length(Psi);
Sum_Psi = cell(1, length(Psi{1}));
for j = 1:Nx
    Sum_Psi{j} = y(1)*Psi{1}{j};
end

for i = 2:n_addend
    for j = 1:Nx
        Sum_Psi{j} = Sum_Psi{j} + y(i)*Psi{i}{j};
    end
end

Sum_Omega = cell(1, length(Omega{1}));
for j = 1:Nx-1
    Sum_Omega{j} = y(1)*Omega{1}{j};
end

for i = 2:n_addend
    for j = 1:Nx-1
        Sum_Omega{j} = Sum_Omega{j} + y(i)*Omega{i}{j};
    end
end

for i = 1: Nx-1
    % [U,S,V] = svd(Sum_Omega{i});
    % S(S<5*eps*S(1)) = 0;
    % C{i} = (Sum_Psi{i}*V)*pinv(S)*U';
    %[Qo, Ro] = qr(Sum_Omega{i}, 0);
    %C{i} = Sum_Psi{i}*pinv(Sum_Omega{i});
    %C{i} = (Sum_Psi{i}/Ro)*Qo';
    warning('off', 'MATLAB:rankDeficientMatrix');
    C{i} = Sum_Psi{i} / Sum_Omega{i};
    warning('on', 'MATLAB:rankDeficientMatrix');
end
C{Nx} = Sum_Psi{Nx};

%DA FINIRE
T = tt_tensor;
% dovremmo completare i fields di T
T.d = Nx;
T.n = Ix;

T.core = zeros(1,1);
for i = 1:Nx
    T.r(i) = ranks(i);
end
T.r(Nx+1) = 1;
T.r = T.r';
ps = cumsum([1;Ix.*T.r(1:Nx).*T.r(2:Nx+1)]);
T.ps = ps;
for i = 1:Nx
    T.core(T.ps(i):T.ps(i+1)-1) = C{i}(:)';
end