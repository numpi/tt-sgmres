function [Omega, Psi, varargout] = STTA_contractions(X, Y, Z, Ycell, Zcell) % compute left and right sketches/contractions of T.

% compute:
% L_k = B_{\leq k}^T X_{\leq k}
% R_k = X_{\geq k} A_{\geq k}^T
% then form Psi and Omega 
% as described in STTA.

% if a third output is given it computes matrices representing the core 
% tensors

% if a fourth output is given it computes the approximated tensor in
% tt-format

Nx = X.d;
Ix = X.n;
Rx = X.r;

%Ny = Y.d;
Iy = Y.n;
Ry = Y.r;

Nz = Z.d;
Iz = Z.n;
Rz = Z.r;


L = cell(Nx-1, 1);

%L{1} =reshape(Y.core(Y.ps(1):Y.ps(2)-1), Iy(1), Ry(2))'* reshape(X.core(X.ps(1):X.ps(2)-1), Ix(1), Rx(2));
L{1} =Ycell{1}'* reshape(X.core(X.ps(1):X.ps(2)-1), Ix(1), Rx(2));
for n = 2:Nx-1
    %proviamo = reshape(Y.core(Y.ps(n):Y.ps(n+1)-1), Ry(n), Iy(n)*Ry(n+1));
    tempLY = L{n-1}'*Ycell{n};%*proviamo;
    reshaped_tempLY = reshape(tempLY,Rx(n)*Ix(n),Ry(n+1));
    L{n} = reshaped_tempLY'*reshape(X.core(X.ps(n):X.ps(n+1)-1), Rx(n)*Ix(n), Rx(n+1));
end


R = cell(Nx-1, 1);

%R{Nx-1} = reshape(X.core(X.ps(Nx):X.ps(Nx+1)-1), Rx(Nx), Ix(Nx)) * reshape(Z.core(Z.ps(Nx):Z.ps(Nx+1)-1), Rz(Nz), Iz(Nz))';
R{Nx-1} = reshape(X.core(X.ps(Nx):X.ps(Nx+1)-1), Rx(Nx), Ix(Nx)) * Zcell{Nz}';
for n = Nx-1:-1:2
    %tempZR = reshape(Z.core(Z.ps(n):Z.ps(n+1)-1), Iz(n)*Rz(n), Rz(n+1))*R{n}';
    tempZR = Zcell{n}*R{n}';
    reshaped_tempZR = reshape(tempZR,Rz(n),Ix(n)*Rx(n+1));
    R{n-1} = reshape(X.core(X.ps(n):X.ps(n+1)-1), Rx(n), Ix(n)*Rx(n+1))*reshaped_tempZR';   
end


Omega = cell(Nx-1,1);

for i = 1:Nx-1
    Omega{i} = L{i}*R{i};
end

Psi = cell(Nx-1,1);

Psi{1} = reshape(X.core(X.ps(1):X.ps(2)-1), Ix(1), Rx(2))*R{1};
for n = 2:Nx-1
    tempLX = L{n-1}*reshape(X.core(X.ps(n):X.ps(n+1)-1), Rx(n), Ix(n)*Rx(n+1));
    Psi{n} = reshape(tempLX, Ry(n)*Ix(n),Rx(n+1))*R{n};
end
Psi{Nx} = L{Nx-1}*reshape(X.core(X.ps(Nx):X.ps(Nx+1)-1), Rx(Nx), Ix(Nx));

varargout = {};

if nargout >= 3

    C = cell(Nx, 1);
    for i = 1: Nx-1
        C{i} = Psi{i}*pinv(Omega{i});
    end
    C{Nx} = Psi{Nx};
    varargout{1} = C;

    if nargout == 4
        T = tt_tensor;
        
        T.d = Nx;
        T.n = Ix;
        
        T.core = zeros(1,1);
        for i = 1:Nx
            T.r(i) = Y.r(i);
        end
        T.r(Nx+1) = 1;
        T.r = T.r';
        ps = cumsum([1;Ix.*T.r(1:Nx).*T.r(2:Nx+1)]);
        T.ps = ps;
        for i = 1:Nx
            T.core(T.ps(i):T.ps(i+1)-1) = C{i}(:)';
        end
        varargout{2} = T;
    end

end