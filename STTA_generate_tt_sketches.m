function [Y, Z, Ycell, Zcell] = STTA_generate_tt_sketches(N, ry, rz)

    d = length(N);

    Y = tt_gaussian(N, ry);
    Z = tt_gaussian(N, rz);
    
    Ny = Y.d;
    Iy = Y.n;
    Ry = Y.r;
    
    Nz = Z.d;
    Iz = Z.n;
    Rz = Z.r;
        
    Ycell = cell(1, d-1);
    Ycell{1} = reshape(Y.core(Y.ps(1):Y.ps(2)-1), Iy(1), Ry(2));
    for n = 2:d-1
        Ycell{n} = reshape(Y.core(Y.ps(n):Y.ps(n+1)-1), Ry(n), Iy(n)*Ry(n+1));
    end
    
    Zcell = cell(1, d-1);
    Zcell{Nz} = reshape(Z.core(Z.ps(Nz):Z.ps(Nz+1)-1), Rz(Nz), Iz(Nz)); 
    for n = d-1:-1:2
        Zcell{n} = reshape(Z.core(Z.ps(n):Z.ps(n+1)-1), Iz(n)*Rz(n), Rz(n+1));
    end
end