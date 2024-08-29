function Sw = KR_sketching(w, S)

d = length(S);
V = cell(1, d);
k = size(S{1}, 1);

for j = 1 : d
    TT = reshape(w.core(w.ps(j) : w.ps(j+1)-1), w.r(j), w.n(j), w.r(j+1));
    TT = permute(TT, [2, 1, 3]); 
    TT = reshape(TT, w.n(j), w.r(j) * w.r(j+1));
    V{j} = S{j} * TT;
    V{j} = reshape(V{j}, k, w.r(j), w.r(j+1));
    V{j} = permute(V{j}, [2, 3, 1]);
end

Sw = zeros(k, 1);

for i = 1 : k
    VV = V{1};
    SWW = squeeze(VV(:,:,i));
    for j = 2 : d
        VV = V{j};
        SWW = SWW * VV(:,:,i);
    end
    Sw(i) = SWW;
end

end