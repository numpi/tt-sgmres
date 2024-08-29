function S = TT_khatri_rao_embedding(k, n, d)
%TT_KHATRI_RAO_EMBEDDING2 

S = cell(1, d);
for j = 1 : d
    S{j} = randn(k, n(j)) / (sqrt(k))^(1/d);
end

end

