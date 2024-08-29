function S = TT_kronecker_embedding(k, n, d)
% create a Khatri_Rao_embedding in TT format, that is we construct a cell
% with each entry equal to a random vector in Kronecker format

rng(19)


cc = cell(d, 1);
for j = 1:d
    cc{j}(1,:,:,1) = randn(k(j), n(j), 1)/(n(j)*k(j));
end
S = cell2core(tt_matrix,cc);
end
