function [tt]=tt_gaussian(I, R)
%
%compute a random TT-tensor with normalized independent Gaussian-distributed entries.
%returns a Random Gaussian TT-tensor with mode sizes I(1), ..., I(N) and TT ranks R(1), ..., R(N + 1).
% ---------------------------
r=R(:);
n=I(:);
tt=tt_tensor;
tt.n=n;
d = length(n);
tt.d=d;
tt.r = r;
ps=cumsum([1;n.*r(1:d).*r(2:d+1)]);
tt.ps = ps;
cr = zeros(ps(d+1),1);
for i=1:d
    cr(ps(i):ps(i+1)-1) = randn(ps(i+1)-ps(i),1)/sqrt(ps(i+1)-ps(i));
end
tt.core=cr;
end
