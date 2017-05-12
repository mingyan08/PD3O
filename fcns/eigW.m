function [lambda2, lambdan] = eigW(W)
% analyze the eigenvalues of W
a = eig(W);
[~, indmax]=max(a);
ind = ones(size(a)); ind(indmax) = 0;
lambda2 = max(a(logical(ind)));
lambdan = min(a);
end