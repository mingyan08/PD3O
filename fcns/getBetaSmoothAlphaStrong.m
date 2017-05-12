function [Lips,mus] = getBetaSmoothAlphaStrong
% 170410 get the belta-smooth an alpha-strong
global n m M
Lips = zeros(n,1);
mus = zeros(n,1);

for i = 1:n
    a = eig(M(:,:,i)'*M(:,:,i));
    % belta
    Lips(i) = max(a);
    % alpha
    mus(i) = min(a);
end
end