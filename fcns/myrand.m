function x = myrand(n, ratio, maxlength)
% First generate randomly the positions of first non-zeros elements in x
% with probability ratio; 
% Second generate the number of non-zeros in succession by using
% the  uniformly distributed in the interval (0,maxlength);
% Third appoint their values as a same value belong to [-3, -2, -1, 0, 1, 2, 3, 4, -4]
% by turn at last. 
x       = zeros(n,1);
rand('state', 3);
a       = rand(n,1);
pos     = find(a < ratio);
m       = length(pos);
rand('state',2);
b       = ceil(rand(m,1) * maxlength);
for i = 1:m
    r   = mod(i,9) - 4;
    x(pos(i):pos(i) + b(i)) = r;
end
