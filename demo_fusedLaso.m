%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fused LASSO penalized problems 
%
%   min    1/2 || A x - a||^2 + mu1 sum_i |x_{i+1}-x_i|+mu2 ||x||_1.   (1)
% x in R^n 
%
%  F(x) = 1/2||Ax-a||^2;  G(x) = ||.||_1;  H(x) = ||.||_1; 
%    (-1  1                   )
%    (   -1  1                )
% B= ( . . . . . . . . . . . .)
%    (                -1  1   )
%    (                   -1  1)
%
%    Contact:
%       Ming Yan yanm @ math.msu.edu
%       Downloadable from https://github.com/mingyan08/PD3O  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc;

r   = 500;  
n   = 10000;    % The data matrix is of size r x n
% ---------------------- Generate random data ----------------------
% for reproducibility
randNum     = 1;
randn('state',(randNum - 1) * 3 + 1);
A           = randn(r, n);                       % the data matrix
x_true      = myrand(n, 0.005, 50);               % sparse varibles values
%the sparsity in both x_true and their successive differences of x_true.
randn('state',(randNum - 1) * 3 + 3);
noise       = randn(r,1);                       % noise       
a           = A * x_true + 0.01 * noise;        % the response
nOnes       = ones(n, 1);
B           = diag(-nOnes, 0) +  diag(nOnes(1:end - 1), 1);
B(end,:)    = [];
B           = sparse(B);

%----------------------- Set optional items ------------------------
beta        = 1/normest(A*A');
mu1         = 200;
mu2         = 20;

cs={'r-', 'g-', 'b-', 'k-', 'm-', 'y-'};

obj         = PrimalDual;  % using the clas PrimalDual
%% Define all its functions:
% Create handles to functions F, G, H, A:
obj.myF = @(x) 0.5 * (A * x - a)' * (A * x - a);
obj.myG = @(x) mu2 * sum(abs(x));
obj.myH = @(y) mu1 * sum(abs(y));
obj.myA = @(x) B * x;
% Create handles to adjoint/gradient/prox:
ST          = @(u, t) u .* max(0, 1 - t ./ abs(u));
obj.myGradF = @(x) A' * (A * x - a);
obj.myProxG = @(x, t) ST(x, t * mu2);
obj.myProxH = @(y, t) ST(y, t * mu1);
obj.myAdjA  = @(y) B' * y;

%% Define the parameters 
obj.gamma   = 1.9 * beta;
obj.lambda  = 0.19/4;

h1  = figure;
h2  = figure;

j = 1;
tic
[x_PD3O, s_PD3O, E_PD3O]  = obj.minimize(zeros(n,1), zeros(n-1,1), 1500, 'PD3O');
time(j) = toc;

figure(h1)
hold on
plot(E_PD3O(:), cs{j});

j = 2;
tic
[x_PDFP, s_PDFP, E_PDFP]  = obj.minimize(zeros(n,1), zeros(n-1,1), 1500, 'PDFP');
time(j) = toc;

figure(h1)
hold on
plot(E_PDFP(:), cs{j});

j = 3;
tic
[x_CV, s_CV, E_CV]        = obj.minimize(zeros(n,1), zeros(n-1,1), 1500, 'CV');
time(j) = toc;

figure(h1)
hold on
plot(E_CV(:), cs{j});

j = 4;
obj.gamma   = 1.9 * beta;
obj.lambda  = 1/4;

tic
[x_PD3O2, s_PD3O2, E_PD3O2]  = obj.minimize(zeros(n,1), zeros(n-1,1), 1500, 'PD3O');
time(j) = toc;

figure(h1)
hold on
plot(E_PD3O2(:), cs{j});

j = 5;
tic
[x_PDFP2, s_PDFP2, E_PDFP2]  = obj.minimize(zeros(n,1), zeros(n-1,1), 1500, 'PDFP');
time(j) = toc;

figure(h1)
hold on
plot(E_PDFP2(:), cs{j});

j = 6;
tic
[x_CV2, s_CV2, E_CV2]        = obj.minimize(zeros(n,1), zeros(n-1,1), 1500, 'CV');
time(6) = toc;

figure(h1)
hold on
plot(E_CV2(:), cs{j});
legend('PD3O1', 'PDFP1', 'Condat1', 'PD3O2', 'PDFP2', 'Condat2')

figure(h2)
subplot(2, 2, 1);
plot(x_true, 'g+')
hold on 
plot(x_CV, 'b')
plot(x_PD3O, 'r')
plot(x_PDFP, 'y')
hold off
xlabel('');
legend('True','Condat','PD3O','PDFP')

subplot(2, 2, 2)
plot(x_true, 'g+')
hold on 
plot(x_CV, 'b')
plot(x_PD3O, 'r')
plot(x_PDFP, 'y')
hold off
xlabel('');
legend('True', 'Condat', 'PD3O', 'PDFP')
axis([3000, 5250, -4, 4])

subplot(2, 2, 3);
plot(x_true, 'g+')
hold on 
plot(x_CV2, 'b')
plot(x_PD3O2, 'r')
plot(x_PDFP2, 'y')
hold off
xlabel('');
legend('True','Condat','PD3O','PDFP')

subplot(2, 2, 4)
plot(x_true, 'g+')
hold on 
plot(x_CV2, 'b')
plot(x_PD3O2, 'r')
plot(x_PDFP2, 'y')
hold off
xlabel('');
legend('True', 'Condat', 'PD3O', 'PDFP')
axis([3000, 5250, -4, 4])