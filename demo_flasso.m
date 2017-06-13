%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fused LASSO penalized problems 
%
%   min    1/2 || A x - a||^2 + mu1 sum_i |x_{i+1}-x_i|+mu2 ||x||_1.   (1)
% x in R^n 
%
%  F(x) = 1/2||Ax-a||^2;  G(x) = mu2 ||.||_1;  H(x) = mu1 ||.||_1; 
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
close all
clear
clc;
addpath('fcns','data','output') 
r   = 500;  
n   = 10000;    % The data matrix is of size r x n
% ---------------------- Generate random data ----------------------
% for reproducibility
randNum     = 1;
%randn('state', (randNum - 1) * 3 + 1);  % old version
rng((randNum - 1) * 3 + 1, 'v5normal');
A           = randn(r, n);                        % the data matrix
x_true      = myrand(n, 0.005, 50);               % sparse varibles values
%randn('state', (randNum - 1) * 3 + 3);  % old version
rng((randNum - 1) * 3 + 3, 'v5normal');
noise       = randn(r,1);                       % noise       
a           = A * x_true + 0.01 * noise;        % the response
nOnes       = ones(n, 1);
B           = diag(-nOnes, 0) +  diag(nOnes(1:end - 1), 1);
B(end,:)    = [];
B           = sparse(B);

%----------------------- Set parameters for the problem ------------------------
beta        = 1/normest(A*A');
mu1         = 200;
mu2         = 20;
obj         = PrimalDual;  % using the clas PrimalDual

%% Define all its functions:
% Create handles to functions F, G, H, A, and H*:
obj.myF = @(x) 0.5 * (A * x - a)' * (A * x - a);
obj.myG = @(x) mu2 * sum(abs(x));
obj.myH = @(y) mu1 * sum(abs(y));
obj.myHA = @(y) 10000*(max(abs(y)) > mu1+1e-5);  % the conjugate function of H
obj.myA = @(x) B * x;
% Create handles to adjoint/gradient/prox:
ST          = @(u, t) u .* max(0, 1 - t ./ abs(u));
obj.myGradF = @(x) A' * (A * x - a);
obj.myProxG = @(x, t) ST(x, t * mu2);
obj.myProxH = @(y, t) ST(y, t * mu1);
obj.myAdjA  = @(y) B' * y;

%% Define the parameters for the algorithm
obj.gamma   = 0.5 * beta;          % two parameters for the primal-dual algorithms
obj.lambda  = 3/16;             % we will choose different gammas in this example

load x_PD3O_flasso.mat           % the true solution that is computed from 10,000 iterations using PD3O
E_min   = min(E_PD3O_min);       % computed from 10,000 iterations using PD3O
obj.myF2 = @(x) (A * x_PD3O_min - a)'* (A * (x - x_PD3O_min)); % linearization at x_PD3O_min (the optimal solution)
iter    = 4000;                  % the number of iterations

obj.input.x     = zeros(n,1);
obj.input.s     = zeros(n-1,1);
obj.input.x_min = x_PD3O_min;
obj.input.s_min = s_PD3O_min;
obj.input.iter  = iter;


%% Run the primal-dual codes
% PD3O %%%
j   = 1;
tic
[x_PD3O, s_PD3O, E_PD3O, out_PD3O]  = obj.minimize('PD3O', 1);
time(j) = toc;

% PDFP %%%
j   = 2;
tic
[x_PDFP, s_PDFP, E_PDFP, out_PDFP]  = obj.minimize('PDFP', 1);
time(j) = toc;

% CV %%%
j   = 3;
tic
[x_CV, s_CV, E_CV, out_CV]        = obj.minimize('CV', 1);
time(j) = toc;

% AFBA %%%
j   = 4;
tic
[x_AFBA, s_AFBA, E_AFBA, out_AFBA]  = obj.minimize('AFBA', 1);
time(4) = toc;

%% choose a different gamma 
obj.gamma   = 1 * beta;
% PD3O %%%
j   = 11;
tic
[x_PD3O2, s_PD3O2, E_PD3O2, out_PD3O2]  = obj.minimize('PD3O', 1);
time(j) = toc;

% PDFP %%%
j   = 12;
tic
[x_PDFP2, s_PDFP2, E_PDFP2, out_PDFP2]  = obj.minimize('PDFP', 1);
time(j) = toc;

% CV %%%
j   = 13;
tic
[x_CV2, s_CV2, E_CV2, out_CV2]        = obj.minimize('CV', 1);
time(j) = toc;

% AFBA %%%
j   = 14;
tic
[x_AFBA2, s_AFBA2, E_AFBA2, out_AFBA2]  = obj.minimize('AFBA', 1);
time(4) = toc;

%% choose another gamma 
obj.gamma   = 2.0 * beta;
% PD3O %%%
j   = 21;
tic
[x_PD3O3, s_PD3O3, E_PD3O3, out_PD3O3]  = obj.minimize('PD3O', 1);
time(j) = toc;

% PDFP %%%
j   = 22;
tic
[x_PDFP3, s_PDFP3, E_PDFP3, out_PDFP3]  = obj.minimize('PDFP', 1);
time(j) = toc;

% CV %%%
j   = 23;
tic
[x_CV3, s_CV3, E_CV3, out_CV3]        = obj.minimize('CV', 1);
time(j) = toc;

% AFBA %%%
j   = 24;
tic
[x_AFBA3, s_AFBA3, E_AFBA3, out_AFBA3]  = obj.minimize('AFBA', 1);
time(4) = toc;

%% Run the primal-dual codes with different settings
obj.gamma   = 1 * beta;          % two parameters for the primal-dual algorithms
obj.lambda  = 1/8;             % we will choose different lambda in this example

% PD3O %%%
j   = 1;
tic
[x_PD3O4, s_PD3O4, E_PD3O4, out_PD3O4]  = obj.minimize('PD3O', 1);
time2(j) = toc;

% PDFP %%%
j   = 2;
tic
[x_PDFP4, s_PDFP4, E_PDFP4, out_PDFP4]  = obj.minimize('PDFP', 1);
time2(j) = toc;

% CV %%%
j   = 3;
tic
[x_CV4, s_CV4, E_CV4, out_CV4]        = obj.minimize('CV', 1);
time2(j) = toc;

% CV %%%
j   = 4;
tic
[x_AFBA4, s_AFBA4, E_AFBA4, out_AFBA4]  = obj.minimize('AFBA', 1);
time2(j) = toc;

%% choose a different gamma 
obj.lambda   = 1.5/8;
% PD3O %%%
j   = 11;
tic
[x_PD3O5, s_PD3O5, E_PD3O5, out_PD3O5]  = obj.minimize('PD3O', 1);
time2(j) = toc;

% PDFP %%%
j   = 12;
tic
[x_PDFP5, s_PDFP5, E_PDFP5, out_PDFP5]  = obj.minimize('PDFP', 1);
time2(j) = toc;

% CV %%%
j   = 13;
tic
[x_CV5, s_CV5, E_CV5, out_CV5]        = obj.minimize('CV', 1);
time2(j) = toc;

% AFBA %%%
j   = 14;
tic
[x_AFBA5, s_AFBA5, E_AFBA5, out_AFBA5]  = obj.minimize('AFBA', 1);
time2(4) = toc;

%% choose another gamma 
obj.lambda   = 1/4;
% PD3O %%%
j   = 21;
tic
[x_PD3O6, s_PD3O6, E_PD3O6, out_PD3O6]  = obj.minimize('PD3O', 1);
time2(j) = toc;

% PDFP %%%
j   = 22;
tic
[x_PDFP6, s_PDFP6, E_PDFP6, out_PDFP6]  = obj.minimize('PDFP', 1);
time2(j) = toc;

% CV %%%
j   = 23;
tic
[x_CV6, s_CV6, E_CV6, out_CV6]        = obj.minimize('CV', 1);
time2(j) = toc;

% CV %%%
j   = 24;
tic
[x_AFBA6, s_AFBA6, E_AFBA6, out_AFBA6]  = obj.minimize('AFBA', 1);
time2(j) = toc;

%% save the final results for plot
save output_flasso.mat out_* E_* x_* s_*