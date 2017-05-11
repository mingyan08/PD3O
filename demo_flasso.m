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
A           = randn(r, n);                       % the data matrix
x_true      = myrand(n, 0.005, 50);               % sparse varibles values
%the sparsity in both x_true and their successive differences of x_true.
%randn('state', (randNum - 1) * 3 + 3);  % old version
rng((randNum - 1) * 3 + 3, 'v5normal');
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

cs={'r--','b--','k--','r-','b-','k-','r-.','b-.','k-.'};

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
obj.gamma   = 1 * beta;          % two parameters for the primal-dual algorithms
obj.lambda  = 0.5/4;             % we will choose different gammas in this example

load x_PD3O_flasso.mat           % the true solution that is computed from 10,000 iterations using PD3O
E_min   = 5.898955369429127e+04; % computed from 10,000 iterations using PD3O
iter    = 1000;                  % the number of iterations

%% Run the primal-dual codes
h1  = figure;
legend_text = {};
% PD3O %%%
j   = 1;
tic
[x_PD3O, s_PD3O, E_PD3O, out]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PD3O', x_PD3O_min, 1);
time(j) = toc;

figure(h1)
semilogy(E_PD3O(:)./E_min-1,cs{j});
ylim([1e-4,1e3])
legend_text(end+1) = {'PD3O-$\gamma_1$'};
% PDFP %%%
j   = 2;
tic
[x_PDFP, s_PDFP, E_PDFP]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PDFP');
time(j) = toc;

figure(h1)
hold on
semilogy(E_PDFP(:)./E_min-1,cs{j});
legend_text(end+1) = {'PDFP-$\gamma_1$'};
% CV %%%
j   = 3;
tic
[x_CV, s_CV, E_CV]        = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'CV');
time(j) = toc;

figure(h1)
hold on
semilogy(E_CV(:)./E_min-1,cs{j});
legend_text(end+1) = {'Condat-Vu-$\gamma_1$'};
%% choose a different gamma 
obj.gamma   = 1.5 * beta;
% PD3O %%%
j   = 4;
tic
[x_PD3O2, s_PD3O2, E_PD3O2]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PD3O');
time(j) = toc;

figure(h1)
hold on
semilogy(E_PD3O2(:)./E_min-1,cs{j});
legend_text(end+1) = {'PD3O-$\gamma_2$'};
% PDFP %%%
j   = 5;
tic
[x_PDFP2, s_PDFP2, E_PDFP2]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PDFP');
time(j) = toc;

figure(h1)
hold on
semilogy(E_PDFP2(:)./E_min-1,cs{j});
legend_text(end+1) = {'PDFP-$\gamma_2$'};
% CV %%%
% j   = 6;
% tic
% [x_CV2, s_CV2, E_CV2]        = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'CV');
% time(j) = toc;
% 
% figure(h1)
% hold on
% semilogy(E_CV2(:)./E_min-1,cs{j});
% legend_text(end+1) = {'Condat-Vu-$\gamma_2$'};

%% choose another gamma 
obj.gamma   = 2.0 * beta;
% PD3O %%%
j   = 7;
tic
[x_PD3O3, s_PD3O3, E_PD3O3]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PD3O');
time(j) = toc;

figure(h1)
hold on
semilogy(E_PD3O3(:)./E_min-1,cs{j});
legend_text(end+1) = {'PD3O-$\gamma_3$'};
% PDFP %%%
j   = 8;
tic
[x_PDFP3, s_PDFP3, E_PDFP3]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PDFP');
time(j) = toc;

figure(h1)
hold on
semilogy(E_PDFP3(:)./E_min-1,cs{j});
legend_text(end+1) = {'PDFP-$\gamma_3$'};
% CV %%%
% j  = 9;
% tic
% [x_CV3, s_CV3, E_CV3]        = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'CV');
% time(j) = toc;
% 
% figure(h1)
% hold on
% semilogy(E_CV3(:)./E_min-1,cs{j});
% legend_text(end+1) = {'Condat-Vu-$\gamma_1$'};

h_legend = legend(legend_text,'Interpreter','latex');
set(h_legend,'FontSize',10);
xlabel('iteration','FontSize',20)
ylabel('$\frac{f-f^*}{f^*}$','Interpreter','LaTex','FontSize',20);

myprint('output/flasso_1',h1)   % print the file in .pdf and .eps

%% Run the primal-dual codes with different settings
obj.gamma   = 1.9 * beta;          % two parameters for the primal-dual algorithms
obj.lambda  = 0.05/4;             % we will choose different lambda in this example

h2  = figure;
legend_text = {};
% PD3O %%%
j   = 1;
tic
[x_PD3O, s_PD3O, E_PD3O, out]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PD3O', x_PD3O_min, 1);
time2(j) = toc;

figure(h2)
semilogy(E_PD3O(:)./E_min-1,cs{j});
ylim([1e-4,1e3])
legend_text(end+1) = {'PD3O-$\gamma_1$'};
% PDFP %%%
j   = 2;
tic
[x_PDFP, s_PDFP, E_PDFP]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PDFP');
time2(j) = toc;

figure(h2)
hold on
semilogy(E_PDFP(:)./E_min-1,cs{j});
legend_text(end+1) = {'PDFP-$\gamma_1$'};
% CV %%%
j   = 3;
tic
[x_CV, s_CV, E_CV]        = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'CV');
time2(j) = toc;

figure(h2)
hold on
semilogy(E_CV(:)./E_min-1,cs{j});
legend_text(end+1) = {'Condat-Vu-$\gamma_1$'};
%% choose a different gamma 
obj.lambda   = 1/4;
% PD3O %%%
j   = 4;
tic
[x_PD3O2, s_PD3O2, E_PD3O2]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PD3O');
time2(j) = toc;

figure(h2)
hold on
semilogy(E_PD3O2(:)./E_min-1,cs{j});
legend_text(end+1) = {'PD3O-$\gamma_2$'};
% PDFP %%%
j   = 5;
tic
[x_PDFP2, s_PDFP2, E_PDFP2]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PDFP');
time2(j) = toc;

figure(h2)
hold on
semilogy(E_PDFP2(:)./E_min-1,cs{j});
legend_text(end+1) = {'PDFP-$\gamma_2$'};
% CV %%%
% j   = 6;
% tic
% [x_CV2, s_CV2, E_CV2]        = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'CV');
% time2(j) = toc;
% 
% figure(h2)
% hold on
% semilogy(E_CV2(:)./E_min-1,cs{j});
% legend_text(end+1) = {'Condat-Vu-$\gamma_2$'};

%% choose another gamma 
obj.lambda   = 1/3;
% PD3O %%%
j   = 7;
tic
[x_PD3O3, s_PD3O3, E_PD3O3]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PD3O');
time2(j) = toc;

figure(h2)
hold on
semilogy(E_PD3O3(:)./E_min-1,cs{j});
legend_text(end+1) = {'PD3O-$\gamma_3$'};
% PDFP %%%
j   = 8;
tic
[x_PDFP3, s_PDFP3, E_PDFP3]  = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'PDFP');
time2(j) = toc;

figure(h2)
hold on
semilogy(E_PDFP3(:)./E_min-1,cs{j});
legend_text(end+1) = {'PDFP-$\gamma_3$'};
% CV %%%
% j  = 9;
% tic
% [x_CV3, s_CV3, E_CV3]        = obj.minimize(zeros(n,1), zeros(n-1,1), iter, 'CV');
% time2(j) = toc;
% 
% figure(h2)
% hold on
% semilogy(E_CV3(:)./E_min-1,cs{j});
% legend_text(end+1) = {'Condat-Vu-$\gamma_1$'};

h_legend = legend(legend_text,'Interpreter','latex');
set(h_legend,'FontSize',10);
xlabel('iteration','FontSize',20)
ylabel('$\frac{f-f^*}{f^*}$','Interpreter','LaTex','FontSize',20);

myprint('output/flasso_2',h2)   % print the file in .pdf and .eps
%% Plot the final results 
figure
plot(x_true, 'g+')
hold on 
plot(x_CV, 'b')
plot(x_PD3O2, 'r')
plot(x_PDFP2, 'y')
hold off
xlabel('');
legend('True','Condat','PD3O','PDFP')

figure
plot(x_true,'g+')
hold on 
plot(x_CV, 'b')
plot(x_PD3O2, 'r')
plot(x_PDFP2, 'y')
hold off
xlabel('');
legend('True','Condat','PD3O','PDFP')
axis([3000,5250,-4,4])