%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CT Reconstruction Problems 
%
%  F(x) = 1/2||Ax-a||^2;  G = iota_C;  H(x) = mu ||.||_12 (TV); 
%  B is the discret TV using the finite difference scheme with the periodic
%  condition
%
%    Contact:
%       Ming Yan yanm @ math.msu.edu
%       Downloadable from https://github.com/mingyan08/PD3O  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc;
addpath('fcns','data','output') 
N       = 128; % The image will be N x N
M       = N;
MN      = M * N;
nt      = 50;  % The number of projections
image   = phantom(N) * 255;
u       = image;

% Load the given Radon transform matrix R with size 9250*16384
load R_128.mat
A       = R;

% Form the simulated measurements
realb   = A * u(:);
%randns=randn(size(realb)); % Generate a normally distributed random matrix
                            % In order to reuse and compare, we save a
                            % results to randns.mat
load randns.mat 
sigma   = 1; % Variance
a       = realb + sigma * randns;

uFBP    = FBP(reshape(a, numel(a)/nt, nt), N);
uFBP    = uFBP * N/2;
x0      = reshape(uFBP, MN, 1);

% Generate the discrete gradient matrix.
% For large scale problems, we can realize it through a function operation.
%[B1,B2]=generate_B(M,N);
[B1, B2]    = generate_B_circular(M, N);
%B=[B1;B2];

mu      = 0.05;

%----------------------- Set optional items ------------------------

cs={'r--','b--','k--','r-','b-','k-','r-.','b-.','k-.'};

obj         = PrimalDual;  % using the clas PrimalDual
%% Define all its functions:
% Create handles to functions F, G, H, A:
obj.myF = @(x) 0.5 * (A * x - a)' * (A * x - a);
obj.myG = @(x) 0;
obj.myH = @(y) mu * sum(sqrt(y(1:MN).^2 + y(MN+1:end).^2));
obj.myA = @(x) [B1 * x; B2 * x];
% Create handles to adjoint/gradient/prox:
ST          = @(u, t) u .* max(0, 1 - t ./ [abs(sqrt(u(1:MN).^2 + u(MN+1:end).^2)); abs(sqrt(u(1:MN).^2 + u(MN+1:end).^2))]);
obj.myGradF = @(x) A' * (A * x - a);
obj.myProxG = @(x, t) min(max(x, 0), 255);
obj.myProxH = @(y, t) ST(y, t * mu);
obj.myAdjA  = @(y) [B1' B2'] * y;
beta    = 1/normest(A*A');

%% Define the parameters 
obj.gamma   = 1*beta;
obj.lambda  = 1/16;
obj.input.iter  = 1000;
obj.input.x     = x0;
obj.input.s     = zeros(2*MN,1);

E_min = 1.300790966333291e+04;

h1  = figure;
j   = 1;
tic
[x_PD3O, s_PD3O, E_PD3O]  = obj.minimize('PD3O', 0);
time(j)=toc;

figure(h1)
semilogy(E_PD3O(:)./E_min-1,cs{j});

j=2;
tic
[x_PDFP, s_PDFP, E_PDFP]  = obj.minimize('PDFP', 0);
time(j)=toc;

figure(h1)
hold on
semilogy(E_PDFP(:)./E_min-1,cs{j});

j=3;
tic
[x_CV, s_CV, E_CV]        = obj.minimize('CV', 0);
time(j)=toc;

figure(h1)
hold on
semilogy(E_CV(:)./E_min-1,cs{j});


j=4;
obj.gamma   = 1.5 * beta;

tic
[x_PD3O2, s_PD3O2, E_PD3O2]  = obj.minimize('PD3O', 0);
time(j)=toc;

figure(h1)
hold on
semilogy(E_PD3O2(:)./E_min-1,cs{j});

j=5;
tic
[x_PDFP2, s_PDFP2, E_PDFP2]  = obj.minimize('PDFP', 0);
time(j)=toc;

figure(h1)
hold on
semilogy(E_PDFP2(:)./E_min-1,cs{j});

j=6;
tic
[x_CV2, s_CV2, E_CV2]        = obj.minimize('CV', 0);
time(j)=toc;

figure(h1)
hold on
semilogy(E_CV2(:)./E_min-1,cs{j});


j=7;
obj.gamma   = 1.9 * beta;

tic
[x_PD3O3, s_PD3O3, E_PD3O3]  = obj.minimize('PD3O', 0);
time(j)=toc;

figure(h1)
hold on
semilogy(E_PD3O3(:)./E_min-1,cs{j});

j=8;
tic
[x_PDFP3, s_PDFP3, E_PDFP3]  = obj.minimize('PDFP', 0);
time(j)=toc;

figure(h1)
hold on
semilogy(E_PDFP3(:)./E_min-1,cs{j});

j=9;
tic
[x_CV3, s_CV3, E_CV3]        = obj.minimize('CV', 0);
time(j)=toc;

figure(h1)
hold on
semilogy(E_CV3(:)./E_min-1,cs{j});

h_legend = legend({'PD3O-$\gamma_1$','PDFP-$\gamma_1$','Condat-Vu-$\gamma_1$','PD3O-$\gamma_2$','PDFP-$\gamma_2$','Condat-Vu-$\gamma_2$','PD3O-$\gamma_3$','PDFP-$\gamma_3$','Condat-Vu-$\gamma_3$'},'Interpreter','latex');
set(h_legend,'FontSize',10);
xlabel('iteration','FontSize',20)
ylabel('$\frac{f-f^*}{f^*}$','Interpreter','LaTex','FontSize',20);

myprint('output/CT_energy_1',h1)