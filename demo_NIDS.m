%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Decentralized Consensus Optimization Problem (NIDS)
%
%
%  F(x) = 1/2||Ax-a||^2;  G(x) = mu2 ||.||_1;  H(x) = iota_0; 
%  A = (I-W)^{1/2}
%
%    Contact:
%       Ming Yan yanm @ math.msu.edu
%       Downloadable from https://github.com/mingyan08/PD3O  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demo_NIDS
    close all
    clear
    clc;

    global n m p M y_ori lam
    addpath('fcns','data','output') 

    n = 40;
    m = 3;
    p = 200;

    L = n;
    per = 4/L;

    % may changed in the following function
    min_mu = 0.1; % set the smallest strongly convex parameter mu in S
    max_Lips = 1; % set the Lipschitz constant

    % lam is the parameter in function R
    % [M, x_ori, y_ori, lam, W] = generateAll(m, p, n, per,...
    %     'withNonsmoothR', min_mu,max_Lips);
    W = generateW(L, per);
    [M, x_ori, y_ori, lam] = generateS(m, p, n,...
        'withNonsmoothR', min_mu, max_Lips);

    [~, lambdan] = eigW(W); % find the smallest eigenvalue of W
    lambda_max = 1/(1-lambdan); 
    % set parameters
    x0      = zeros(n,p);
    x_star  = x_ori;     % true solution
    x_star_norm = norm(x_star, 'fro');
    % Set the parameter for the solver
    W2      = eye(n) - W;
    [U,S,V] = svd(W2);
    W       = U * sqrt(S) * V';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start using the PrimalDual class
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj            =  PrimalDual;  % using the class PrimalDual

    %% Define all its functions:
    % Create handles to functions F, G, H, A:
    obj.myF        = @(x) feval(@funS, x);
    obj.myG        = @(x) feval(@funR, x);
    obj.myA        = @(x) W * x;
    % Create handles to adjoint/gradient/prox:
    obj.myGradF    = @(x) feval(@funGradS, x);
    obj.myProxG    = @(x,t) feval(@funProxR, x, t);
    obj.myProxH    = @(y, t) zeros(size(y));
    obj.myAdjA     = @(y) W' * y;

    %% Define the parameters 
    obj.gamma   = 1;                % two parameters for the primal-dual algorithms
    obj.lambda  = 0.5 * lambda_max;              % we will choose different gammas in this example

    cs={'r--','b--','k--','r-','b-','k-','r-.','b-.','k-.'};

    %% Define the parameters 
    iter    = 10000;                  % the number of iterations

    %% Run the primal-dual codes
    h1  = figure;
    legend_text = {};
    % PD3O %%%
    j   = 1;
    tic
    [x_PD3O, s_PD3O, E_PD3O, out_PD3O]  = obj.minimize(x0, x0, iter, 'PD3O', x_star, 1);
    time(j) = toc;

    figure(h1)
    semilogy(out_PD3O.LS(:)/x_star_norm,cs{j});
    ylim([1e-4,1])
    legend_text(end+1) = {'PD3O-$\gamma_1$'};
    % PDFP %%%
    j   = 2;
    tic
    [x_PDFP, s_PDFP, E_PDFP, out_PDFP]  = obj.minimize(x0, x0, iter, 'PDFP', x_star, 1);
    time(j) = toc;

    figure(h1)
    hold on
    semilogy(out_PDFP.LS(:)/x_star_norm,cs{j});
    legend_text(end+1) = {'PDFP-$\gamma_1$'};
    % CV %%%
    j   = 3;
    tic
    [x_CV, s_CV, E_CV, out_CV]        = obj.minimize(x0, x0, iter, 'CV', x_star, 1);
    time(j) = toc;

    figure(h1)
    hold on
    semilogy(out_CV.LS(:)/x_star_norm, cs{j});
    legend_text(end+1) = {'Condat-Vu-$\gamma_1$'};
    %% choose a different gamma 
    obj.gamma   = 1.5;
    % PD3O %%%
    j   = 4;
    tic
    [x_PD3O2, s_PD3O2, E_PD3O2, out_PD3O2]  = obj.minimize(x0, x0, iter, 'PD3O', x_star, 1);
    time(j) = toc;

    figure(h1)
    hold on
    semilogy(out_PD3O2.LS(:)/x_star_norm, cs{j});
    legend_text(end+1) = {'PD3O-$\gamma_2$'};
    % PDFP %%%
    j   = 5;
    tic
    [x_PDFP2, s_PDFP2, E_PDFP2, out_PDFP2]  = obj.minimize(x0, x0, iter, 'PDFP', x_star, 1);
    time(j) = toc;

    figure(h1)
    hold on
    semilogy(out_PDFP2.LS(:)/x_star_norm,cs{j});
    legend_text(end+1) = {'PDFP-$\gamma_2$'};

    %% choose another gamma 
    obj.gamma   = 2.0;
    % PD3O %%%
    j   = 7;
    tic
    [x_PD3O3, s_PD3O3, E_PD3O3, out_PD3O3]  = obj.minimize(x0, x0, iter, 'PD3O', x_star, 1);
    time(j) = toc;

    figure(h1)
    hold on
    semilogy(out_PD3O3.LS(:)/x_star_norm,cs{j});
    legend_text(end+1) = {'PD3O-$\gamma_3$'};
    % PDFP %%%
    j   = 8;
    tic
    [x_PDFP3, s_PDFP3, E_PDFP3, out_PDFP3]  = obj.minimize(x0, x0, iter, 'PDFP', x_star, 1);
    time(j) = toc;

    figure(h1)
    hold on
    semilogy(out_PDFP3.LS(:)/x_star_norm,cs{j});
    legend_text(end+1) = {'PDFP-$\gamma_3$'};

    h_legend = legend(legend_text,'Interpreter','latex');
    set(h_legend,'FontSize',10);
    xlabel('iteration','FontSize',20)
    ylabel('$\frac{\|x-x^*\|}{\|x^*\|}$','Interpreter','LaTex','FontSize',20);

    myprint('output/NIDS_1',h1)   % print the file in .pdf and .eps

    %% Run the primal-dual codes with different settings
    obj.gamma   = 1.9;          % two parameters for the primal-dual algorithms
    obj.lambda  = 0.05 * lambda_max;             % we will choose different lambda in this example

    h2  = figure;
    legend_text = {};
    % PD3O %%%
    j   = 1;
    tic
    [x_PD3O, s_PD3O, E_PD3O, out_PD3O]  = obj.minimize(x0, x0, iter, 'PD3O', x_star, 1);
    time2(j) = toc;

    figure(h2)
    semilogy(out_PD3O.LS(:)/x_star_norm,cs{j});
    ylim([1e-4,1])
    legend_text(end+1) = {'PD3O-$\gamma_1$'};
    % PDFP %%%
    j   = 2;
    tic
    [x_PDFP, s_PDFP, E_PDFP, out_PDFP]  = obj.minimize(x0, x0, iter, 'PDFP', x_star, 1);
    time2(j) = toc;

    figure(h2)
    hold on
    semilogy(out_PDFP.LS(:)/x_star_norm, cs{j});
    legend_text(end+1) = {'PDFP-$\gamma_1$'};
    % CV %%%
    j   = 3;
    tic
    [x_CV, s_CV, E_CV, out_CV]        = obj.minimize(x0, x0, iter, 'CV', x_star, 1);
    time2(j) = toc;

    figure(h2)
    hold on
    semilogy(out_CV.LS(:)/x_star_norm, cs{j});
    legend_text(end+1) = {'Condat-Vu-$\gamma_1$'};
    %% choose a different gamma 
    obj.lambda   = 0.5 * lambda_max;
    % PD3O %%%
    j   = 4;
    tic
    [x_PD3O2, s_PD3O2, E_PD3O2, out_PD3O2]  = obj.minimize(x0, x0, iter, 'PD3O', x_star, 1);
    time2(j) = toc;

    figure(h2)
    hold on
    semilogy(out_PD3O2.LS(:)/x_star_norm,cs{j});
    legend_text(end+1) = {'PD3O-$\gamma_2$'};
    % PDFP %%%
    j   = 5;
    tic
    [x_PDFP2, s_PDFP2, E_PDFP2, out_PDFP2]  = obj.minimize(x0, x0, iter, 'PDFP', x_star, 1);
    time2(j) = toc;

    figure(h2)
    hold on
    semilogy(out_PDFP2.LS(:)/x_star_norm,cs{j});
    legend_text(end+1) = {'PDFP-$\gamma_2$'};

    %% choose another gamma 
    obj.lambda   = 1.0 * lambda_max;
    % PD3O %%%
    j   = 7;
    tic
    [x_PD3O3, s_PD3O3, E_PD3O3, out_PD3O3]  = obj.minimize(x0, x0, iter, 'PD3O', x_star, 1);
    time2(j) = toc;

    figure(h2)
    hold on
    semilogy(out_PD3O3.LS(:)/x_star_norm,cs{j});
    legend_text(end+1) = {'PD3O-$\gamma_3$'};
    % PDFP %%%
    j   = 8;
    tic
    [x_PDFP3, s_PDFP3, E_PDFP3, out_PDFP3]  = obj.minimize(x0, x0, iter, 'PDFP', x_star, 1);
    time2(j) = toc;

    figure(h2)
    hold on
    semilogy(out_PDFP3.LS(:)/x_star_norm,cs{j});
    legend_text(end+1) = {'PDFP-$\gamma_3$'};


    h_legend = legend(legend_text,'Interpreter','latex');
    set(h_legend,'FontSize',10);
    xlabel('iteration','FontSize',20)
    ylabel('$\frac{\|x-x^*\|}{\|x^*\|}$','Interpreter','LaTex','FontSize',20);

    myprint('output/NIDS_2',h2)   % print the file in .pdf and .eps
end

function a = funGradS(x)
    global n p M y_ori
    a = zeros(n, p);
    for j = 1:n
        a(j,:) = (M(:,:,j)' * (M(:,:,j) * (x(j,:))' - y_ori(:,j)))';
    end
end

function a = funR(x)
    global n lam
    a = 0;
    for j = 1:n
        a = a + lam * norm(x(j,:), 1);
    end
end

function a = funS(x)
    global n M y_ori
    a = 0;
    for j = 1:n
        a   = a + 0.5 * sum((M(:,:,j) * (x(j,:))' - y_ori(:,j)).^2);
    end
end

function a = funProxR(x,t)
    global n p lam
    a = zeros(n, p);
    for j = 1:n
        a(j,:)  = (wthresh(x(j,:), 's', t*lam))';
    end
end