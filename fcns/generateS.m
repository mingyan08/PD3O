function [M_f,x_f,y_f,lam] = generateS(m,p,n,...
    method,mu,Lips,Mrate,modifyM_num)

if (nargin < 6)
    Lips = 1;
end
if (nargin < 7)
    Mrate = 1;
end
if (nargin < 8)
    modifyM_num = 4;
end
L_f=1;


% eg_num=1;
M_mean=0;
M_var=1;

max_iter = 20000*2;

rng(8,'twister');
% x_F in [-15,15], only around 1/4 elements are non-zeros
max_x_F = 15; min_x_F = -15;
x_F = ((max_x_F-min_x_F).*rand(p,1) + min_x_F) .* (rand(p,1)<0.25);

% size(M_f)=(m,p,n);
M_f=random('Normal',M_mean,M_var,m,p,n);

for i=1:n
    M_f(:,:,i)=((sqrt(L_f)/norm(M_f(:,:,i),2))*M_f(:,:,i)')';
end;

switch modifyM_num
    case 4
        [M_f, M_f_T] = modifyM4(m,p,n,M_f,mu,Lips,Mrate);
    case 7
        [M_f, M_f_T] = modifyM7(m,p,n,M_f,mu,Lips,Mrate);
    otherwise
        disp('other value')
end

M_F=(reshape(M_f_T,p,n*m))';

y_F = M_F*x_F + (max_x_F-min_x_F)/10000*rand(n*m,1);
y_f=reshape(y_F,m,n);

switch method
    case {'NIDS','PG-EXTRA','PG-EXTRAF',...
            'PGEXTRAF', 'PGEXTRA', 'withNonsmoothR'}
        lam = 1/n;
        oneOverLips=1/norm(M_F,2)^2;
        % fista
        t = 1;  x = x_F; y = x;
        %       GradG(y) = M_F'*(M_F*y-y_F)
        %       xnew = ProxF( y - 1/L*GradG(y), 1/L );
        %       1/L = oneOverLips
        for i=1:max_iter
%           the 3rd input: lam * n * oneOverLips = oneOverLips;
            xnew = wthresh(y-oneOverLips*M_F'*(M_F*y-y_F),'s',oneOverLips);
            tnew = (1+sqrt(1+4*t^2))/2;
            y = xnew + (t-1)/(tnew)*(xnew-x);
            x = xnew; t = tnew;
        end
        x_CS = x;

    case {'NIDSS','EXTRAF', 'EXTRA','withoutNonsmoothR'}
        lam = 0;
        x_CS = M_F\y_F;
end

% lambda_vec=lambda_ori/n;
x_f=(x_CS*ones(1,n))';

end

function [M_f, M_f_T] = modifyM4(m,p,n,M_f,mu,Lips,Mrate)
% change based on even or odd
if (nargin < 6)
    Lips = 1;
end

if (Lips < mu)
    error(message('the strongly convex parameter mu must < the Lipschitz constant'));
end

min_singular = sqrt(mu);
max_singular = sqrt(Lips);
M_f_T = zeros(p,m,n);
% fix random seed
rng(8,'twister');
for i=1:n
    [ss,vv,dd] = svd(M_f(:,:,i));
    d = diag(vv);
    % change min(singular)
    num_d = length(d);
    d = (max_singular-min_singular).*rand((num_d-2),1) + min_singular;
    d = [max_singular; d; min_singular];
    
    % only change even
    if mod(i,2) == 0
        M_f(:,:,i) = ss(:,1:num_d)*(sqrt(Mrate)*diag(d))*dd(:,1:num_d)';
    else
        M_f(:,:,i) = ss(:,1:num_d)*diag(d)*dd(:,1:num_d)';
    end
    
    M_f_T(:,:,i) = M_f(:,:,i)';
end
end

function [M_f, M_f_T] = modifyM7(m,p,n,M_f,mu,Lips,Mrate)
% change based on even or odd
if (nargin < 6)
    Lips = 1;
end

if (Lips < mu)
    error(message('the strongly convex parameter mu must < the Lipschitz constant'));
end

min_singular = sqrt(mu);
max_singular = sqrt(Lips);
M_f_T = zeros(p,m,n);
% fix random seed
rng(8,'twister');
for i=1:n
    [ss,vv,dd] = svd(M_f(:,:,i));
    d = diag(vv);
    % change min(singular)
    num_d = length(d);
    d = (max_singular-min_singular).*rand((num_d-2),1) + min_singular;
    d = [max_singular; d; min_singular];
    
    % only change even
    if mod(i,16) == 0
        M_f(:,:,i) = ss(:,1:num_d)*(sqrt(Mrate)*diag(d))*dd(:,1:num_d)';
    else
        M_f(:,:,i) = ss(:,1:num_d)*diag(d)*dd(:,1:num_d)';
    end
    
    M_f_T(:,:,i) = M_f(:,:,i)';
end
end





